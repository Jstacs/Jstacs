/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.scoringFunctions;

import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * This scoring function implements a cyclic Markov model of arbitrary order and periodicity for any sequence length.
 * The scoring function uses the parametrization of Meila.
 * 
 * @author Jens Keilwagen
 */
public class CMMScoringFunction extends AbstractVariableLengthScoringFunction implements SamplingScoringFunction
{
	private boolean freeParams, plugIn, optimize, optimizeFrame;

	private int order, period, starts, initFrame;
	private int[] powers;
	private double logGammaSum;
	private double[] frameHyper;
	
	// the index can be computed using the powers array
	private double[][][] params, probs, logNorm;
	private double[] frameLogScore, frameParams, frameProbs, partDer;
	private double[][][] hyper;
	private double logFrameNorm;
	private int[][][] counter, distCounter;
	private int[][] offset;

	/**
	 * This method returns the hyper-parameters for a model given some a-priori probabilities.
	 * 
	 * @param alphabetSize the size of the alphabet
	 * @param length the expected sequence length
	 * @param ess the equivalent sample size (ess) of the model
	 * @param frameProb the a-priori probabilities for each frame
	 * @param prob the a-priori probabilities for each frame and order
	 * 
	 * @return specific hyper-parameters
	 */
	public static double[][][] getHyperParams( int alphabetSize, int length, double ess, double[] frameProb, double[][][] prob ) {
		if( alphabetSize <= 0 ) {
			throw new IllegalArgumentException();
		}
		if( length <= 0 ) {
			throw new IllegalArgumentException();
		}
		if( ess <= 0 ) {
			throw new IllegalArgumentException();
		}
		if( frameProb.length != prob.length ) {
			throw new IllegalArgumentException();
		}
		
		int order = prob[0].length-1;
		double[][][] hyper = new double[prob.length][prob[0].length][];
		double[] currentHyper = new double[(int)Math.pow( alphabetSize, order )], nextHyper = new double[currentHyper.length];
		for( int startFrame = 0; startFrame < frameProb.length; startFrame++ ) {
			currentHyper[0] = ess*frameProb[startFrame];
			for( int ord, y, l = 0, f = startFrame; l < length; l++ ) {
				ord = Math.min( order, l );
				y = (int) Math.pow( alphabetSize, ord );
				if( hyper[f][ord] == null ) {
					hyper[f][ord] = new double[y*alphabetSize];
				}
				for( int x = 0; x < y; x++ ) {
					for( int v, a = 0; a < alphabetSize; a++ ){
						v = x*alphabetSize+a;
						hyper[f][ord][v] += currentHyper[x]*prob[f][ord][v];
						nextHyper[v%nextHyper.length] += currentHyper[x]*prob[f][ord][v];
					}
				}
				f = (f+1)%frameProb.length;
				System.arraycopy( nextHyper, 0, currentHyper, 0, currentHyper.length );
				Arrays.fill( nextHyper, 0 );
			}
		}
		
		return hyper;
	}
	
	private static double[][][] getHyper( int period, int alphabetSize, double[] sumOfHyper ) {
		double[][][] hyper = new double[period][sumOfHyper.length][];
		for( int a, o = 0; o < sumOfHyper.length; o++ ) {
				a = (int) Math.pow( alphabetSize, o+1 );
				double h = sumOfHyper[o] / period / a;
				for( int p = 0; p < period; p++ ) {
					hyper[p][o] = new double[a];
					Arrays.fill( hyper[p][o], h );
				}
		}
		return hyper;
	}
	
	private static double[] getHyper( int period, double ess ) {
		double[] hyper = new double[period];
		Arrays.fill( hyper, ess/period );
		return hyper;
	}
	
	/**
	 * The main constructor.
	 * 
	 * @param alphabets the alphabet container
	 * @param order the oder of the model (has to be non-negative)
	 * @param period the period
	 * @param classEss the ess of the class
	 * @param sumOfHyperParams the sum of the hyper parameter for each order (length has to be <code>order</code>+1, each entry has to be non-negative), the sum also sums over the period
	 * @param plugIn a switch which enables to used the MAP-parameters as plug-in parameters
	 * @param optimize a switch which enables to optimize or fix the parameters
	 * @param starts the number of recommended starts
	 * @param initFrame the frame which should be used for plug-in initialization, negative for random initialization
	 * 
	 * @see #getHyperParams(int, int, double, double[], double[][][])
	 * @see #CMMScoringFunction(AlphabetContainer, double[], double[][][], boolean, boolean, int, int)
	 */
	public CMMScoringFunction( AlphabetContainer alphabets, int order, int period, double classEss, double[] sumOfHyperParams, boolean plugIn, boolean optimize, int starts, int initFrame )
	{
		this( alphabets, getHyper( period, classEss ), getHyper( period, (int)alphabets.getAlphabetLengthAt( 0 ), sumOfHyperParams ), plugIn, optimize, starts, initFrame );
	}
	
	/**
	 * This constructor allows to create an instance with specific hyper-parameters for all conditional distributions.
	 * 
	 * @param alphabets the alphabet container
	 * @param frameHyper the hyper-parameters for the frame, the length of this array also defines the period of the model
	 * @param hyper the hyper-parameters for each frame 
	 * @param plugIn a switch which enables to used the MAP-parameters as plug-in parameters
	 * @param optimize a switch which enables to optimize or fix the parameters
	 * @param starts the number of recommended starts
	 * @param initFrame the frame which should be used for plug-in initialization, negative for random initialization
	 */
	public CMMScoringFunction( AlphabetContainer alphabets, double[] frameHyper, double[][][] hyper, boolean plugIn, boolean optimize, int starts, int initFrame ) {
		super( alphabets );
		this.order = hyper[0].length-1;
		this.period = frameHyper.length;		
		createArrays();
		
		this.frameHyper = new double[period];
		this.hyper = new double[period][order+1][];
		for( int p = 0; p < period; p++ ) {
			if( frameHyper[p] < 0 )
			{
				throw new IllegalArgumentException( "The ess for the class has to be non-negative." );
			}
			this.frameHyper[p] = frameHyper[p];
			for( int o = 0; o <= order; o++ ) {
				if( hyper[p][o].length != Math.pow(powers[1],o+1) ) {
					throw new IllegalArgumentException();//TODO
				}
				this.hyper[p][o] = new double[hyper[p][o].length];
				for( int i = 0; i < this.hyper[p][o].length; i++ ) {
					if( hyper[p][o][i] < 0 )
					{
						throw new IllegalArgumentException( "The ess for the class has to be non-negative." );
					}
					this.hyper[p][o][i] = hyper[p][o][i];
				}
			}
		}

		frameParams = new double[period];
		Arrays.fill( frameParams, -Math.log( period ) );
		Arrays.fill( frameProbs, 1d / (double) period );
		params = new double[period][order+1][];
		double uniform = 1d / (double) powers[1], logUniform = Math.log( uniform );
		for( int i, p = 0; p < period; p++ )
		{
			for( i = 0; i <= order; i++ )
			{
				params[p][i] = new double[powers[i+1]];
				probs[p][i] = new double[powers[i+1]];
				logNorm[p][i] = new double[powers[i]];
				Arrays.fill( params[p][i], logUniform );
				Arrays.fill( probs[p][i], uniform );
			}
		}
		this.plugIn = plugIn;
		this.optimize = optimize;
		this.optimizeFrame = true;
		if( starts <= 0 )
		{
			throw new IllegalArgumentException( "The number of starts has to be positive." );
		}
		this.starts = starts;
		setFreeParams( false );
		computeConstantsOfLogPrior();
		if( initFrame < period )
		{
			this.initFrame = initFrame;
		}
		else
		{
			throw new IllegalArgumentException( "Check initFrame." );
		}
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}.
	 * 
	 * @param source the xml representation
	 * 
	 * @throws NonParsableException if the representation could not be parsed.
	 */
	public CMMScoringFunction( StringBuffer source ) throws NonParsableException
	{
		super( source );
	}

	private void createArrays()
	{
		powers = new int[order+2];
		powers[0] = 1;
		powers[1] = (int) alphabets.getAlphabetLengthAt( 0 );
		int i;
		for( i = 2; i < powers.length; i++ )
		{
			powers[i] = powers[i-1]*powers[1];
		}
		probs = new double[period][order+1][];
		logNorm = new double[period][order+1][];
		frameProbs = new double[period];
		counter = new int[period][period][powers[order+1]];
		distCounter = new int[period][period][powers[order]];
		offset = new int[period][order+2];
		frameLogScore = new double[period];
		partDer = new double[powers[1]];
	}
	
	public CMMScoringFunction clone() throws CloneNotSupportedException
	{
		CMMScoringFunction clone = (CMMScoringFunction) super.clone();
		// the powers do not have to be cloned
		clone.frameHyper = frameHyper.clone();
		clone.hyper = ArrayHandler.clone( hyper );
		
		
		clone.params = new double[period][order+1][];
		clone.probs = new double[period][order+1][];
		clone.logNorm = new double[period][order+1][];
		clone.offset = new int[period][];
		clone.counter = new int[period][period][];
		clone.distCounter = new int[period][period][];
		for( int o, p = 0; p < period; p++ )
		{
			for( o = 0; o <= order; o++ )
			{
				clone.params[p][o] = params[p][o].clone();
				clone.probs[p][o] = probs[p][o].clone();
				clone.logNorm[p][o] = logNorm[p][o].clone();
			}
			clone.offset[p] = offset[p].clone();
			for( o = 0; o < period; o++ )
			{
				clone.counter[p][o] = counter[p][o].clone();
				clone.distCounter[p][o] = distCounter[p][o].clone();
			}
		}
		clone.frameParams = frameParams.clone();
		clone.frameProbs = frameProbs.clone();
		clone.frameLogScore = frameLogScore.clone();
		clone.partDer = partDer.clone();
		return clone;
	}

	public String getInstanceName()
	{
		return "cMM(" + order + ", " + period + ")";
	}
	
	private void fillFrameLogScores( Sequence seq, int start, int length ) {
		int l = 0, indexOld, indexNew = 0, o = Math.min( order, length ), p;
		for( p = 0; p < period; p++ )
		{
			frameLogScore[p] = frameParams[p];
		}
		for( ; l < o; l++ )
		{
			indexOld = indexNew;
			indexNew = indexOld*powers[1] + seq.discreteVal( start++ );
			for( p = 0; p < period; p++ )
			{
				frameLogScore[p] += params[(p+l)%period][l][indexNew] - logNorm[(p+l)%period][l][indexOld];
			}
		}
		for( ; l < length; l++ )
		{
			indexOld = indexNew % powers[order];
			indexNew = indexOld * powers[1] + seq.discreteVal( start++ );
			for( p = 0; p < period; p++ )
			{
				frameLogScore[p] += params[(p+l)%period][order][indexNew] - logNorm[(p+l)%period][order][indexOld];
			}
		}
	}
	
	public double getLogScore( Sequence seq, int start, int length )
	{
		fillFrameLogScores(seq, start, length);
		return Normalisation.getLogSum( frameLogScore ) - logFrameNorm;
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int start, int length, IntList indices, DoubleList dList )
	{
		if( optimize )
		{
			int l = 0, indexOld, indexNew = 0, h, o = Math.min( order, length ), index, z, p;
			for( p = 0; p < period; p++ )
			{
				for( z = 0; z < period; z++ )
				{
					Arrays.fill( counter[p][z], 0 );
					Arrays.fill( distCounter[p][z], 0 );
				}
				frameLogScore[p] = frameParams[p];
			}
			int stop = powers[1] - (freeParams?1:0);
			
			// start probablities
			for( ; l < o; l++ )
			{
				indexOld = indexNew;
				z = indexOld * powers[1];
				indexNew = z + seq.discreteVal( start++ );
				h = z - (freeParams?indexOld:0);
				
				for( p = 0; p < period; p++ )
				{
					frameLogScore[p] += params[(p+l)%period][l][indexNew] - logNorm[(p+l)%period][l][indexOld];
					for( index = 0; index < stop; index++ )
					{
						indices.add( offset[p][l] + h + index );
						if( z + index == indexNew )
						{
							dList.add( 1 - probs[p][l][z + index] );
						}
						else
						{
							dList.add( - probs[p][l][z + index] );
						}
					}
				}
			}
			//counting the usage of transition probability parameters for the sequence
			for( ; l < length; l++ )
			{
				indexOld = indexNew % powers[order];
				indexNew = indexOld * powers[1] + seq.discreteVal( start++ );
				for( p = 0; p < period; p++ )
				{
					frameLogScore[p] += params[(p+l)%period][order][indexNew] - logNorm[(p+l)%period][order][indexOld];
					distCounter[p][(p+l)%period][indexOld]++;
					counter[p][(p+l)%period][indexNew]++;
				}
			}
			
			//computing the score
			double erg = Normalisation.logSumNormalisation( frameLogScore, 0, period, frameLogScore, 0 ) - logFrameNorm;
			
			
			//computing the gradient and the score

			//start probs
			indexOld = 0;
			for( l = 0; l < o; l++ )
			{
				for( p = 0; p < period; p++ )
				{
					indexNew = indexOld + stop;
					dList.multiply( indexOld, indexNew, frameLogScore[p] );
					indexOld = indexNew;
				}
			}
			
			boolean used;
			//transition probs
			for( l = 0; l < distCounter[0][0].length; l++ )
			{
				h = l*(powers[1]-(freeParams?1:0));
				o = l*powers[1];
				for( z = 0; z < period; z++ )
				{
					Arrays.fill( partDer, 0 );
					used = false;
					for( p = 0; p < period; p++ )
					{
						if( distCounter[p][z][l] > 0 )
						{
							used = true;
							for( index = 0; index < stop; index++ )
							{
								partDer[index] += (frameLogScore[p] *( counter[p][z][o+index] - distCounter[p][z][l]*probs[z][order][o+index] ));
							}
						}
					}
					if( used )
					{
						for( index = 0; index < stop; index++ )
						{
							indices.add( offset[z][order] + h + index  );
							dList.add( partDer[index] );
						}
					}
				}
			}
			
			if( optimizeFrame )
			{
				for( p = 0; p < period-(freeParams?1:0); p++ )
				{
					indices.add( p );
					dList.add( frameLogScore[p] - frameProbs[p] );
				}
			}
			return erg;
		}
		else
		{
			return getLogScore( seq, start, length );
		}
	}

	public int getNumberOfParameters()
	{
		return offset[period-1][order+1];
	}

	public void setParameters( double[] params, int start )
	{
		if( optimize )
		{
			int j, n, index, o, p, stop;
			if( optimizeFrame )
			{
				stop = period - (freeParams?1:0);
				logFrameNorm = 0;
				for( p = 0; p < stop; p++ )
				{
					frameParams[p] = params[start++];
					frameProbs[p] = Math.exp( frameParams[p] );
					logFrameNorm += frameProbs[p];
				}
				if( stop < period )
				{
					frameProbs[p] = Math.exp( frameParams[p] );
					logFrameNorm += frameProbs[p];
				}
				for( p = 0; p < period; p++ )
				{
					frameProbs[p] /= logFrameNorm;
				}
				logFrameNorm = Math.log( logFrameNorm );
			}
			stop = powers[1] - (freeParams?1:0);
			for( p = 0; p < period; p++ )
			{
				for( o = 0; o <= order; o++ )
				{
					for( index = n = 0; n < logNorm[p][o].length; n++ )
					{
						logNorm[p][o][n] = 0d;
						for( j = 0; j < stop; j++, start++ )
						{
							this.params[p][o][index+j] = params[start]; 
							probs[p][o][index+j] = Math.exp( this.params[p][o][index+j] );
							logNorm[p][o][n] += probs[p][o][index+j]; 
						}
						if( j < powers[1] )
						{
							probs[p][o][index+j] = Math.exp( this.params[p][o][index+j] );
							logNorm[p][o][n] += probs[p][o][index+j];
						}
						for( j = 0; j < powers[1]; j++, index++ )
						{
							this.probs[p][o][index] /= logNorm[p][o][n]; 
						}
						logNorm[p][o][n] = Math.log( logNorm[p][o][n] );
					}
				}
			}
		}
	}

	public StringBuffer toXML()
	{
		StringBuffer b = new StringBuffer( 10000 );
		XMLParser.appendObjectWithTags( b, length, "length" );
		XMLParser.appendObjectWithTags( b, alphabets, "alphabets" );
		XMLParser.appendObjectWithTags( b, order, "order" );
		XMLParser.appendObjectWithTags( b, period, "period" );
		XMLParser.appendObjectWithTags( b, frameHyper, "frameEss" );
		XMLParser.appendObjectWithTags( b, hyper, "hyper" );
		XMLParser.appendObjectWithTags( b, frameParams, "frameParams" );
		for( int p = 0; p < period; p++ )
		{
			XMLParser.appendObjectWithTagsAndAttributes( b, params[p], "params", "frame=\"" + p + "\"" );
		}
		XMLParser.appendObjectWithTags( b, plugIn, "plugIn" );
		XMLParser.appendObjectWithTags( b, optimize, "optimize" );
		XMLParser.appendObjectWithTags( b, optimizeFrame, "optimizeFrame" );
		XMLParser.appendObjectWithTags( b, starts, "starts" );
		XMLParser.appendObjectWithTags( b, freeParams, "freeParams" );
		XMLParser.appendObjectWithTags( b, initFrame, "initFrame" );
		XMLParser.addTags( b, getClass().getSimpleName() );
		return b;
	}

	public double[] getCurrentParameterValues()
	{
		int l = optimize?offset[period-1][order+1]:0;
		double[] erg = new double[l];
		if( optimize )
		{
			int stop, p, i = 0, j, index, o;
			if( optimizeFrame )
			{
				stop = period - (freeParams?1:0);
				for( p = 0; p < stop; p++, i++ )
				{
					erg[i] = frameParams[p];
				}
			}
			stop = powers[1] - (freeParams?1:0);
			for( p = 0; p < period; p++ )
			{
				for( o = 0; o <= order; o++ )
				{
					for( index = 0; index < params[p][o].length; index += powers[1] )
					{
						for( j = 0; j < stop; j++, i++ )
						{
							erg[i] = params[p][o][index+j];
						}
					}
				}
			}
		}
		return erg;
	}
	
	public void initializeFunction( int index, boolean freeParams, Sample[] data, double[][] weights )
	{
		if( optimize && plugIn && data != null && data[index] != null )
		{
			int anz = data[index].getNumberOfElements();
			Sequence seq;
			double w = 1;
			boolean externalWeights = weights != null && weights [index] != null;
			int rMax = initFrame >= 0 ? 1 : 3;
			double[][] frameP = new double[anz][period];
			if( initFrame < 0 ) {
				initializeFunctionRandomly( freeParams );
			}
			for( int r = 0; r < rMax; r++ ) {
				for( int i = 0; i < anz; i++ ) {
					if( initFrame < 0 ) {
						fillFrameLogScores(data[index].getElementAt(i), 0, length);
						Normalisation.logSumNormalisation(frameLogScore, 0, period, frameP[i], 0);
					} else {
						Arrays.fill( frameP[i], 0 );
						frameP[i][initFrame] = 1;
					}
				}
				
				//preparation
				int len, o, indexOld, indexNew, p;
				for( p = 0; p < period; p++ )
				{
					for( o = 0; o <= order; o++ )
					{
						System.arraycopy( hyper[p][o], 0, params[p][o], 0, hyper[p][o].length );
						for( int idx = 0, n= 0; idx < params[p][o].length; n++, idx += powers[1] )
						{
							logNorm[p][o][n] = 0;
							for( int j = 0; j < powers[1]; j++ )
							{
								logNorm[p][o][n] += hyper[p][o][idx+j];
							}
						}
					}
				}
				if( optimizeFrame )
				{
					System.arraycopy( frameHyper, 0, frameProbs, 0, period );
					logFrameNorm = getEss();
				}
				
				//counting
				for( int l, i = 0; i < anz; i++ )
				{
					seq = data[index].getElementAt(i);
					len = seq.getLength();
					o = Math.min( len, order );
					indexNew = 0;
					if( externalWeights )
					{
						w = weights[index][i];
					}
					
					if( optimizeFrame )
					{
						for( p = 0; p < period; p++ )
						{
							frameP[i][p] *= w;
							frameProbs[p] += frameP[i][p];
						}
						logFrameNorm += w;
					}
					
					for( l = 0 ; l < o; l++ )
					{
						indexOld = indexNew;
						indexNew = indexOld*powers[1] + seq.discreteVal( l );
						for( p = 0; p < period; p++ )
						{
							probs[(p+l)%period][l][indexNew] += frameP[i][p];
							logNorm[(p+l)%period][l][indexOld] += frameP[i][p];
						}
					}
					for( ; l < len; l++ )
					{
						indexOld = indexNew % powers[order];
						indexNew = indexOld * powers[1] + seq.discreteVal( l );
						for( p = 0; p < period; p++ )
						{
							probs[(p+l)%period][order][indexNew] += frameP[i][p];
							logNorm[(p+l)%period][order][indexOld] += frameP[i][p];
						}
					}
				}
				
				//computing freqs and parameters
				if( optimizeFrame )
				{
					for( p = 0; p < period; p++ )
					{
						frameProbs[p] /= logFrameNorm;
						frameParams[p] = Math.log( frameProbs[p] );
					}
					logFrameNorm = 0;
				}
				for( p = 0; p < period; p++ )
				{
					for( o = 0; o <= order; o++ )
					{
						for( indexOld = indexNew = 0; indexOld < logNorm[p][o].length; indexOld++ )
						{
							for( len = 0; len < powers[1]; len++, indexNew++ )
							{
								probs[p][o][indexNew] /= logNorm[p][o][indexOld];
								params[p][o][indexNew] = Math.log( probs[p][o][indexNew]);
							}
							logNorm[p][o][indexOld] = 0;
						}
					}
				}
			}
		}
		else
		{
			initializeFunctionRandomly( freeParams );
		}
		setFreeParams( freeParams );	
	}
	
	public void initializeFunctionRandomly( boolean freeParams )
	{
		if( optimize )
		{
			int o, normCounter, paramCounter, len, p;
			DirichletMRGParams hyper;
			double[] freq, currentHyper = new double[powers[1]];
			if( optimizeFrame )
			{
				hyper = new DirichletMRGParams( frameHyper );
				freq = DirichletMRG.DEFAULT_INSTANCE.generate( period, hyper );
				for( p = 0; p < period; p++ )
				{
					frameProbs[p] = freq[p];
					frameParams[p] = Math.log( freq[p] );
				}
				logFrameNorm = 0;
			}
			freq = new double[powers[1]];
			for( p = 0; p < period; p++ )
			{
				for( o = 0; o <= order; o++ )
				{
					for( normCounter = paramCounter = 0; normCounter < logNorm[p][o].length; normCounter++ )
					{
						logNorm[p][o][normCounter] = 0;
						for( len = 0; len < powers[1]; len++ )
						{
							currentHyper[len] = this.hyper[p][o][paramCounter+len];
						}
						hyper = new DirichletMRGParams( currentHyper );
						DirichletMRG.DEFAULT_INSTANCE.generate( freq, 0, powers[1], hyper );
						for( len = 0; len < powers[1]; len++, paramCounter++ )
						{
							probs[p][o][paramCounter] = freq[len];
							params[p][o][paramCounter] = Math.log( freq[len] );
						}
					}
				}
			}
			setFreeParams( freeParams );
		}
	}

	protected void fromXML( StringBuffer xml ) throws NonParsableException
	{
		StringBuffer b = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		length = XMLParser.extractObjectForTags( b, "length", int.class );
		alphabets = XMLParser.extractObjectForTags( b, "alphabets", AlphabetContainer.class );
		order = XMLParser.extractObjectForTags( b, "order", int.class );
		period = XMLParser.extractObjectForTags( b, "period", int.class );
		createArrays();
		StringBuffer help = XMLParser.extractForTag( b, "frameEss" );
		if( help == null ) {
			//old
			frameHyper = getHyper( period, XMLParser.extractObjectForTags( b, "classEss", double.class )/period );
			hyper = getHyper( period, powers[1], XMLParser.extractObjectForTags( b, "sumOfHyperParams", double[].class ) );
			
		} else {
			//new
			XMLParser.addTags( help, "frameEss" );
			frameHyper = (double[]) XMLParser.extractObjectForTags( help, "frameEss" );
			hyper = (double[][][]) XMLParser.extractObjectForTags( b, "hyper" );
		}
		frameParams = XMLParser.extractObjectForTags( b, "frameParams", double[].class );
		logFrameNorm = 0;
		int j, n, index, o, p = 0;
		params = new double[period][][];
		Map<String, String> map = new TreeMap<String, String>();
		for( ; p < period; p++ )
		{
			frameProbs[p] = Math.exp( frameParams[p] );
			logFrameNorm += frameProbs[p];
			map.clear();
			map.put( "frame", ""+p );
			params[p] = XMLParser.extractObjectAndAttributesForTags( b, "params", null, map, double[][].class );
			for(o = 0; o <= order; o++ )
			{
				probs[p][o] = new double[params[p][o].length];
				logNorm[p][o] = new double[powers[o]];
				for( n = index = 0; n < logNorm[p][o].length; n++ )
				{
					logNorm[p][o][n] = 0d;
					for( j = 0; j < powers[1]; j++ )
					{
						probs[p][o][index+j] = Math.exp( params[p][o][index+j] );
						logNorm[p][o][n] += probs[p][o][index+j]; 
					}
					for( j = 0; j < powers[1]; j++, index++ )
					{
						this.probs[p][o][index] /= logNorm[p][o][n]; 
					}
					logNorm[p][o][n] = Math.log( logNorm[p][o][n] );
				}
			}
		}
		for( p = 0; p < period; p++ )
		{
			frameProbs[p] /= logFrameNorm;
		}
		logFrameNorm = Math.log( logFrameNorm );
		
		plugIn = XMLParser.extractObjectForTags( b, "plugIn", boolean.class );
		optimize = XMLParser.extractObjectForTags( b, "optimize", boolean.class );
		optimizeFrame = XMLParser.extractObjectForTags( b, "optimizeFrame", boolean.class );
		starts = XMLParser.extractObjectForTags( b, "starts", int.class );
		setFreeParams( XMLParser.extractObjectForTags( b, "freeParams", boolean.class ) );
		initFrame = XMLParser.extractObjectForTags( b, "initFrame", int.class );
		computeConstantsOfLogPrior();
	}
	
	private void setFreeParams( boolean freeParams )
	{
		this.freeParams = freeParams;
		if( optimize )
		{
			for( int o, p = 0; p < period; p++ )
			{
				if( p == 0 )
				{
					offset[0][0] = optimizeFrame?period -(freeParams?1:0):0;
				}
				else
				{
					offset[p][0] = offset[p-1][order+1];
				}
				for( o = 0; o <= order; o++ )
				{
					offset[p][o+1] = offset[p][o] + params[p][o].length - (freeParams?powers[o]:0);
				}
			}
		}
		else
		{
			offset[period-1][order+1] = 0;
		}
	}

	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index )
	{
		if( index < offset[0][0] ) {
			return period;
		} else {
			int p = 0;
			while( index >= offset[p][order+1] ) {
				p++;
			}
			int o = 1;
			while( index >= offset[p][o] ) {
				o++;
			}
			return powers[o-1];
		}
	}

	public double getLogNormalizationConstant( int length )
	{
		return 0;
	}
	
	public double getLogPartialNormalizationConstant( int parameterIndex, int length ) throws Exception
	{
		if( parameterIndex < offset[period-1][order+1] )
		{
			return Double.NEGATIVE_INFINITY;
		}
		else
		{
			throw new IndexOutOfBoundsException();
		}
	}

	public double getEss()
	{
		double ess = 0;
		for( int f = 0; f < frameHyper.length; f++ ) {
			ess += frameHyper[f];
		}
		return ess;
	}
	
	public String toString()
	{
		DiscreteAlphabet abc =(DiscreteAlphabet) alphabets.getAlphabetAt( 0 ); 
		int i = 0, o, p, index, l = (int)abc.length();
		StringBuffer info = new StringBuffer( (int) Math.pow(l,order) * l * period * 15 );
		
		String[] sym = new String[l];
		l--;
		for( ; i <= l; i++ )
		{
			sym[i] = abc.getSymbolAt(i);
		}
		int[] context = new int[order+1];
		for( p = 0; p < period; p++ )
		{
			info.append( "frame " + p + ": p("+p+") = " + frameProbs[p] + "\n" );
			for( i = 0; i <= l; i++ )
			{
				info.append( "\t" + sym[i] );
			}
			info.append( "\n" );
			for( o = 0; o <= order; o++ )
			{
				info.append( "P(X_" + o );
				for( i = 0; i < o; i++ )
				{
					if( i == 0 )
					{
						info.append( "|" );
					}
					else
					{
						info.append( " " );
					}
					info.append( "X_" + i );
				}
				info.append( ")\n" );
				Arrays.fill( context, 0 );
				for( index = 0; index < probs[p][o].length; )
				{
					for( i = 0; i < o; i++ )
					{
						info.append( sym[context[i]] );
					}
					for( i = 0; i <= l; i++, index++ )
					{
						info.append( "\t" + probs[p][o][index] + "\t("+hyper[p][o][index]+")" );
					}
					info.append( "\n" );
					
					i = o-1;
					while( i >= 0 && context[i] == l )
					{
						context[i] = 0;
						i--;
					}
					if( i >= 0 )
					{
						context[i]++;
					}
				}
				info.append( "\n" );
			}
		}
		return info.toString();
		/**/
	}
	
	public double getLogPriorTerm()
	{
		if( optimize )
		{
			double classESS = getEss(), sum, val = -classESS*logFrameNorm;//XXX -?
			int A = (int)alphabets.getAlphabetLengthAt( 0 );
			for( int f = 0; f < params.length; f++ ) {
				 val += frameParams[f] * frameHyper[f];
				 for( int o = 0; o < params[f].length; o++ ) {
					 sum = 0;
					 for( int a = 0, i = 0; i < params[f].length; i++ ) {
						 val += params[f][o][i] * hyper[f][o][i];
						 sum += hyper[f][o][i];
						 a++;
						 if( a == A ) {
							 val -= logNorm[f][o][i/A] * sum;
							 sum = 0;
							 a = 0;
						 }
					 }
				 }
			}
			return val + logGammaSum;
		}
		else
		{
			return 0;
		}
	}
	
	private double getLogPriorTerm( int offset )
	{
		if( optimize )
		{
			double classESS = getEss(), sum, val = -classESS*logFrameNorm;//XXX -?
			int A = (int)alphabets.getAlphabetLengthAt( 0 );
			for( int f = 0; f < params.length; f++ ) {
				 val += frameParams[f] * frameHyper[(f+offset)%period];
				 for( int o = 0; o < params[f].length; o++ ) {
					 sum = 0;
					 for( int a = 0, i = 0; i < params[f].length; i++ ) {
						 val += params[f][o][i] * hyper[(f+offset)%period][o][i];
						 sum += hyper[(f+offset)%period][o][i];
						 a++;
						 if( a == A ) {
							 val -= logNorm[f][o][i/A] * sum;
							 sum = 0;
							 a = 0;
						 }
					 }
				 }
			}
			return val + logGammaSum;
		}
		else
		{
			return 0;
		}
	}

	private void computeConstantsOfLogPrior()
	{
		double classESS = getEss(), sum = 0;
		int A = (int)alphabets.getAlphabetLengthAt( 0 );
		logGammaSum = Gamma.logOfGamma( classESS );
		for( int f = 0; f < params.length; f++ ) {
			 logGammaSum -= Gamma.logOfGamma( frameHyper[f] );
			 for( int o = 0; o < params[f].length; o++ ) {
				 for( int a = 0, i = 0; i < params[f][o].length; i++ ) {
					 logGammaSum -= Gamma.logOfGamma( hyper[f][o][i] );
					 sum += hyper[f][o][i];
					 a++;
					 if( a == A ) {
						 logGammaSum = Gamma.logOfGamma( sum );
						 sum = 0;
						 a = 0;
					 }
				 }
			 }
		}
	}
	
	public void addGradientOfLogPriorTerm( double[] grad, int start )
	{
		if( optimize )
		{
			int j, p, o, index;
			if( optimizeFrame )
			{
				double classESS = getEss();
				j = period - (freeParams?1:0);
				for( p = 0; p < j; p++, start++ )
				{
					grad[start] += frameHyper[p] - classESS * frameProbs[p];
				}
			}
			
			double sum = 0;
			int stop = powers[1] - (freeParams?1:0);
			for( p = 0; p < params.length; p++ ) {
				for( o = 0; o < params[p].length; o++ ) { 
					for( index = 0; index < params[p][o].length; index += powers[1] )
					{
						sum = 0;
						for( j = 0; j < powers[1]; j++ )
						{
							sum += hyper[p][o][index+j];
						}
						for( j = 0; j < stop; j++, start++ )
						{
							grad[start] += hyper[p][o][index+j] - sum*probs[p][o][index+j]; 
						}
					}
				}
			}
			//System.out.println( start );System.exit(0);
		}
	}
	
	public boolean isNormalized()
	{
		return true;
	}

	public boolean isInitialized()
	{
		return true;
	}
	
	public int getNumberOfRecommendedStarts()
	{
		return starts;
	}
	
	/**
	 * This method enables the user to choose whether the parameters should be optimized or not.
	 *  
	 * @param optimize the switch for optimization of the parameters
	 */
	public void setParameterOptimization( boolean optimize )
	{
		this.optimize = optimize;
		setFreeParams( freeParams );
	}
	
	/**
	 * This method enables the user to choose whether the frame parameters should be optimized or not.
	 *  
	 * @param optimize the switch for optimization of the frame parameters
	 */
	public void setFrameParameterOptimization( boolean optimize )
	{
		this.optimizeFrame = optimize;
		setFreeParams( freeParams );
	}
	
	/*
	public void show()
	{
		for( int p = 0; p < period; p++ )
		{
			System.out.println( p + "\t" + Arrays.toString( offset[p] ) );
		}
	}*/
	
	public void setStatisticForHyperparameters( int[] length, double[] weight ) throws Exception
	{
		/*
		if( weight.length != length.length )
		{
			throw new IllegalArgumentException( "The length of both arrays (length, weight) have to be identical." );
		}
		Arrays.fill( sumOfHyperParams, 0 );
		for( int l, i = 0; i < length.length; i++ )
		{
			if( weight[i] < 0 || length[i] < 0 )
			{
				throw new IllegalArgumentException( "check length and weight for entry " + i );
			}
			else
			{
				for( l = 0; l < length[i] && l < order; l++ )
				{
					sumOfHyperParams[l] += weight[i];
				}
				if( order < length[i] )
				{
					sumOfHyperParams[order] += (length[i]-order) * weight[i];
				}
			}
		}
		computeConstantsOfLogPrior();
		*/
		//XXX
		throw new RuntimeException();

	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.SamplingScoringFunction#getSamplingGroups(int)
	 */
	@Override
	public int[][] getSamplingGroups(int parameterOffset) {
		if( optimize ) {
			int[][] res = new int[optimizeFrame?1:0 + (offset[period-1][order+1] % powers[1])][];
			int i = 0;
			if( optimizeFrame ) {
				res[i] = new int[period];
				for( int j = 0; j < res[i].length; j++, parameterOffset++ ) {
					res[i][j] = parameterOffset;
				}
				i++;
			}
			for( ; i < res.length; i++ ) {
				res[i] = new int[powers[1]];
				for( int j = 0; j < res[i].length; j++, parameterOffset++ ) {
					res[i][j] = parameterOffset;
				}
			}
			return res;
		} else {
			return new int[0][];
		}
	}
}
