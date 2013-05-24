/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.differentiable.mixture;

import java.text.NumberFormat;
import java.util.Arrays;

import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.ComplementableDiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.Mutable;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.NormalizedDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;

/**
 * This class enables the user to search on both strand. So the motif can be found on the forward or on the reverse
 * complementary strand.
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.data.alphabets.ComplementableDiscreteAlphabet
 * @see de.jstacs.data.AlphabetContainer#isReverseComplementable()
 */
public class StrandDiffSM extends AbstractMixtureDiffSM implements Mutable
{
	/**
	 * This enum defines the different types of plug-in initialization of a {@link StrandDiffSM}.
	 * 
	 * @author Jens Keilwagen
	 */
	public enum InitMethod{
		/**
		 * This value indicates that the model is initialized with the data on the forward strand.
		 */
		INIT_FORWARD_STRAND,
		/**
		 * This value indicates that the model is initialized with the data on the backward strand.
		 */
		INIT_BACKWARD_STRAND,
		/**
		 * This value indicates that the model is initialized with the data, where the weight for both strand of each sequence is chosen randomly.
		 */
		INIT_BOTH_STRANDS;
	}
	
	private InitMethod initMethod;	
	private double forwardPartOfESS;

	private StrandDiffSM( DifferentiableStatisticalModel function, int starts, boolean optimizeHidden,
			boolean plugIn, double forwardPartOfESS, InitMethod initMethod ) throws CloneNotSupportedException, WrongAlphabetException
	{
		super( function.getLength(), starts, 2, optimizeHidden, plugIn, new DifferentiableStatisticalModel[]{ function } );
		if( !function.getAlphabetContainer().isReverseComplementable() )
		{
			throw new WrongAlphabetException( "The given AlphabetContainer can not be used for building a reverse complement." );
		}
		if( forwardPartOfESS < 0 || forwardPartOfESS > 1 )
		{
			throw new IllegalArgumentException( "The part of the ESS for the forward strand has to be in [0,1]." );
		}
		this.forwardPartOfESS = forwardPartOfESS;
		this.initMethod = initMethod;
		computeLogGammaSum();
	}

	/**
	 * This constructor creates a StrandDiffSM that optimizes the usage of each strand.
	 * 
	 * @param function
	 *            the DifferentiableSequenceScore
	 * @param forwardPartOfESS
	 *            the part of the full ESS that should be used as hyperparameter for the forward strand
	 * @param starts
	 *            the number of starts the should be done in an optimization
	 * @param plugIn
	 *            whether the initial parameters for an optimization should be related to the data or randomly drawn
	 * @param initMethod
	 *            only used if <code>plugIn==true</code><br>
	 *            whether the initial parameters for an optimization should be related to the data of the forward strand,
	 *            the backward strand or both strands
	 *            
	 * @throws CloneNotSupportedException if <code>function</code> could not be cloned
	 * @throws WrongAlphabetException if the alphabet of <code>function</code> is not {@link de.jstacs.data.AlphabetContainer#isReverseComplementable()} and, hence, cannot be used for a strand mixture
	 * 
	 * @see StrandDiffSM.InitMethod
	 */
	public StrandDiffSM( DifferentiableStatisticalModel function, double forwardPartOfESS, int starts,
			boolean plugIn, InitMethod initMethod ) throws CloneNotSupportedException, WrongAlphabetException
	{
		this( function, starts, true, plugIn, forwardPartOfESS, initMethod );
	}

	/**
	 * This constructor creates a StrandDiffSM that has a fixed frequency for the strand usage.
	 * 
	 * @param function
	 *            the DifferentiableSequenceScore
	 * @param starts
	 *            the number of starts the should be done in an optimization
	 * @param plugIn
	 *            whether the initial parameters for an optimization should be related to the data or randomly drawn
	 * @param initMethod
	 *            only used if <code>plugIn==true</code><br>
	 *            whether the initial parameters for an optimization should be related to the data of the forward strand,
	 *            the backward strand or both strands
	 * @param forward
	 *            the probability of a motif to be on the forward strand
	 *            
	 * @throws CloneNotSupportedException if <code>function</code> could not be cloned
	 * @throws WrongAlphabetException if the alphabet of <code>function</code> is not {@link de.jstacs.data.AlphabetContainer#isReverseComplementable()} and, hence, cannot be used for a strand mixture
	 * 
	 * @see StrandDiffSM.InitMethod
	 */
	public StrandDiffSM( DifferentiableStatisticalModel function, int starts, boolean plugIn, InitMethod initMethod, double forward )
			throws CloneNotSupportedException, WrongAlphabetException
	{
		this( function, starts, false, plugIn, forward, initMethod );
		if( forward < 0 || forward > 1 )
		{
			throw new IllegalArgumentException( "The value for forward is no probability." );
		}
		setForwardProb(forward);
	}
	
	/**
	 * This method can be used to set the forward strand probability.
	 * 
	 * @param forward the forward strand probability in (0,1)
	 */
	protected void setForwardProb( double forward )
	{
		hiddenParameter[0] = Math.log( forward );
		hiddenParameter[1] = Math.log( 1d - forward );
		setHiddenParameters( hiddenParameter, 0 );
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}.
	 * 
	 * @param xml the xml representation
	 * 
	 * @throws NonParsableException if the representation could not be parsed.
	 */
	public StrandDiffSM( StringBuffer xml ) throws NonParsableException
	{
		super( xml );
	}

	protected double getLogNormalizationConstantForComponent( int i )
	{
		return function[0].getLogNormalizationConstant();
	}

	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception
	{
		if( isNormalized() )
		{
			return Double.NEGATIVE_INFINITY;
		}
		else
		{
			if( Double.isNaN( norm ) )
			{
				precomputeNorm();
			}
			int[] ind = getIndices( parameterIndex );
			if( ind[0] == 1 )
			{
				return logHiddenPotential[ind[1]] + function[0].getLogNormalizationConstant();
			}
			else
			{
				return Normalisation.getLogSum( logHiddenPotential )
						+ function[ind[0]].getLogPartialNormalizationConstant( ind[1] );
			}
		}
	}

	public double getHyperparameterForHiddenParameter( int index )
	{
		switch( index )
		{
			case 0:
				return forwardPartOfESS * function[0].getESS();
			case 1:
				return (1d - forwardPartOfESS) * function[0].getESS();
			default:
				throw new IndexOutOfBoundsException();
		}
	}

	/**
	 * This methoth returns the a-priori probability for the forward strand.
	 *
	 * @return the a-priori probability for the forward strand
	 */
	public double getForwardProbability(){
		double d = hiddenPotential[0]+hiddenPotential[1];
		return (hiddenPotential[0]/d);
	}
	
	public double getESS()
	{
		return function[0].getESS();
	}

	protected void initializeUsingPlugIn( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception
	{
		DataSet myData = data[index];
		
		double[] stat = new double[2];
		switch( initMethod )
		{
			case INIT_BOTH_STRANDS:
				double p;
				if( optimizeHidden )
				{
					double[] h = new double[2];
					if( getESS() == 0 )
					{
						h[0] = h[1] = 1;
					}
					else
					{
						h[0] = getHyperparameterForHiddenParameter( 0 );
						h[1] = getHyperparameterForHiddenParameter( 1 );
					}			
					p = DirichletMRG.DEFAULT_INSTANCE.generate( 2, new DirichletMRGParams( h ) )[0]; 
				}
				else
				{
					p = hiddenPotential[0] / (hiddenPotential[0] + hiddenPotential[1]);
				}
				if( myData != null )
				{
					Sequence[] seqs = new Sequence[myData.getNumberOfElements()];
					double w = 1;
					//randomly;
					for( int i = 0; i < seqs.length; i++ )
					{
						if( weights != null && weights[index] != null ) {
							w = weights[index][i];
						}
						if( r.nextDouble() < p )
						{
							stat[0] += w;
							seqs[i] = myData.getElementAt( i );
						}
						else
						{
							stat[1] += w;
							seqs[i] = myData.getElementAt( i ).reverseComplement();
						}
					}
					data[index] = new DataSet( "randomly strand scrambled", seqs );
					function[0].initializeFunction( index, freeParams, data, weights );
					if( optimizeHidden )
					{
						computeHiddenParameter( stat, true );
					}
					//specific
					Arrays.fill( stat, 0 );
					for( int i = 0; i < seqs.length; i++ )
					{
						if( weights != null && weights[index] != null ) {
							w = weights[index][i];
						}
						if( getIndexOfMaximalComponentFor( myData.getElementAt( i ), 0 ) == 0 )
						{
							stat[0] += w;
							seqs[i] = myData.getElementAt( i );
						}
						else
						{
							stat[1] += w;
							seqs[i] = myData.getElementAt( i ).reverseComplement();
						}
					}
					data[index] = new DataSet( "strand scrambled", seqs );					
				}
				else
				{
					data [index] = null;
				}
				break;
			case INIT_BACKWARD_STRAND:
				Sequence[] rcs = new Sequence[myData.getNumberOfElements()];
				for( int i = 0; i < rcs.length; i++ )
				{
					rcs[i] = myData.getElementAt( i ).reverseComplement();
				}
				data[index] = new DataSet( "backward strand", rcs );
			default:
				stat[0] = forwardPartOfESS;
				stat[1] = 1d-forwardPartOfESS;
				break;
		}
		function[0].initializeFunction( index, freeParams, data, weights );
		data[index] = myData;
		if( optimizeHidden )
		{
			computeHiddenParameter( stat, true );
		}
	}

	public String getInstanceName()
	{
		String erg = "strand-mixture(" + function[0].getInstanceName();
		if( !optimizeHidden )
		{
			erg += ", " + Arrays.toString( hiddenPotential );
		}
		return erg + ")";
	}

	protected void fillComponentScores( Sequence seq, int start )
	{
		componentScore[0] = logHiddenPotential[0] + function[0].getLogScoreFor( seq, start );
		try
		{
			if( length != 0 )
			{
				componentScore[1] = logHiddenPotential[1]
						+ function[0].getLogScoreFor( seq.reverseComplement(), seq.getLength() - start - length );
			}
			else
			{
				if( start == 0 )
				{
					componentScore[1] = logHiddenPotential[1] + function[0].getLogScoreFor( seq.reverseComplement(), 0 );
				}
				else
				{
					//FIXME
					//componentScore[1] =;
					throw new Exception( "strand scoring for variable length function" );
				}
			}
		}
		catch( Exception doesNotHappen )
		{
			RuntimeException r = new RuntimeException( doesNotHappen.getClass().getName() + ": " + doesNotHappen.getMessage() );
			r.setStackTrace( doesNotHappen.getStackTrace() );
			throw r;
		}
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer )
	{
		iList[0].clear();
		dList[0].clear();
		componentScore[0] = logHiddenPotential[0]
				+ function[0].getLogScoreAndPartialDerivation( seq, start, iList[0], dList[0] );
		iList[1].clear();
		dList[1].clear();
		try
		{
			if( length != 0 )
			{
				componentScore[1] = logHiddenPotential[1]
				                    + function[0].getLogScoreAndPartialDerivation( seq.reverseComplement(), seq.getLength() - start
				                   		- length, iList[1], dList[1] );
			}
			else
			{
				if( start == 0 )
				{
					componentScore[1] = logHiddenPotential[1]
					                    + function[0].getLogScoreAndPartialDerivation( seq.reverseComplement(), 0, iList[1], dList[1] );
				}
				else
				{
					//FIXME
					//componentScore[1] =;
					throw new Exception( "strand scoring for variable length function" );
				}
			}
		}
		catch( Exception doesNotHappen )
		{
			RuntimeException r = new RuntimeException( doesNotHappen.getClass().getName() + ": " + doesNotHappen.getMessage() );
			r.setStackTrace( doesNotHappen.getStackTrace() );
			throw r;
		}
		double logScore = Normalisation.logSumNormalisation( componentScore, 0, 2, componentScore, 0 );
		int j, i = 0;
		for( ; i < logHiddenPotential.length; i++ )
		{
			for( j = 0; j < iList[i].length(); j++ )
			{
				indices.add( iList[i].get( j ) );
				partialDer.add( componentScore[i] * dList[i].get( j ) );
			}
		}
		i = paramRef[2] - paramRef[1];
		for( j = 0; j < i; j++ )
		{
			indices.add( paramRef[1] + j );
			partialDer.add( componentScore[j] - (isNormalized()?hiddenPotential[j]:0) );
		}
		return logScore;
	}

	protected StringBuffer getFurtherInformation()
	{
		StringBuffer erg = new StringBuffer( 100 );
		XMLParser.appendObjectWithTags( erg, forwardPartOfESS, "forwardPartOfESS" );
		XMLParser.appendObjectWithTags( erg, initMethod, "initMethod" );
		return erg;
	}

	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException
	{
		forwardPartOfESS = XMLParser.extractObjectForTags( xml, "forwardPartOfESS", double.class );
		initMethod = XMLParser.extractObjectForTags( xml, "initMethod", InitMethod.class );
	}
	
	//this is needed to enable the invocation of this method by other classes of this package
	protected void init( boolean freeParams )
	{
		super.init( freeParams );
	}
	
	public String toString( NumberFormat nf )
	{
		StringBuffer erg = new StringBuffer( 1500 );
		double d = hiddenPotential[0]+hiddenPotential[1];
		erg.append( "forward: " + nf.format(hiddenPotential[0]/d) + "\n" );
		erg.append( "reverse: " + nf.format(hiddenPotential[1]/d) + "\n\n" );
		erg.append( function[0].toString(nf) );
		return erg.toString();
	}
	
	public boolean modify( int offsetLeft, int offsetRight )
	{
		boolean res = false;
		if( function[0] instanceof Mutable )
		{
			res = ((Mutable) function[0]).modify(offsetLeft, offsetRight) ;
			if( res ) {
				length = function[0].getLength();
				init( freeParams );
				norm = Double.NaN;
			}
		}
		return res;
	}
	
	/**
	 * This method computes the reverse complement distributions for given conditional distributions.
	 * This method is used to determine the context of a motif.
	 * 
	 * @param abc the alphabet
	 * @param condDistr the conditional distribution
	 * 
	 * @return the complement of the conditional distribution that can be used for computing a combing conditional distribution
	 */
	public static double[][][] getReverseComplementDistributions( ComplementableDiscreteAlphabet abc, double[][][] condDistr ) {
		int l = (int)abc.length(), idx = 0, i, o, ord = condDistr.length, h, anz = condDistr[ord-1].length*l;
		double[][][] result = new double[ord][][];
		for( o = 0; o < ord; o++ ) {
			result[o] = new double[condDistr[o].length][l];
		}
		int[] assign = new int[ord];
		double joint;
		for( ; idx < anz; idx++ ) {
			
			//compute assignment
			i = idx;
			for( o = ord-1; o >= 0; o-- ) {
				assign[o] = i % l;
				i /= l;
			}

			// compute joint distribution
			h = 0;
			joint = 1;
			for( o = 0; o < ord; o++ ) {
				joint *= condDistr[o][h][assign[o]];
				h = (h*l) + assign[o];
			}
			
			//complement
			for( o = 0; o < ord; o++ ) {
				assign[o] = abc.getComplementaryCode( assign[o] );
			}
			
			// set reverted
			h = 0;
			for( o = 0; o < ord; o++ ) {
				result[o][h][assign[ord-1-o]] += joint;
				h = (h*l) + assign[ord-1-o];
			}
		}
		
		// marginalize
		for( o = 1; o < ord; o++ ) {
			for( h = 0; h < result[o].length; h++ ) {
				joint = 0;
				for( i = 0; i < l; i++ ) {
					joint += result[o][h][i]; 
				}
				for( i = 0; i < l; i++ ) {
					result[o][h][i] /= joint; 
				}
			}
		}
		
		return result;
	}

	/**
	 * This method returns the preferred {@link Strand} for a given subsequence.
	 * 
	 * @param seq the sequence
	 * @param startPos the start position
	 * 
	 * @return the {@link Strand} of this subsequence
	 * 
	 * @see AbstractMixtureDiffSM#getIndexOfMaximalComponentFor(Sequence, int)
	 */
	public Strand getStrand( Sequence seq, int startPos )
	{
		return getIndexOfMaximalComponentFor( seq, startPos )==0?Strand.FORWARD:Strand.REVERSE;
	}
	
	/**
	 * Check whether a {@link DifferentiableStatisticalModel} is a {@link StrandDiffSM}.
	 * 
	 * @param nsf the original {@link DifferentiableStatisticalModel}
	 * 
	 * @return <code>true</code> if the {@link DifferentiableStatisticalModel} is a {@link StrandDiffSM}
	 */
	public static boolean isStrandModel( DifferentiableStatisticalModel nsf )
	{
		if( nsf instanceof NormalizedDiffSM )
		{
			return ((NormalizedDiffSM) nsf).isStrandModel();
		}
		else
		{
			return nsf instanceof StrandDiffSM;
		}
	}
}
