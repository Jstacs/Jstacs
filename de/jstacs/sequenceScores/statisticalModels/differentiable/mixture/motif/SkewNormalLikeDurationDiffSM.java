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

package de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.motif;

import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.Map.Entry;

import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.RandomNumberGenerator;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * This class implements a skew normal like discrete truncated distribution.
 * 
 * @author Jens Keilwagen
 */
public class SkewNormalLikeDurationDiffSM extends DurationDiffSM
{
	private boolean trainMean, trainPrecision, trainSkew;

	private double par0, par1, par2, hyperMeanMean, hyperMeanStdev, hyperPrec1, hyperPrec2, hyperSkewMean, hyperSkewStdev;

	private double priorC, partDerMu, mu, sigma, prec, logNorm, partDerLogNormPar0, partDerLogNormPar1, partDerLogNormPar2;
	private double[] logScore, densDivCDF;
	private int starts;

	
	private SkewNormalLikeDurationDiffSM( int min, int max, double ess, boolean trainMean, double param0, boolean trainPrecision, double param1, boolean trainSkew, double param2, int starts )
	{
		super( min, max, ess );
		setParameters( param0, param1, param2 );
		this.trainMean = trainMean;
		this.trainPrecision = trainPrecision;
		this.trainSkew = trainSkew;
		if( starts < 1 ) {
			throw new IllegalArgumentException( "The number of starts has to be positive." );
		}
		this.starts = starts;
	}
	
	/**
	 * This is the main constructor if the parameters are fixed.
	 * 
	 * @param min the minimal value
	 * @param max the maximal value
	 * @param param0 the fixed parameter value for the first parameter (mean)
	 * @param param1 the fixed parameter value for the second parameter (precision)
	 * @param param2 the fixed parameter value for the third parameter (skew)
	 */
	public SkewNormalLikeDurationDiffSM( int min, int max, double param0, double param1, double param2 )
	{
		this( min, max, 0, false, param0, false, param1, false, param2, 1 );
	}
	
	/**
	 * This is the constructor that allows the most flexible handling of the parameters.
	 * 
	 * @param min the minimal value
	 * @param max the maximal value
	 * @param trainMean a switch whether to optimize the first parameter
	 * @param hyperMeanMean the mean hyper parameter for the first parameter 
	 * @param hyperMeanSigma the standard deviation hyper parameter for the first parameter
	 * @param trainPrecision a switch whether to optimize the second parameter
	 * @param hyperPrec1 the first hyper parameter for the precision (first parameter of the transformed gamma density);
	 * 		this is value is used to determine the ess: <code>hyperPrec1 = 0.5*ess</code>
	 * @param hyperPrec2 the second hyper parameter for the precision (second parameter of the transformed gamma density)
	 * @param trainSkew a switch whether to optimize the third parameter
	 * @param hyperSkewMean the mean hyper parameter for the third parameter
	 * @param hyperSkewStdev the standard deviation hyper parameter for the third parameter
	 * @param starts the number of recommended starts
	 */
	//ess =^= 2*hyperPrec1
	public SkewNormalLikeDurationDiffSM( int min, int max, boolean trainMean, double hyperMeanMean, double hyperMeanSigma,
			boolean trainPrecision, double hyperPrec1, double hyperPrec2, boolean trainSkew, double hyperSkewMean, double hyperSkewStdev, int starts )
	{
		this( min, max, 2*hyperPrec1, trainMean, 0, trainPrecision, -2 * Math.log((max-min)/4d), trainSkew, 0, starts );
		if( ess > 0 ) {
			if( hyperMeanSigma <= 0 )
			{
				throw new IllegalArgumentException( "The prior of the mean parameter is wrongly specified. (check the second parameter: " + hyperMeanSigma + ")" );
			}
			this.hyperMeanMean = hyperMeanMean;
			this.hyperMeanStdev = hyperMeanSigma;
			
			if( hyperPrec1 <= 0 || hyperPrec2 <= 0 )
			{
				throw new IllegalArgumentException( "The prior of the precision parameter is wrongly specified. (" + hyperPrec1 + ", " + hyperPrec2 + ")" );
			}
			this.hyperPrec1 = hyperPrec1;
			this.hyperPrec2 = hyperPrec2;
			
			if( hyperSkewStdev <= 0 )
			{
				throw new IllegalArgumentException( "The prior of the skew parameter is wrongly specified. (check the second parameter: " + hyperSkewStdev + ")" );
			}
			this.hyperSkewMean = hyperSkewMean;
			this.hyperSkewStdev = hyperSkewStdev;
			precomputePriorConstants();
		}
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link SkewNormalLikeDurationDiffSM} out of a {@link StringBuffer}.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public SkewNormalLikeDurationDiffSM( StringBuffer source ) throws NonParsableException
	{
		super( source );
	}

	public SkewNormalLikeDurationDiffSM clone() throws CloneNotSupportedException {
		SkewNormalLikeDurationDiffSM clone = (SkewNormalLikeDurationDiffSM) super.clone();
		if( logScore != null ) {
			clone.logScore = logScore.clone();
			clone.densDivCDF = densDivCDF.clone();
		}
		return clone;
	}
	
	private static RandomNumberGenerator randNumGen = new RandomNumberGenerator();
	
	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception
	{
		if( data[index].getAlphabetContainer().checkConsistency( alphabets ) ) {
		
			double w = 1;
			int i = 0;
			Integer val;
			double[] weightForVal;
			Hashtable<Integer, double[]> hash = new Hashtable<Integer, double[]>();
			DiscreteAlphabet abc = (DiscreteAlphabet) alphabets.getAlphabetAt( 0 );
			for( ; i < data[index].getNumberOfElements(); i++ ) {
				val = new Integer( abc.getSymbolAt( data[index].getElementAt(i).discreteVal(0) ) );
				if( weights != null && weights[index] != null ) {
					w = weights[index][i];
				}
				weightForVal = hash.get( val );
				if( weightForVal == null ) {
					hash.put( val, new double[]{w} );
				} else {
					weightForVal[0] += w;
				}
			}
			
			Set<Entry<Integer,double[]>> s= hash.entrySet();
			Iterator<Entry<Integer, double[]>> it = s.iterator();
			int[] len = new int[s.size()];
			double[] lenWeights = new double[len.length];
			Entry<Integer, double[]> current;
			i = 0;
			while( it.hasNext() ) {
				current = it.next();
				len[i] = current.getKey();
				lenWeights[i++] = current.getValue()[0];
			}
			
			adjust( len, lenWeights );
		} else {
			System.out.println( "Warning: Try to initialize " + getClass().getName() + " with data over another AlphabetContainer." );
			initializeFunctionRandomly( freeParams );
		}
	}
	
	public void adjust( int[] length, double[] weight )
	{
		double mu = hyperMeanMean, sum = 0, precision = 0, c;
		int i = 0;
		for( ; i < length.length; i++ )
		{
			if( weight[i] > 0 )
			{
				mu += weight[i]*length[i];
				sum += weight[i];
			} else if( Double.isNaN( weight[i] ) ) {
				throw new IllegalArgumentException( "Check the " + i + "-th weight (for length " + length[i] + ")" );
			}
		}
		mu /= (sum+1);
		
		for( i = 0; i < length.length; i++ )
		{
			c = length[i] - mu;
			precision += weight[i]*c*c;
		}
		precision = ( 0.5*sum + hyperPrec1 ) / ( 0.5*precision + hyperPrec2 );
		
		//System.out.println( "adjust gauss (" + sum + "): " + mu + "\t" + precision + "\t(" + (hyperPrec1/hyperPrec2) +")");
	
		//compute parameters			
		precision = Math.log( precision );
		
		mu = (mu-min)/delta;//TODO
		mu = Math.log( mu/(1-mu) );
		
		setParameters( new double[]{mu, precision, 0}, 0 );
		
		//System.out.print( this + "\t" );
	}
	
	
	private final static double V = 6.72;
	
	public void initializeFunctionRandomly( boolean freeParams ) throws Exception
	{
		double[] init = new double[getNumberOfParameters()];
		int i = 0;
		double drawn;
		if( trainMean ) {
			//draw uniform
			init[i] = r.nextDouble();
			
			init[i] = V*init[i]-V/2d;
			i++;
		}
		if( trainPrecision ) {
			//draw a big precision in some random sense
			/*
			drawn = max / 4d;
			drawn += drawn * this.r.nextDouble();
			init[i] = -Math.log( 2 * drawn * drawn );
			*/
			//draw from the (truncated) prior
			
			double d = delta*delta, alpha = ess>0? hyperPrec1 : 1, beta = ess>0 ? hyperPrec2 : 10000 ;
			do
			{
				drawn = randNumGen.nextGamma(alpha,1d/beta);
			}while( drawn > 16d/d );//TODO
			init[i] = Math.log( drawn );
			//System.out.println( drawn + " " + Math.exp( init[i] ) + " " + init[i] );
			i++;
		}
		if( trainSkew ) {
			init[i] = hyperSkewMean + r.nextGaussian()*hyperSkewStdev*hyperSkewStdev;
		}

		setParameters( init, 0 );
		//System.out.println( this );
	}

	protected void fromXML( StringBuffer rep ) throws NonParsableException
	{
		StringBuffer xml = XMLParser.extractForTag( rep, getInstanceName() );
		super.fromXML(xml);
		trainMean = XMLParser.extractObjectForTags( xml, "trainMean", boolean.class );
		trainPrecision = XMLParser.extractObjectForTags( xml, "trainPrecision", boolean.class );
		trainSkew = XMLParser.extractObjectForTags( xml, "trainSkew", boolean.class );
		setParameters( XMLParser.extractObjectForTags( xml, "par0", double.class ), XMLParser.extractObjectForTags( xml, "par1", double.class ), XMLParser.extractObjectForTags( xml, "par2", double.class ) );
		hyperMeanMean = XMLParser.extractObjectForTags( xml, "hyperMeanMean", double.class );
		try{//
		hyperMeanStdev = XMLParser.extractObjectForTags( xml, "hyperMeanStdev", double.class );
		} catch( NonParsableException n ) {
			hyperMeanStdev = 250;
		}
		hyperPrec1 = XMLParser.extractObjectForTags( xml, "hyperPrec1", double.class );
		hyperPrec2 = XMLParser.extractObjectForTags( xml, "hyperPrec2", double.class );
		hyperSkewMean = XMLParser.extractObjectForTags( xml, "hyperSkewMean", double.class );
		hyperSkewStdev = XMLParser.extractObjectForTags( xml, "hyperSkewStdev", double.class );
		precomputePriorConstants();
		
		starts = XMLParser.extractObjectForTags( xml, "starts", int.class );
	}

	public String getInstanceName()
	{
		return getClass().getSimpleName();
	}

	public double[] getCurrentParameterValues() throws Exception
	{
		double[] init = new double[getNumberOfParameters()];
		int i = 0;
		if( trainMean )
		{
			init[i++] = par0;
		}
		if( trainPrecision ) {
			init[i++] = par1;
		}
		if( trainSkew )
		{
			init[i] = par2;
		}
		return init;
	}
	
	public double getLogScore( int... values )
	{
		return logScore[values[0]-min] - logNorm;
	}

	public double getLogScoreAndPartialDerivation( IntList indices, DoubleList partialDer, int... values )
	{
		double z = (values[0] - mu)/sigma, h = z + densDivCDF[values[0]-min] * -par2;
		int i = 0;
		if( trainMean ) {
			indices.add( i++ );
			partialDer.add( -partDerLogNormPar0 + partDerMu/sigma*h );
		}
		if( trainPrecision ) {
			indices.add( i++ );
			partialDer.add( -partDerLogNormPar1 + 0.5 * -z * h );
		}
		if( trainSkew ) {
			indices.add( i++ );
			partialDer.add( -partDerLogNormPar2 + densDivCDF[values[0]-min] * z );
		}
		return logScore[values[0]-min] - logNorm;
	}

	public int getNumberOfParameters()
	{
		return (trainMean ? 1: 0) + (trainPrecision ? 1 : 0) + (trainSkew ? 1: 0);
	}

	public void setParameters( double[] params, int start )
	{
		setParameters( trainMean ? params[start] : par0,
		               trainPrecision ? params[start + (trainMean ? 1 : 0)] : par1,
		               trainSkew ? params[start + (trainMean ? 1 : 0) + (trainPrecision ? 1 : 0)] : par2 );
	}
	
	/**
	 * this method can be used to set the parameters even if the parameters are not allowed to be optimized.
	 *  
	 * @param par0 the first parameter (for the mean or maximum)
	 * @param par1 the second parameter (for the precision)
	 * @param par2 the third parameter (for the skew)
	 */
	public void setParameters( double par0, double par1, double par2 )
	{
		double z, zSq, expCurrent, partDerPhiPart, phi, diff;
		
		this.par0 = par0;
		expCurrent = Math.exp(par0);
		mu = min + delta * (0.01*par0 + expCurrent/(1+expCurrent));
		partDerMu = delta * (0.01 + expCurrent/((1+expCurrent)*(1+expCurrent)));
		
		this.par1 = par1;
		prec = Math.exp( par1 );
		sigma = 1d/Math.sqrt( prec ); 
		
		this.par2 = par2;
		
		if( logScore == null || logScore.length != delta+1 ) {
			logScore = new double[delta+1];
			densDivCDF = new double[delta+1];
		}
		
		partDerLogNormPar0 = partDerLogNormPar1 = partDerLogNormPar2 = 0;
		for( int j = 0; j < logScore.length; j++ )
		{
			diff = (min + j - mu);
			z = diff/sigma;
			zSq = prec*diff*diff;
			
			partDerPhiPart = ONE_DIV_BY_SQRT_OF_2_TIMES_PI * Math.exp( -0.5 * par2 * par2 * zSq );
			densDivCDF[j] = par2==0 ? Math.log(0.5) : CDFOfNormal.getLogCDF( par2*z );
			phi = Math.exp( densDivCDF[j] );
			
			logScore[j] = -0.5 * zSq + densDivCDF[j];
		
			expCurrent = Math.exp( -0.5 * zSq );
			
			/*
			partDerLogNormPar0 += expCurrent * ( z * phi + partDerPhiPart * -par2 );
			partDerLogNormPar1 += expCurrent * 0.5 * z * ( -z * phi + partDerPhiPart * par2 );
			partDerLogNormPar2 += expCurrent * partDerPhiPart * z;
			*/
			
			partDerLogNormPar2 += expCurrent * partDerPhiPart * z;
			
			expCurrent *= ( z * phi + partDerPhiPart * -par2 );
			partDerLogNormPar0 += expCurrent;
			partDerLogNormPar1 -= expCurrent * z;
			
			densDivCDF[j] = ONE_DIV_BY_SQRT_OF_2_TIMES_PI * Math.exp( -0.5*par2*par2*zSq - densDivCDF[j] );
		}	
		logNorm = Normalisation.getLogSum( logScore );
		
		expCurrent = Math.exp( logNorm );
		
		partDerLogNormPar0 = partDerLogNormPar0 * partDerMu/sigma / expCurrent;
		partDerLogNormPar1 = partDerLogNormPar1 * 0.5 / expCurrent;
		partDerLogNormPar2 /= expCurrent;
	}
	
	private static final double ONE_DIV_BY_SQRT_OF_2_TIMES_PI = 1d / Math.sqrt( 2 * Math.PI );

	public StringBuffer toXML()
	{
		StringBuffer xml = super.toXML();
		XMLParser.appendObjectWithTags( xml, trainMean, "trainMean" );
		XMLParser.appendObjectWithTags( xml, trainPrecision, "trainPrecision" );
		XMLParser.appendObjectWithTags( xml, trainSkew, "trainSkew" );
		XMLParser.appendObjectWithTags( xml, par0, "par0" );
		XMLParser.appendObjectWithTags( xml, par1, "par1" );
		XMLParser.appendObjectWithTags( xml, par2, "par2" );
		XMLParser.appendObjectWithTags( xml, hyperMeanMean, "hyperMeanMean" );
		XMLParser.appendObjectWithTags( xml, hyperMeanStdev, "hyperMeanStdev" );
		XMLParser.appendObjectWithTags( xml, hyperPrec1, "hyperPrec1" );
		XMLParser.appendObjectWithTags( xml, hyperPrec2, "hyperPrec2" );
		XMLParser.appendObjectWithTags( xml, hyperSkewMean, "hyperSkewMean" );
		XMLParser.appendObjectWithTags( xml, hyperSkewStdev, "hyperSkewStdev" );
		XMLParser.appendObjectWithTags( xml, starts, "starts" );
		XMLParser.addTags( xml, getInstanceName() );
		return xml;
	}

	protected String getRNotation( String distributionName )
	{
		return "l = " + min + ":" + max + "; " +
			distributionName + " = exp( -0.5 * (l -" + mu + ")^2/" + sigma + "^2 - " + logNorm + " ) * pnorm(" + par2 +"*(l-" + mu +")/" + sigma + ");";
	}

	public double getLogPriorTerm()
	{
		double val = priorC, h;
		if( ess > 0 ) {
			if( trainMean ) {
				h = (mu-hyperMeanMean)/hyperMeanStdev;
				val -= 0.5 * h * h; 
				h = Math.exp( par0 );
				val += Math.log( 0.01 + h/((1+h)*(1+h)) );
			}
			if( trainPrecision ) {
				val += hyperPrec1 * par1 - prec * hyperPrec2;
			}
			if( trainSkew ) {
				h = (par2-hyperSkewMean)/hyperSkewStdev;
				val -= 0.5 * h * h; 
			}
		}
		return val;
	}

	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception
	{
		if( ess > 0 ) {
			if( trainMean ) {
				double h = Math.exp( par0 );
				double h1 = 1+h;
				grad[start++] += -(mu-hyperMeanMean)/(hyperMeanStdev*hyperMeanStdev) * partDerMu + h*(1-h) / h1 / (0.01*h1*h1 + h);
			}
			if( trainPrecision ) {
				grad[start++] += hyperPrec1 - prec * hyperPrec2;
			}
			if( trainSkew ) {
				grad[start++] += -(par2-hyperSkewMean)/(hyperSkewStdev*hyperSkewStdev);
			}
		}
	}

	public boolean isInitialized()
	{
		return true;
	}
	
	public boolean isNormalized()
	{
		return true;
	}

	public void initializeUniformly()
	{
		setParameters( new double[]{0, Double.NEGATIVE_INFINITY, 0}, 0 );
	}
	
	public void modify( int delta ) {
		super.modify( delta );
		precomputePriorConstants();
		setParameters( par0, par1, par2 );
	}
	
	private void precomputePriorConstants() {
		priorC = 0;
		if( trainMean ) {
			priorC += -Math.log( Math.sqrt(2*Math.PI) * hyperMeanStdev ) + Math.log(delta);
		}
		if( trainPrecision ) {
			priorC += hyperPrec1 * Math.log(hyperPrec2) - Gamma.logOfGamma( hyperPrec1 );
		}
		if( trainSkew ) {
			priorC -= Math.log( Math.sqrt( 2 * Math.PI ) * hyperSkewStdev );
		}
	}
	
	public int getNumberOfRecommendedStarts() {
		return starts;
	}
}
