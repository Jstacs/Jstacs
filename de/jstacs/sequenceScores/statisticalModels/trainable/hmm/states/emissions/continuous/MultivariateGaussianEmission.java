package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.continuous;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import javax.naming.OperationNotSupportedException;

import Jama.Matrix;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;

/**
 * Multivariate Gaussian emission density for a Hidden Markov Model. 
 * 
 * @author Michael Seifert, Jens Keilwagen
 *
 */
public class MultivariateGaussianEmission implements Emission
{
	/**
	 * Dimensions
	 */
	int dim;
	
	/**
	 * Initial mean vector.
	 */
	double[] initialMean;
	
	/**
	 * Mean vector.
	 */
	double[] mean;
	
	/**
	 * Initial standard deviations.
	 */
	double[] initialSds;
	
	/**
	 * Standard deviations.
	 */
	double[] sds;
	
	/**
	 * Initial correlations between pairwise dimensions.
	 */
	double[][] initialCorrelation;
	
	/**
	 * Correlations between pairwise dimensions.
	 */
	double[][] correlation;
	
	/**
	 * Inverse covariance matrix: internal parameterization.
	 */
	private Matrix inverseCov;
	

	//Prior hyper-parameters
	
	/**
	 * A priori mean (nhi).
	 */
	private double[] aprioriMean;
	
	/**
	 * Scaling factor for apriori mean (epsilon).
	 */
	private double scaleMean;
	
	/**
	 * Shape of inverse-wishart prior (r).
	 */
	private double shapeSd;
	
	/**
	 * Scale matrix of inverse-wishart prior (Omega).
	 */
	private double[][] scaleSd;

	//Count statistics and helpers
	
	/**
	 * Contains the emission sequences and corresponding gammas (state-posteriors) required for the estimation of the standard deviation.
	 */
	//protected Hashtable<Sequence,double []> gammas;
	protected HashMap<Sequence,double []> gammas;
	
	/**
	 * Sum of gammas for estimation of mean vector and covariance matrix.
	 */
	private double sumOfGammas;
	
	/**
	 * Sum of gamma_t(i)*emission_t
	 * - sumOfGammaWeightedEmissions[ 0 ] -> first dimension
	 * - sumOfGammaWeightedEmissions[ 1 ] -> second dimension
	 * - ...
	 * - sumOfGammaWeightedEmissions[ this.dim-1 ] -> this.dim dimension
	 */
	private double[] sumOfGammaWeightedEmissions;
	
	
	private AlphabetContainer con;
	
	private static final String TAG = "MultivariateGaussianEmission";
	
	double[] emission;

	
	/**
	 * Creates a Multivariate Gaussian emission density.
	 * 
	 * @param mean mean vector
	 * @param sds standard deviation vector
	 * @param correlation correlation matrix
	 * @param scaleMean scale of a priori mean
	 * @param aprioriMean a priori mean
	 * @param shapeSd scale for scaleSd (a covariance matrix)
	 * @param scaleSd a priori covariance matrix
	 */
	public MultivariateGaussianEmission( double[] mean, double[] sds, double[][] correlation, double scaleMean, double aprioriMean[], double shapeSd, double[][] scaleSd )
	{
		this.dim = mean.length;
		this.emission = new double[ this.dim ];
		
		//Initially specified parameters
		this.initialMean = mean.clone();
		this.initialSds  = sds.clone();
		this.initialCorrelation = correlation.clone(); 
		
		//Working parameters
		this.mean = mean.clone();
		this.sds  = sds.clone();
		this.correlation = correlation.clone();
		this.inverseCov  = getInverseCovarianceMatrix();
		
		//Prior parameters
		this.scaleMean   = scaleMean;
		this.aprioriMean = aprioriMean.clone();
		this.shapeSd     = shapeSd;
		this.scaleSd     = scaleSd.clone();
		
		//Emission sequences and corresponding gammas for estimation of sd
		this.gammas = new HashMap<Sequence, double[]>();
		
		//Alphabet container
		this.con = new AlphabetContainer( new ContinuousAlphabet() );
		
	}


	/**
	 * Creates a {@link MultivariateGaussianEmission} from its XML representation.
	 * @param xml the XML representation.
	 * @throws NonParsableException if the XML representation could not be parsed
	 */
	public MultivariateGaussianEmission( StringBuffer xml ) throws NonParsableException
	{
		this.fromXML( xml );
		//this.gammas = new Hashtable<Sequence, double[]>();
		this.gammas = new HashMap<Sequence, double[]>();
		this.resetStatistic();	
	}

	
	
	public void addToStatistic(boolean forward, int startPos, int endPos, double weight, Sequence seq) throws OperationNotSupportedException
	{
		if(!forward){
			throw new OperationNotSupportedException();
		}

		
		//For HashMap
		boolean contained = this.gammas.containsKey( seq );
		if( !contained )
		{
			this.gammas.put( seq, new double[ seq.getLength() ] );
		}		
		double[] weights = this.gammas.get( seq );		

		
		//startPos == endPos (at least in my models)
		for( int pos = startPos; pos <= endPos; pos++ )
		{
			//Save weight
			weights[ pos ] += weight;

			seq.fillContainer( emission, pos);			
			
			//Sum of gammas
			this.sumOfGammas += weight;
			
			//Sum of gamma weighted emissions
			for( int d = 0; d < this.dim; d++ )
			{
				this.sumOfGammaWeightedEmissions[ d ] += weight * emission[ d ];
			}		
		}		

	}
	
	@Override
	public void joinStatistics(Emission... emissions) {
		for(int i=0;i<emissions.length;i++){
			if(emissions[i] != this){
				MultivariateGaussianEmission c = (MultivariateGaussianEmission) emissions[i];
				
				sumOfGammas += c.sumOfGammas;
				for( int d = 0; d < this.dim; d++ )
				{
					this.sumOfGammaWeightedEmissions[ d ] += c.sumOfGammaWeightedEmissions[ d ];
				}
				
				/*simple way*/
				//this.gammas.putAll(c.gammas);
				/*complex secure way*/			
				Iterator<Entry<Sequence,double[]>> it = c.gammas.entrySet().iterator();
				while( it.hasNext() ) {
					Entry<Sequence,double[]> e = it.next();
					double[] v = gammas.get(e.getKey());
					if( v == null ) {
						gammas.put( e.getKey(), e.getValue().clone() );
					} else {
						double[] w = e.getValue();
						for( int j = 0; j < w.length; j++ ) {
							v[j] += w[j];
						}
						gammas.put( e.getKey(), v );
					}
				}
				/**/
			}
		}
		
		//copy
		for(int i=0;i<emissions.length;i++){
			if(emissions[i] != this){
				MultivariateGaussianEmission c = (MultivariateGaussianEmission) emissions[i];
				
				c.sumOfGammas = sumOfGammas;
				for( int d = 0; d < this.dim; d++ )
				{
					c.sumOfGammaWeightedEmissions[ d ] = this.sumOfGammaWeightedEmissions[ d ];
				}
				
				c.gammas.clear();
				/*simple way*/
				//c.gammas.putAll( gammas );
				/*complex secure way*/
				/*
				Iterator<Entry<Sequence,double[]>> it = c.gammas.entrySet().iterator();
				while( it.hasNext() ) {
					Entry<Sequence,double[]> e = it.next();
					c.gammas.put(e.getKey(), e.getValue().clone());
				}/**/
				/*secure, efficient way*/
				c.resetGammas();
				Iterator<Entry<Sequence,double[]>> it = gammas.entrySet().iterator();
				while( it.hasNext() ) {
					Entry<Sequence,double[]> e = it.next();
					double[] v = c.gammas.get(e.getKey());
					if( v == null ) {
						c.gammas.put(e.getKey(), e.getValue().clone());
					} else {
						System.arraycopy(e.getValue(), 0, v, 0, v.length);
					}
				}
			}
		}
	}
	
	public void estimateFromStatistic()
	{
		//MAP Estimate mean vector
		for( int d = 0; d < this.dim; d++ )
		{
			this.mean[ d ] = ( this.sumOfGammaWeightedEmissions[ d ] + this.aprioriMean[ d ] * this.scaleMean ) / ( this.sumOfGammas + this.scaleMean );
		}

		//Estimate covariance matrix in dependency of the new mean
		Iterator<Entry<Sequence, double[]>> keyIterator = this.gammas.entrySet().iterator();
		
		
		//MAP: add pseudocounts
		//cov_numerator

		Matrix cov_numerator   = new Matrix( this.scaleSd );		
		Matrix dummy = new Matrix( new double[][] { this.mean } );
		dummy = dummy.minus( new Matrix( new double[][]{ this.aprioriMean } ) );
		cov_numerator = cov_numerator.plus( dummy.transpose().times( dummy ).times( this.scaleMean ) );
		
		//cov_denominator
		double cov_denominator = this.sumOfGammas + this.shapeSd - this.dim;

		do
		{
			Entry<Sequence, double[]> entry = keyIterator.next();			
			Sequence seq = entry.getKey();
			double[] weights = entry.getValue();
			
			int T = weights.length;
			for( int t = 0; t < T; t++ )
			{
				seq.fillContainer( emission, t );
				
				dummy = new Matrix( new double[][] { emission } );
				dummy = dummy.minus( new Matrix( new double[][]{ this.mean } ) );
				cov_numerator = cov_numerator.plus( dummy.transpose().times( dummy ).times( weights[ t ] ) );		
			}
			
		} while( keyIterator.hasNext() );

		
		//Covariance matrix
		Matrix cov = cov_numerator.times( 1 / cov_denominator );

		for( int d = 0; d < this.dim; d++ )
		{
			this.sds[ d ] = Math.sqrt( cov.get( d, d ) );			
		}
		this.correlation = getCorrelations( cov, this.sds );
		
		//Inverse covariance matrix internally used parameterization.				
		this.inverseCov = cov.inverse();
	}

	public AlphabetContainer getAlphabetContainer()
	{
		return this.con;
	}

	@Override
	public double getLogPriorTerm()
	{
		double res = 0;
		

		res += (this.shapeSd - this.dim) / 2 * Math.log( this.inverseCov.det() );
		
		Matrix dummy = new Matrix( new double[][] { this.mean } );
		dummy = dummy.minus( new Matrix( new double[][] { this.aprioriMean } ) );
		Matrix dummyRes = dummy.times( this.inverseCov ).times( dummy.transpose() );
		
		res -= ( this.scaleMean / 2 ) * dummyRes.get( 0, 0 );
		
		dummy = new Matrix( this.scaleSd );
		dummyRes = dummy.times( this.inverseCov );
		
		res -= 0.5 * dummyRes.trace();
		
		return res;
	}

	public double getLogProbFor(boolean forward, int startPos, int endPos, Sequence seq) throws OperationNotSupportedException
	{
		double res = 0;

		for(int pos = startPos; pos <= endPos;pos++)
		{			
			//Get corresponding emission vector
			seq.fillContainer( emission, pos);
			
			res -= Math.log( Math.sqrt( Math.pow( 2 * Math.PI, this.dim ) ) ); //pre-computed would be faster
			res += 0.5 * Math.log( this.inverseCov.det() ); //pre-computed would be faster
			Matrix dummy = new Matrix( new double[][] { emission } );
			dummy = dummy.minus( new Matrix( new double[][] { this.mean } ) );
			dummy = dummy.times( this.inverseCov ).times( dummy.transpose() );
			res -= 0.5 * dummy.get( 0, 0 );
		}		

		return res;
	}

	@Override
	public String getNodeLabel(double weight, String name, NumberFormat nf) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getNodeShape(boolean forward) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void initializeFunctionRandomly()
	{
		//Working parameters
		this.mean = this.initialMean.clone();
		this.sds  = this.initialSds.clone();
		this.correlation = this.initialCorrelation.clone();
		this.inverseCov  = getInverseCovarianceMatrix();	
	}

	public void resetStatistic()
	{
		//Reset count statistics
		this.sumOfGammas = 0;
		if( this.sumOfGammaWeightedEmissions == null )
		{
			this.sumOfGammaWeightedEmissions = new double[ this.dim ];
		}
		else
		{
			Arrays.fill(this.sumOfGammaWeightedEmissions, 0 );
		}
		resetGammas();
	}
	
	private void resetGammas() {
		//For HashMap
		if( !this.gammas.isEmpty() )
		{
			Iterator<Entry<Sequence, double[]>> keyIterator = this.gammas.entrySet().iterator();
			
			do
			{
				Entry<Sequence, double[]> entry = keyIterator.next();			
				Sequence seq = entry.getKey();
				double[] weights = entry.getValue();
				Arrays.fill( weights, 0 );		
				
			} while( keyIterator.hasNext() );
		}
	
	}

	@Override
	public StringBuffer toXML()
	{
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags( buf, con, "alphabet" );
		XMLParser.appendObjectWithTags( buf, this.initialMean, "initialMean" );
		XMLParser.appendObjectWithTags( buf, this.mean, "mean" );
		XMLParser.appendObjectWithTags( buf, this.initialSds, "initialSds" );
		XMLParser.appendObjectWithTags( buf, this.sds, "sds" );
		XMLParser.appendObjectWithTags( buf, this.initialCorrelation, "initialCorrelation" );
		XMLParser.appendObjectWithTags( buf, this.correlation, "correlation" );
		XMLParser.appendObjectWithTags( buf, this.aprioriMean, "aprioriMean" );
		XMLParser.appendObjectWithTags( buf, this.scaleMean, "scaleMean" );
		XMLParser.appendObjectWithTags( buf, this.shapeSd, "shapeSd" );
		XMLParser.appendObjectWithTags( buf, this.scaleSd, "scaleSd" );
		XMLParser.addTags( buf, TAG );
		
		return buf;
	}
	
	protected void fromXML( StringBuffer xml ) throws NonParsableException
	{
		xml = XMLParser.extractForTag( xml, TAG );
				
		this.con = XMLParser.extractObjectForTags( xml, "alphabet", AlphabetContainer.class );
		this.initialMean   = (double[]) XMLParser.extractObjectForTags( xml, "initialMean" );
		this.mean          = (double[]) XMLParser.extractObjectForTags( xml, "mean" );
		this.initialSds    = (double[]) XMLParser.extractObjectForTags( xml, "initialSds" );
		this.sds           = (double[]) XMLParser.extractObjectForTags( xml, "sds" );
		this.initialCorrelation  = (double[][]) XMLParser.extractObjectForTags( xml, "initialCorrelation" );
		this.correlation         = (double[][]) XMLParser.extractObjectForTags( xml, "correlation" );		
		this.aprioriMean         = (double[]) XMLParser.extractObjectForTags( xml, "aprioriMean" );
		this.scaleMean           = XMLParser.extractObjectForTags( xml, "scaleMean", double.class );		
		this.shapeSd             = XMLParser.extractObjectForTags( xml, "shapeSd", double.class );
		this.scaleSd             = (double[][]) XMLParser.extractObjectForTags( xml, "scaleSd" );
		
		this.dim = mean.length;
		this.emission = new double[ this.dim ];
		
		//Re-construct inverse-covariance matrix
		this.inverseCov = this.getInverseCovarianceMatrix();
	}
	
	/*
	 * String representation of multivariate Gaussian emission.
	 * 
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission#toString(java.text.NumberFormat)
	 */
	public String toString( NumberFormat nf )
	{
		String res = "- Means  = ";
		for( int i = 0; i < this.dim; i++ )
		{
			res += nf.format(this.mean[ i ]) + "\t";
		}
		res += "\n\n";
		
		for( int i = 0; i < this.dim; i++ )
		{
			res += "- Standard dev. = " + nf.format(this.sds[ i ]) + "\n";
		}
		res += "\n\n";		
		
		for( int i = 0; i < this.dim; i++ )
		{
			for( int j = i + 1; j < this.dim; j++ )
			{
				res += "- Correlation(" + ( i + 1 ) + "," + ( j + 1 ) + ")  = ";
				res += nf.format(this.correlation[ i ][ j ]) + "\n";
			}
		}		
				
		return res;
	}
	
	public MultivariateGaussianEmission clone() throws CloneNotSupportedException
	{
		MultivariateGaussianEmission clone = (MultivariateGaussianEmission) super.clone();
		clone.correlation = ArrayHandler.clone(correlation);
		clone.emission = emission==null?null:emission.clone();
		if( gammas != null ) {
			clone.gammas = new HashMap<Sequence,double[]>();
			Iterator<Entry<Sequence,double[]>> it = gammas.entrySet().iterator();
			while( it.hasNext() ) {
				Entry<Sequence,double[]> e = it.next();
				clone.gammas.put( e.getKey(), e.getValue().clone() );
			}
		} else {
			clone.gammas = null;
		}
		clone.initialCorrelation = ArrayHandler.clone(initialCorrelation);
		clone.initialMean = initialMean==null?null:initialMean.clone();
		clone.initialSds = initialSds==null?null:initialSds.clone();
		clone.inverseCov = inverseCov==null? null : new Matrix( inverseCov.getArray() );
		clone.scaleSd = ArrayHandler.clone(scaleSd);
		clone.sds = sds==null?null:sds.clone();
		clone.sumOfGammaWeightedEmissions = sumOfGammaWeightedEmissions==null?null:sumOfGammaWeightedEmissions.clone();
		return clone;
	}
	
	/**
	 * Computes the inverse-covariance of a multivariate Gaussian.
	 *  
	 * @return inverse-covariance matrix
	 */
	private Matrix getInverseCovarianceMatrix()
	{			
		return getCovarianceMatrix().inverse();
		
	}

	/**
	 * Computes the covariance matrix of a multivariate Gaussian based on given standard deviations and correlations.
	 * 
	 * @return covariance matrix
	 */
	private Matrix getCovarianceMatrix()
	{
		Matrix dummy = new Matrix( this.dim, this.dim, 0 );
		
		//Set variances
		for( int i = 0; i < this.dim; i++ )
		{
			dummy.set( i, i, Math.pow( this.sds[ i ], 2 ) );
		}
		
		//Set covariances
		for( int i = 0; i < this.dim; i++ )
		{
			for( int j = i + 1; j < this.dim; j++ )
			{
				double covVar = this.sds[ i ] * this.sds[ j ] * this.correlation[ i ][ j ];
				dummy.set( i, j, covVar );
				dummy.set( j, i, covVar );
			}
		}
		
		return dummy;
	}
	
	/**
	 * Computes the correlation parameters included in a given covariance matrix and its corresponding standard deviations.
	 * 
	 * @param cov covariance matrix
	 * 
	 * @param sds standard deviations
	 * 
	 * @return correlations
	 */
	private double[][] getCorrelations( Matrix cov, double[] sds )
	{
		double[][] res = new double[ this.dim ][ this.dim ];
		
		for( int i = 0; i < this.dim; i++ )
		{
			for( int j = i + 1; j < this.dim; j++ )
			{
				res[ i ][ j ] = cov.get( i, j ) / sds[ i ] / sds[ j ];
			}
		}
		
		return res;
	}


	@Override
	public void setParameters(Emission t) throws IllegalArgumentException {
		if( !t.getClass().equals( getClass() ) || ((MultivariateGaussianEmission)t).dim != dim ) {//TODO more?
			throw new IllegalArgumentException( "The transitions are not comparable." );
		}
		MultivariateGaussianEmission c = (MultivariateGaussianEmission) t;
		
		System.arraycopy(c.mean, 0, mean, 0, mean.length);
		System.arraycopy(c.sds, 0, sds, 0, sds.length);
		inverseCov = c.inverseCov.copy();
		for( int i = 0; i < correlation.length; i++ ) {
			System.arraycopy(c.correlation[i], 0, correlation[i], 0, correlation[i].length);
		}
		
	}
}
