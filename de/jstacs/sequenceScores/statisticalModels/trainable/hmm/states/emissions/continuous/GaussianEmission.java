package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.continuous;

import java.text.NumberFormat;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.random.RandomNumberGenerator;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * Emission for continuous values following a Gaussian distribution. The Gaussian density is parameterized in terms of mean and (log) precision.
 * The prior of the Gaussian density is a normal-gamma density parameterized in terms of shape and rate.
 * 
 * @author Jan Grau
 *
 */
public class GaussianEmission implements DifferentiableEmission {

	private static RandomNumberGenerator rand = new RandomNumberGenerator();
	
	private double ess;
	private double priorMu;
	private double priorAlpha;
	private double priorBeta;
	
	private double mu;
	private double logPrecision;
	private double precision;
	
	private double mean;
	private double meansq;
	private double n;
	
	private double logNorm;
	
	private boolean transformed;
	
	private int offset;
	
	private AlphabetContainer con;

	/**
	 * Creates a {@link GaussianEmission} which can be used for maximum likelihood.
	 * 
	 * @param con the alphabet of the emissions
	 */
	public GaussianEmission( AlphabetContainer con ) {
		this( con, 0, 0, 0, 0, false );
	}
	
	/**
	 * Creates a {@link GaussianEmission} with normal-gamma prior by directly defining the hyper-parameters of the prior.
	 * 
	 * @param con the alphabet of the emissions
	 * @param ess the equivalent sample size of the normal-gamma prior
	 * @param priorMu the a-priori mean of the normal part of the prior
	 * @param priorAlpha the shape parameter of the gamma part of the prior
	 * @param priorBeta the rate parameter of the gamma part of the prior
	 * @param transformed use the transformed Gaussian density, i.e. exponential precision, ignored for numerical optimization
	 */
	public GaussianEmission( AlphabetContainer con,double ess, double priorMu, double priorAlpha, double priorBeta, boolean transformed ) {
		this.con = con;
		this.ess = ess;
		this.priorMu = priorMu;
		this.priorAlpha = priorAlpha;
		this.priorBeta = priorBeta;
		this.transformed = transformed;
	}
	
	/**
	 * Creates a {@link GaussianEmission} with normal-gamma prior by defining the expected precision and the expected standard deviation of the precision, i.e. via the 
	 * expectation and variance of the gamma part of the normal-gamma prior.
	 * 
	 * @param ess the equivalent sample size of the normal-gamma prior
	 * @param con the alphabet of the emissions
	 * @param priorMu the a-priori mean of the normal part of the prior
	 * @param expectedPrecision the expected value of the precision
	 * @param sdPrecision the expected standard deviation of the precision
	 * @param transformed use the transformed Gaussian density, i.e. exponential precision, ignored for numerical optimization
	 */
	public GaussianEmission( double ess, AlphabetContainer con, double priorMu, double expectedPrecision, double sdPrecision, boolean transformed ) {
		this(con,ess,priorMu,( expectedPrecision/(2.0*sdPrecision*sdPrecision)+Math.sqrt( (expectedPrecision/(2.0*sdPrecision*sdPrecision))*(expectedPrecision/(2.0*sdPrecision*sdPrecision)) + 1.0/(sdPrecision*sdPrecision) ) )*expectedPrecision + 1.0 ,expectedPrecision/(2.0*sdPrecision*sdPrecision)+Math.sqrt( (expectedPrecision/(2.0*sdPrecision*sdPrecision))*(expectedPrecision/(2.0*sdPrecision*sdPrecision)) + 1.0/(sdPrecision*sdPrecision) ), transformed);
	}

	/**
	 * Creates a {@link GaussianEmission} from its XML representation.
	 * @param xml the XML representation.
	 * @throws NonParsableException if the XML representation could not be parsed
	 */
	public GaussianEmission( StringBuffer xml ) throws NonParsableException {
		fromXML( xml );
	}
	
	public GaussianEmission clone() throws CloneNotSupportedException {
		GaussianEmission clone = (GaussianEmission) super.clone();
		return clone;
	}

	
	public void addGradientOfLogPriorTerm( double[] gradient, int offset ) {
		if( ess > 0 ) {
			double val = (mu-priorMu);
			double gradmu = ess*precision*val;
			gradient[this.offset+offset] -= gradmu;
			gradient[this.offset+offset + 1] += 0.5 - 0.5*gradmu*val + priorAlpha - priorBeta*precision;
		}
	}

	
	public void joinStatistics(Emission... emissions){
		for(int i=0;i<emissions.length;i++){
			if(emissions[i] != this){
				mean += ((GaussianEmission)emissions[i]).mean;
				meansq += ((GaussianEmission)emissions[i]).meansq;
				n += ((GaussianEmission)emissions[i]).n;
			}
		}
		for(int i=0;i<emissions.length;i++){
			((GaussianEmission)emissions[i]).mean = this.mean;
			((GaussianEmission)emissions[i]).meansq = this.meansq;
			((GaussianEmission)emissions[i]).n = this.n;
		}
	}
	
	public void addToStatistic( boolean forward, int startPos, int endPos, double weight, Sequence seq ) throws OperationNotSupportedException {
		if(!forward){
			throw new OperationNotSupportedException();
		}
		double w;
		for(;startPos<=endPos;startPos++){
			w = weight * seq.continuousVal( startPos ); 
			mean += w;
			meansq += w*seq.continuousVal( startPos );
			n += weight;
		}
	}

	
	public void estimateFromStatistic() {
		if(ess == 0){
			if( n == 0 ) { //=> mean = meanSq = 0
				//if no data has been seen the distribution degenerates to a uniform distribution
				//in this case the mean can be chosen arbitrary and precision has to be 0
				n = 1;
			}
			mu = mean / n;
			precision = meansq/n - mu*mu;
			if( precision != 0 ) {
				precision = 1.0/precision;
			} else {
				//using L'Hospital we obtain a uniform distribution that can be parameterized by precision=0 
				precision = 0;
			}	
		}else{
			mu = (mean + ess*priorMu)/(n + ess);
			double s = meansq - 2.0*mu*mean + n*mu*mu;
			double s0 = (mu-priorMu)*(mu-priorMu);
			if(transformed){
				precision = (n + 2.0*priorAlpha + 1)/(s + 2.0*priorBeta + ess*s0);
			}else{
				precision = (n + 2.0*priorAlpha - 1)/(s + 2.0*priorBeta + ess*s0);
			}
		}
		logPrecision = Math.log( precision );
		precompute();
	}

	
	public void fillCurrentParameter( double[] params ) {
		params[offset] = mu;        
		params[offset+1] = logPrecision;
	}

	
	public void setParameter( double[] params, int offset ) {
		mu = params[offset];
		logPrecision = params[offset+1];
		precision = Math.exp( logPrecision );
		precompute();
	}

	
	public int setParameterOffset( int offset ) {
		this.offset = offset;
		return offset + 2;
	}

	/**
	 * This method is internally used by the constructor {@link #GaussianEmission(StringBuffer)}.
	 * 
	 * @param xml the {@link StringBuffer} containing the xml representation of an instance
	 * 
	 * @throws NonParsableException if the {@link StringBuffer} is not parsable
	 * 
	 * @see #GaussianEmission(StringBuffer)
	 */
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		con = XMLParser.extractObjectForTags( xml, "alphabet", AlphabetContainer.class );
		ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
		priorMu = XMLParser.extractObjectForTags( xml, "priorMu", double.class );
		priorAlpha = XMLParser.extractObjectForTags( xml, "priorAlpha", double.class );
		priorBeta = XMLParser.extractObjectForTags( xml, "priorBeta", double.class );
		mu = XMLParser.extractObjectForTags( xml, "mu", double.class );
		logPrecision = XMLParser.extractObjectForTags( xml, "logPrecision", double.class );
		precision = Math.exp( logPrecision );
		mean = XMLParser.extractObjectForTags( xml, "mean", double.class );
		meansq = XMLParser.extractObjectForTags( xml, "meansq", double.class );
		transformed = XMLParser.extractObjectForTags( xml, "transformed", boolean.class );
		n = XMLParser.extractObjectForTags( xml, "n", double.class );
		offset = XMLParser.extractObjectForTags( xml, "offset", int.class );
		precompute();
	}

	
	public double getLogPriorTerm() {
		if( ess > 0 ) {
			double val = mu - priorMu;
			
			return 0.5*( Math.log( ess/(2.0*Math.PI) ) + logPrecision  - ess*precision*val*val )
						+ priorAlpha*Math.log( priorBeta ) - Gamma.logOfGamma( priorAlpha ) + priorAlpha*logPrecision - priorBeta*precision;
		} else {
			return 0;
		}
	}

	
	public double getLogProbAndPartialDerivationFor( boolean forward, int startPos, int endPos, IntList indices, DoubleList partDer,
			Sequence seq ) throws OperationNotSupportedException {
		if(!forward){
			throw new OperationNotSupportedException();
		}
		double score = 0;
		double derivmu = 0;
		double derivprc = 0;
		for(;startPos<=endPos;startPos++){
			double val = seq.continuousVal( startPos )-mu;

			derivmu += precision*val;
			derivprc += 0.5*(1.0 - precision*val*val);

			score += logNorm - 0.5*val*val*precision;
		}
		indices.add( offset );
		partDer.add( derivmu );
		indices.add( offset + 1 );
		partDer.add( derivprc );
		return score;
	}

	
	public double getLogProbFor( boolean forward, int startPos, int endPos, Sequence seq ) throws OperationNotSupportedException {
		double score = 0;
		for(;startPos<=endPos;startPos++){
			double val = seq.continuousVal( startPos ) - mu;
			score += logNorm - 0.5*val*val*precision;
		}
		return score;
	}

	
	public void initializeFunctionRandomly() {
		if( ess == 0 ) {
			//TODO ML case?
			precision = rand.nextGamma( 1.0, 1.0 );
			mu = rand.nextGaussian()/(precision) + priorMu;
		} else {
			precision = rand.nextGamma( priorAlpha, 1.0/priorBeta );
			mu = rand.nextGaussian()/(ess*precision) + priorMu;
		}
		logPrecision = Math.log( precision );
		
		precompute();
	}

	/**
	 * This method precomputes some normalization constant.
	 */
	protected void precompute() {
		logNorm = 0.5*(logPrecision - Math.log( 2.0*Math.PI ));
	}

	
	public void resetStatistic() {
		mean = meansq = n = 0;
	}

	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags( buf, con, "alphabet" );
		XMLParser.appendObjectWithTags( buf, ess, "ess" );
		XMLParser.appendObjectWithTags( buf, priorMu, "priorMu" );
		XMLParser.appendObjectWithTags( buf, priorAlpha, "priorAlpha" );
		XMLParser.appendObjectWithTags( buf, priorBeta, "priorBeta" );
		XMLParser.appendObjectWithTags( buf, mu, "mu" );
		XMLParser.appendObjectWithTags( buf, logPrecision, "logPrecision" );
		XMLParser.appendObjectWithTags( buf, mean, "mean" );
		XMLParser.appendObjectWithTags( buf, meansq, "meansq" );
		XMLParser.appendObjectWithTags( buf, n, "n" );
		XMLParser.appendObjectWithTags( buf, transformed, "transformed" );
		XMLParser.appendObjectWithTags( buf, offset, "offset" );
		XMLParser.addTags( buf, getClass().getSimpleName() );
		return buf;
	}

	public AlphabetContainer getAlphabetContainer() {
		return con;
	}
	
	public String toString() {
		return "p = sqrt(" + precision + "/(2*pi)) * exp( -0.5 * " + precision + " * (x - " + mu + ")^2 );\n";
	}

	@Override
	public String getNodeShape( boolean forward ) {
		return "\"box\"";
	}

	@Override
	public String getNodeLabel( double weight, String name, NumberFormat nf ) {
		return "\""+name+"\"";
	}

	@Override
	public void fillSamplingGroups( int parameterOffset, LinkedList<int[]> list ) {
		list.add(new int[]{parameterOffset+offset, parameterOffset+offset+1});
	}

	@Override
	public int getNumberOfParameters() {
		return 2;
	}

	@Override
	public int getSizeOfEventSpace() {
		return 0;
	}
	
	@Override
	public void setParameters(Emission t) throws IllegalArgumentException {
		if( !t.getClass().equals( getClass() ) ) {
			throw new IllegalArgumentException( "The transitions are not comparable." );
		}
		GaussianEmission tt = (GaussianEmission) t;
		mu = tt.mu;
		logPrecision = tt.logPrecision;
		precision = tt.precision;
		logNorm = tt.logNorm;
	}
	
}