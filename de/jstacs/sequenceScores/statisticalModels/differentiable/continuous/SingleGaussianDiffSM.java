package de.jstacs.sequenceScores.statisticalModels.differentiable.continuous;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Random;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.DinucleotideProperty;
import de.jstacs.data.sequences.ArbitrarySequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.random.RandomNumberGenerator;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

public class SingleGaussianDiffSM extends AbstractDifferentiableStatisticalModel {

	private static RandomNumberGenerator rand = new RandomNumberGenerator();
	
	private double ess;
	private double priorMu;
	private double priorAlpha;
	private double priorBeta;
	
	private double mu;
	private double logPrecision;
	private double precision;
	
	private boolean initialized;
	private boolean alwaysInitRandomly;
	
	private double logNorm;
	
	private DinucleotideProperty prop;

	private boolean fixMu;
	
	public SingleGaussianDiffSM(double ess, double priorMu, double priorAlpha, double priorBeta, AlphabetContainer alphabet, boolean alwaysInitRandomly, boolean fixMu, double mu ){
		super(alphabet,1);
		this.fixMu = fixMu;
		this.mu = mu;
		this.ess = ess;
		this.priorMu = priorMu;
		this.priorAlpha = priorAlpha;
		this.priorBeta = priorBeta;
		this.initialized = false;
		this.alwaysInitRandomly = alwaysInitRandomly;
	}
	
	public SingleGaussianDiffSM(AlphabetContainer alphabet, double ess, double priorMu, double expectedPrecision, double sdPrecision, boolean alwaysInitRandomly){
		this(ess,priorMu,( expectedPrecision/(2.0*sdPrecision*sdPrecision)+Math.sqrt( (expectedPrecision/(2.0*sdPrecision*sdPrecision))*(expectedPrecision/(2.0*sdPrecision*sdPrecision)) + 1.0/(sdPrecision*sdPrecision) ) )*expectedPrecision + 1.0 ,expectedPrecision/(2.0*sdPrecision*sdPrecision)+Math.sqrt( (expectedPrecision/(2.0*sdPrecision*sdPrecision))*(expectedPrecision/(2.0*sdPrecision*sdPrecision)) + 1.0/(sdPrecision*sdPrecision) ), alphabet, alwaysInitRandomly, false, 0);
	}
	
	public SingleGaussianDiffSM(StringBuffer buf) throws NonParsableException{
		super(buf);
	}
	
	public SingleGaussianDiffSM clone() throws CloneNotSupportedException{
		return (SingleGaussianDiffSM) super.clone();
	}
	
	
	public DataSet emitDataSet( int numberOfSequences, int... seqLength) throws Exception {
		Sequence[] seqs = new Sequence[numberOfSequences];
		Random r = new Random();
		int l = seqLength[0];
		double sd = Math.sqrt(1/precision);
		for( int i = 0; i < numberOfSequences; i++ ) {
			double[] v;
			if( seqLength.length > 1 ) {
				v = new double[seqLength[i]];
			} else {
				v = new double[l];
			}
			for( int k = 0; k < v.length; k++ ) {
				v[k]=r.nextDouble()*sd+mu;
			}
			seqs[i] = new ArbitrarySequence(alphabets, v);
		}
		return new DataSet( "sampled from " + getInstanceName(), seqs );
	}
	
	public void setDinucleotideProperty(DinucleotideProperty prop){
		this.prop = prop;
	}
	
	@Override
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		length = 1;
		xml = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		alphabets = XMLParser.extractObjectForTags( xml, "alphabet", AlphabetContainer.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
		priorMu = XMLParser.extractObjectForTags( xml, "priorMu", double.class );
		priorAlpha = XMLParser.extractObjectForTags( xml, "priorAlpha", double.class );
		priorBeta = XMLParser.extractObjectForTags( xml, "priorBeta", double.class );
		fixMu = XMLParser.extractObjectForTags( xml, "fixMu", boolean.class );
		mu = XMLParser.extractObjectForTags( xml, "mu", double.class );
		logPrecision = XMLParser.extractObjectForTags( xml, "logPrecision", double.class );
		precision = Math.exp( logPrecision );
		initialized = XMLParser.extractObjectForTags( xml, "initialized", boolean.class );
		prop = XMLParser.extractObjectForTags( xml, "prop", DinucleotideProperty.class );
		precomputeNormalization();
	}
	
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags( buf, alphabets, "alphabet" );
		XMLParser.appendObjectWithTags( buf, ess, "ess" );
		XMLParser.appendObjectWithTags( buf, priorMu, "priorMu" );
		XMLParser.appendObjectWithTags( buf, priorAlpha, "priorAlpha" );
		XMLParser.appendObjectWithTags( buf, priorBeta, "priorBeta" );
		XMLParser.appendObjectWithTags( buf, fixMu, "fixMu" );
		XMLParser.appendObjectWithTags( buf, mu, "mu" );
		XMLParser.appendObjectWithTags( buf, logPrecision, "logPrecision" );
		XMLParser.appendObjectWithTags( buf, initialized, "initialized" );
		XMLParser.appendObjectWithTags( buf, prop, "prop" );
		XMLParser.addTags( buf, getClass().getSimpleName() );
		return buf;
	}

	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception {
		double val = (mu-priorMu);
		double gradmu = ess*precision*val;
		if( !fixMu ) {
			grad[start++] -= gradmu;
		}
		grad[start] += 0.5 - 0.5*gradmu*val + priorAlpha - priorBeta*precision;
		
	}

	public double getESS() {
		return ess;
	}

	public double getLogPriorTerm() {//TODO fixMu?
		double val = mu - priorMu;
		
		return 0.5*( Math.log( ess/(2.0*Math.PI) ) + logPrecision  - ess*precision*val*val )
					+ priorAlpha*Math.log( priorBeta ) - Gamma.logOfGamma( priorAlpha ) + priorAlpha*logPrecision - priorBeta*precision;
	}

	public double getLogNormalizationConstant() {
		return 0;
	}

	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception {
		return Double.NEGATIVE_INFINITY;
	}

	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index ) {
		return 0;
	}

	public double[] getCurrentParameterValues() throws Exception {
		if( fixMu ) {
			return new double[]{logPrecision};
		} else {
			return new double[]{mu,logPrecision};
		}
	}

	public String getInstanceName() {
		return getClass().getSimpleName()+" with "+mu+" and "+precision;
	}

	public double getLogScoreFor( Sequence seq, int start ) {
		if(Double.isNaN( seq.continuousVal( start )) ){
			return 0;
		}
		double val = 0;
		if(prop == null){
			val = seq.continuousVal( start ) - mu;
		}else{
			val = prop.getProperty( seq, start ) - mu;
		}
		return logNorm - 0.5*val*val*precision;
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer ) {
		if(Double.isNaN( seq.continuousVal( start )) ){
			return 0;
		}
		double val = 0;
		if(prop == null){
			val = seq.continuousVal( start ) - mu;
		}else{
			val = prop.getProperty( seq, start ) - mu;
		}
		
		int idx = 0;
		if( !fixMu ) {
			indices.add( idx++ );
			partialDer.add( precision*val );
		}
		indices.add( idx );
		partialDer.add( 0.5*(1.0 - precision*val*val) );
		
		return logNorm - 0.5*val*val*precision;
	}

	public int getNumberOfParameters() {
		return fixMu?1:2;
	}

	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		if(alwaysInitRandomly){
			initializeFunctionRandomly( freeParams );
		}else{
			double x = 0.0;
			double xsq = 0.0;
			double norm = 0.0;
			for(int i=0;i<data[index].getNumberOfElements();i++){
				double w = weights == null || weights[index] == null ? 1.0 : weights[index][i];
				//double temp = data[index].getElementAt( i ).continuousVal( 0 );
				if(!Double.isNaN( data[index].getElementAt( i ).continuousVal( 0 ) )){
					double temp = 0;
					if(prop == null){
						temp = data[index].getElementAt( i ).continuousVal( 0 ) - mu;
					}else{
						temp = prop.getProperty( data[index].getElementAt( i ), 0 ) - mu;
					}
					x += w*temp;
					xsq += w*temp*temp;
					norm += w;
				}
			}

			double var = xsq/norm - (x/norm)*(x/norm);
			if( !fixMu ) {
				mu = (x + ess*priorMu)/(norm + ess);
				if(Double.isNaN( mu ) && (norm == 0 || Double.isNaN( norm ))){
					mu = priorMu;
				}
			}
			
			double betap = priorBeta + 0.5*norm*var + (ess*norm*(x/norm - priorMu)*(x/norm - priorMu))/(2*(ess+norm));

			precision = (2*priorAlpha + norm - 1.0)/(2.0*betap);			
			if(Double.isNaN( precision ) && (norm == 0 || Double.isNaN( norm ))){
				precision = (2*priorAlpha-1.0)/(2.0*priorBeta);
			}
			
			logPrecision = Math.log( precision );
			
		//	System.out.println("plugin: prec: "+precision+", mu: "+mu);
			
			precomputeNormalization();
			initialized = true;
		}
	}

	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {
		precision = rand.nextGamma( priorAlpha, 1.0/priorBeta );
		logPrecision = Math.log( precision );
		if( !fixMu ) {
			mu = rand.nextGaussian()/precision + priorMu;
		}
		//System.out.println("prec: "+precision+", mu: "+mu);
		precomputeNormalization();
		initialized = true;
	}

	public boolean isInitialized() {
		return initialized;
	}

	public void setParameters( double[] params, int start ) {
		if( !fixMu ) {
			mu = params[start++];
		}
		logPrecision = params[start];
		precision = Math.exp( logPrecision );
		precomputeNormalization();
	}

	private void precomputeNormalization(){
		logNorm = 0.5*(logPrecision - Math.log( 2.0*Math.PI ));
	}
	
	public String toString(NumberFormat nf){
		return getClass().getSimpleName()+" with mu="+nf.format( mu )+" precision="+nf.format( precision )+(prop == null ? "" : " prop="+prop.name());
	}

}
