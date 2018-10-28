package de.jstacs.sequenceScores.statisticalModels.differentiable.continuous;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Random;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.continuous.gamma.GammaPriorFunction;
import de.jstacs.sequenceScores.statisticalModels.differentiable.continuous.gamma.NumericalIntegration;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;


public class ExpGammaDiffSM extends AbstractDifferentiableStatisticalModel {

	private static Random r = new Random();
	private double[] alphas;
	private double[] betas;
	private boolean isInitialized;
	private boolean plugin;
	private double ess;
	private double[] mua;
	private double[] mug;
	
	private double norm;
	private double priorNorm;
	private double[] alphaNorms;
	
	public ExpGammaDiffSM(int length, double ess, double[] mua, double[] mug, boolean plugin){
		super(new AlphabetContainer(new ContinuousAlphabet()),length);
		alphas = new double[length];
		betas = new double[length];
		this.alphaNorms = new double[length];
		isInitialized = false;
		this.plugin = plugin;
		this.ess = ess;
		this.mua = mua.clone();
		this.mug = mug.clone();
		//precomputePriorNormalization();
		priorNorm = 1.0;
	}
	
	public ExpGammaDiffSM(StringBuffer buf) throws NonParsableException{
		super(buf);
	}
	
	
	public ExpGammaDiffSM clone() throws CloneNotSupportedException{
		ExpGammaDiffSM clone = (ExpGammaDiffSM) super.clone();
		clone.alphas = alphas.clone();
		clone.betas = betas.clone();
		clone.mua = mua.clone();
		clone.mug = mug.clone();
		clone.alphaNorms = alphaNorms.clone();
		return clone;
	}
	

	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags( buf, alphabets, "alphabets" );
		XMLParser.appendObjectWithTags( buf, length, "length" );
		XMLParser.appendObjectWithTags( buf, alphas, "alphas" );
		XMLParser.appendObjectWithTags( buf, betas, "betas" );
		XMLParser.appendObjectWithTags( buf, isInitialized, "isInitialized" );
		XMLParser.appendObjectWithTags( buf, plugin, "plugin" );
		XMLParser.appendObjectWithTags( buf, ess, "ess" );
		XMLParser.appendObjectWithTags( buf, mua, "mua" );
		XMLParser.appendObjectWithTags( buf, mug, "mug" );
		XMLParser.appendObjectWithTags( buf, norm, "norm" );
		XMLParser.appendObjectWithTags( buf, priorNorm, "priorNorm" );
		XMLParser.appendObjectWithTags( buf, alphaNorms, "alphaNorms" );
		XMLParser.addTags( buf, getClass().getSimpleName() );
		return buf;
	}
	
	@Override
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		alphabets = XMLParser.extractObjectForTags( xml, "alphabets", AlphabetContainer.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		alphas = XMLParser.extractObjectForTags( xml, "alphas", double[].class );
		betas = XMLParser.extractObjectForTags( xml, "betas", double[].class );
		isInitialized = XMLParser.extractObjectForTags( xml, "isInitialized", boolean.class );
		plugin = XMLParser.extractObjectForTags( xml, "plugin", boolean.class );
		ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
		mua = XMLParser.extractObjectForTags( xml, "mua", double[].class );
		mug = XMLParser.extractObjectForTags( xml, "mug", double[].class );
		norm = XMLParser.extractObjectForTags( xml, "norm", double.class );
		priorNorm = XMLParser.extractObjectForTags( xml, "priorNorm", double.class );
		alphaNorms = XMLParser.extractObjectForTags( xml, "alphaNorms", double[].class );
	}

	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception {
		for(int i=0;i<alphas.length;i++){
			//normalization for class parameters
			//grad[start+i] += getPartialNormalizationConstant( i );
			//grad[start+alphas.length+i] += getPartialNormalizationConstant( i+alphas.length );
			
			//alphas
			try{
			grad[start+i] += ess*( -alphas[i]*digamma( alphas[i] ) + Math.log( betas[i] )*alphas[i] + alphas[i]*Math.log( mug[i] ) ) + 1;
			}catch(StackOverflowError e){
				System.out.println(alphas[i]);
				throw e;
			}
			//betas
			grad[start+alphas.length+i] += ess*( alphas[i] - mua[i]*betas[i] ) + 1;
			
		}
	}

	public double getESS() {
		return ess;
	}

	public double getLogPriorTerm() {
		double val = 0;
		for(int i=0;i<alphas.length;i++){
			double temp = -Gamma.logOfGamma( alphas[i] )
			+ Math.log( betas[i] )*alphas[i] 
			                              - mua[i]*betas[i] 
			                                             + alphas[i]*Math.log( mug[i] );
			temp *= ess;
			temp += Math.log( alphas[i]) + Math.log( betas[i] );
			val += temp;
		}
		
		return val-priorNorm;
	}

	
	public double getLogNormalizationConstant(){
		return 0;
	}
	
	public double getLogPartialNormalizationConstant(int parameterIndex){
		return Double.NEGATIVE_INFINITY;
	}

	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index ) {
		return 1;
	}

	public double[] getCurrentParameterValues() throws Exception {
		double[] pars = new double[alphas.length+betas.length];
		for(int i=0;i<alphas.length;i++){
			pars[i] = Math.log( alphas[i] );
			pars[i + alphas.length] = Math.log( betas[i] );
		}
		return pars;
	}

	public String getInstanceName() {
		return getClass().getSimpleName();
	}
	
	private void precomputePriorNormalization() {
		priorNorm = 0;
		if(ess > 0){
			try{
				for(int i=0;i<alphas.length;i++){
					double temp = NumericalIntegration.getIntegralByNestedIntervals( new GammaPriorFunction(mua[i], mug[i], ess ), 1E-10, 1E-2, 0.001 );
					priorNorm += Math.log( temp );
				}
			}catch(Exception e){
				e.printStackTrace();
			}
		}
	}
	
	private void precomputeNormalization(){
		norm = 0;
		for(int i=0;i<alphas.length;i++){
			norm += alphas[i]*Math.log( betas[i] ) - Gamma.logOfGamma( alphas[i] );
			alphaNorms[i] = alphas[i]*Math.log( betas[i] ) - alphas[i]*digamma( alphas[i] );
		}
	}

	public double getLogScoreFor( Sequence seq, int start ) {
		double val = norm;
		for(int i=0;i<alphas.length;i++){
			double cv = seq.continuousVal( i+start )+1E-10;
			val += (alphas[i]-1.0)*Math.log( cv ) - betas[i]*cv;
			if(Double.isInfinite(val)|| Double.isNaN(val)){
				System.out.println(val+" "+Arrays.toString(alphas)+" "+Arrays.toString(betas)+" "+cv);
			}
		}
		
		return val;
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer ) {
		double val = norm;
		for(int i=0;i<alphas.length;i++){
			double cv = seq.continuousVal( i+start )+1E-10;
			val += (alphas[i]-1.0)*Math.log( cv ) - betas[i]*cv;
			indices.add( i );
			partialDer.add( alphaNorms[i] + 
					alphas[i]*Math.log( cv ) );
			indices.add( i+alphas.length );
			partialDer.add( alphas[i]
					-betas[i]*cv );
		}
		return val;
	}

	public int getNumberOfParameters() {
		return alphas.length + betas.length;
	}

	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {

		if(plugin){
			Arrays.fill( alphas, 0 );
			Arrays.fill( betas, 0 );
			double norm = 0;
			for(int i=0;i<data[index].getNumberOfElements();i++){
				Sequence seq = data[index].getElementAt( i );
				double w = weights == null || weights[index] == null ? 1.0 : weights[index][i];
				for(int j=0;j<seq.getLength();j++){
					alphas[j] += w*seq.continuousVal( j );
					betas[j] += w*seq.continuousVal( j )*seq.continuousVal( j );
				}
				norm += w;
			}
			for(int i=0;i<alphas.length;i++){
				alphas[i] /= norm;
				betas[i] /= norm;
				double theta = (betas[i]/alphas[i] - alphas[i]);
				double k = ((betas[i]-alphas[i]*alphas[i]) - alphas[i])/(theta*(theta-1.0));
				alphas[i] = k;
				if(Double.isInfinite( alphas[i] )){
					alphas[i] = 1.0;
				}
				betas[i] = 1.0/theta;
				if(Double.isInfinite( betas[i] )){
					betas[i] = 1.0;
				}
			}

			precomputeNormalization();
			isInitialized = true;
		}else{
			initializeFunctionRandomly( freeParams );
			//Arrays.fill( alphas, 1 );
			//Arrays.fill( betas, 1 );
		}
	}

	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {

		for(int i=0;i<alphas.length;i++){
			alphas[i] = r.nextDouble()*10 + Double.MIN_VALUE;
			betas[i] = r.nextDouble();
		}
		precomputeNormalization();
		isInitialized = true;
	}

	public boolean isInitialized() {
		return isInitialized;
	}

	public void setParameters( double[] params, int start ) {
		for(int i=0;i<alphas.length;i++){
			alphas[i] = Math.exp( params[start+i] );
			betas[i] = Math.exp( params[start+alphas.length+i] );
		}
		precomputeNormalization();
	}
	
	public String toString(NumberFormat nf){
		String str = getClass().getSimpleName()+"\nalphas: "+Arrays.toString( alphas )+"\n";
		str = str + "betas: "+Arrays.toString( betas )+"\n";
		return str;
	}
	
	private double digamma (double x) {
		  if(Double.isNaN( x ) || Double.isInfinite( x )){
			  System.out.flush();
			  System.out.println(Arrays.toString( alphas ));
			  System.out.println(Arrays.toString( betas ));
			  System.out.println(Arrays.toString( mua ));
			  System.out.println(Arrays.toString( mug ));
			  System.out.println(norm+", "+priorNorm+","+Arrays.toString( alphaNorms ));
			  System.out.flush();
			  throw new RuntimeException("argument of digamma is "+x);
		  }
	      final double C7[][] = {
	       {1.3524999667726346383e4, 4.5285601699547289655e4, 4.5135168469736662555e4,
	        1.8529011818582610168e4, 3.3291525149406935532e3, 2.4068032474357201831e2,
	        5.1577892000139084710, 6.2283506918984745826e-3},
	       {6.9389111753763444376e-7, 1.9768574263046736421e4, 4.1255160835353832333e4,
	          2.9390287119932681918e4, 9.0819666074855170271e3,
	          1.2447477785670856039e3, 6.7429129516378593773e1, 1.0}
	      };
	      final double C4[][] = {
	       {-2.728175751315296783e-15, -6.481571237661965099e-1, -4.486165439180193579,
	        -7.016772277667586642, -2.129404451310105168},
	       {7.777885485229616042, 5.461177381032150702e1,
	        8.929207004818613702e1, 3.227034937911433614e1, 1.0}
	      };

	      double prodPj = 0.0;
	      double prodQj = 0.0;
	      double digX = 0.0;

	      if (x >= 3.0) {
	         double x2 = 1.0 / (x * x);
	         for (int j = 4; j >= 0; j--) {
	            prodPj = prodPj * x2 + C4[0][j];
	            prodQj = prodQj * x2 + C4[1][j];
	         }
	         digX = Math.log (x) - (0.5 / x) + (prodPj / prodQj);

	      } else if (x >= 0.5) {
	         final double X0 = 1.46163214496836234126;
	         for (int j = 7; j >= 0; j--) {
	            prodPj = x * prodPj + C7[0][j];
	            prodQj = x * prodQj + C7[1][j];
	         }
	         digX = (x - X0) * (prodPj / prodQj);

	      } else {
	         double f = (1.0 - x) - Math.floor (1.0 - x);
	         digX = digamma (1.0 - x) + Math.PI / Math.tan (Math.PI * f);

	      }

	      return digX;
	   }

}
