package de.jstacs.classifiers.performanceMeasures;

import java.util.Arrays;

import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.utils.ToolBox;

/**
 * {@link PerformanceMeasure} using Pearson or Spearman correlation between prediction scores and
 * weighted class labels. The weighted class labels may be logit-transformed before computing correlaton.
 * 
 * @author Jan Grau
 *
 */
public class CorrelationCoefficient extends AbstractNumericalTwoClassPerformanceMeasure {

	/**
	 * The type of correlation used.
	 * @author Jan Grau
	 *
	 */
	public enum Method{
		/**
		 * Spearman correlation
		 */
		SPEARMAN,
		/**
		 * Pearson correlation
		 */
		PEARSON;
	}
	
	private Method method;
	private boolean logit;
	
	/**
	 * Creates a new {@link CorrelationCoefficient} using Spearman correlation and the raw
	 * weighted labels.
	 */
	public CorrelationCoefficient(){
		this(Method.SPEARMAN,false);
	}
	
	/**
	 * Creates a new {@link CorrelationCoefficient} using the suppled type of correlation
	 * and, optionally, logit transformation of weighted labels.
	 * @param method the type of correlation
	 * @param logit if weights should be logit transformed.
	 */
	public CorrelationCoefficient(Method method, boolean logit) {
		super();
		this.method = method;
		this.logit = logit;
	}

	
	
	
	@Override
	public String getName() {
		return "Correlation ("+method.name()+")";
	}
	
	


	@Override
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] weightsClass0, double[] sortedScoresClass1,
			double[] weightsClass1 ) {
		double[] temp = new double[sortedScoresClass0.length+sortedScoresClass1.length];
		if(weightsClass0 == null || weightsClass1 == null){
			Arrays.fill( temp, 0, sortedScoresClass0.length, 1 );
			Arrays.fill( temp, sortedScoresClass0.length, temp.length, 0 );
		}else{
			for(int i=0;i<weightsClass0.length;i++){
				if(logit){
					temp[i] = Math.log( weightsClass0[i]/(1.0-weightsClass0[i]) );
				}else{
					temp[i] = weightsClass0[i];
				}
			}
			for(int i=0;i<weightsClass1.length;i++){
				if(logit){
					temp[i+weightsClass0.length] = Math.log( (1.0-weightsClass1[i])/weightsClass1[i] );
				}else{
					temp[i+weightsClass0.length] = weightsClass1[i];
				}
			}
		}
		double[] temp2 = new double[sortedScoresClass0.length+sortedScoresClass1.length];
		System.arraycopy( sortedScoresClass0, 0, temp2, 0, sortedScoresClass0.length );
		System.arraycopy( sortedScoresClass1, 0, temp2, sortedScoresClass0.length, sortedScoresClass1.length );
		double cor = 0;
		try{
			if(method == Method.PEARSON){
				cor = ToolBox.pearsonCorrelation( temp, temp2 );
			}else{
				cor = ToolBox.spearmanCorrelation( sortedScoresClass0, weightsClass0 );
			}
			return new NumericalResultSet( new NumericalResult( getName(), "", cor ) );
		}catch(Exception e){
			throw new RuntimeException(e);
		}
	}

}
