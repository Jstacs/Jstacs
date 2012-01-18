package de.jstacs.classifiers.performanceMeasures;

import de.jstacs.results.NumericalResultSet;

/**
 * This interface indicates that a Performance measure returns numerical results.
 * 
 * <b>Any class that implements this interface should be an extension of {@link AbstractPerformanceMeasure}.</b>
 * 
 * @see NumericalResultSet
 * @see de.jstacs.results.NumericalResult
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public interface NumericalPerformanceMeasure {

	/**
	 * This method allows to compute the performance measure of given sorted score ratios.
	 * 
	 * <b>This method can only be used for binary classifiers.</b>
	 * 
	 * @param sortedScoresClass0 the sorted score ratios of class 0
	 * @param sortedScoresClass1 the sorted score ratios of class 1
	 *  
	 * @return a numerical result set containing the results of the performance measure
	 * 
	 * @see java.util.Arrays#sort(double[])
	 */
	public NumericalResultSet compute(double[] sortedScoresClass0, double[] sortedScoresClass1);

	/**
	 * This method allows to compute the performance measure of given class specific scores.
	 * 
	 * @param classSpecificScores the scores; first dimension = data sets, second dimension = sequences of the data set, third dimension classes of the classifier
	 *  
	 * @return a numerical result set containing the results of the performance measure
	 */
	public NumericalResultSet compute(double[][][] classSpecificScores);
}
