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
