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

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;

/**
 * This class implements the positive predictive value for a fixed sensitivity.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class PositivePredictiveValueForFixedSensitivity extends TwoClassAbstractPerformanceMeasure implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link PositivePredictiveValueForFixedSensitivity} with empty parameter values.
	 */
	public PositivePredictiveValueForFixedSensitivity() {
		super();
		try{
			parameters.add( new SimpleParameter( DataType.DOUBLE, "Sensitivity", "The fixed sensitivity for the positive predictive value.", true, new NumberValidator<Double>(0d,1d),0.95 ) );
		}catch(ParameterException doesnothappen){}
	}
	
	/**
	 * Constructs a new instance of the performance measure {@link PositivePredictiveValueForFixedSensitivity} with given <code>sensitivity</code>.
	 * 
	 * @param sensitivity the sensitivity for which the positive predictive value should be computed
	 * 
	 * @throws Exception if the internal parameters can not be created or the value can not be set
	 */
	public PositivePredictiveValueForFixedSensitivity(double sensitivity) throws Exception {
		this();
		getParameterAt( 0 ).setValue( sensitivity );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link PositivePredictiveValueForFixedSensitivity} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link PositivePredictiveValueForFixedSensitivity} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public PositivePredictiveValueForFixedSensitivity( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Positive predictive value for fixed sensitivity";
	}

	@Override
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
		double sensitivity = (Double)getParameterAt( 0 ).getValue();
		int d = sortedScoresClass1.length, j = (int)Math.ceil( ( 1 - sensitivity ) * ( sortedScoresClass0.length - 1 ) ), i = d - 1;
		double threshold = sortedScoresClass0[j];
		while( j >= 0 && sortedScoresClass0[j] == threshold ) {
			j--;
		}
		// => (j+1) false negatives
		// => (scoresClass0.length-1-j) true positives
		j = sortedScoresClass0.length - 1 - j; // true positives
		while( i >= 0 && sortedScoresClass1[i] >= threshold ) {
			i--;
		}
		// => (i+1) true negatives
		// => (d-1-i) false positives
		i = d - 1 - i; // false positives
		return new NumericalResultSet(new NumericalResult[]{
				new NumericalResult( getName() + " of "+sensitivity, "", ( (double)j ) / (double)( i + j ) ),
				new NumericalResult( "Threshold for the " + getName().toLowerCase() + " of "+sensitivity, "", threshold ),
		});
	}
	
	public NumericalResultSet compute( double[][][] classSpecificScores ) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}

}
