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
import de.jstacs.utils.ToolBox;

/**
 * This class implements the false positive rate for a fixed sensitivity.
 * The false positive rate is defined as {@latex.inline $\\frac{FP}{TN+FP}$} and the sensitivity is defined as {@latex.inline $\\frac{TP}{TP+FN}$}. 
 * The classification threshold for computing false positive rate is chosen such that the classifier yields at least the specified sensitivity.
 * 
 * This measure corresponds to a specific point on the {@link ROCCurve}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class FalsePositiveRateForFixedSensitivity extends AbstractNumericalTwoClassPerformanceMeasure implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link FalsePositiveRateForFixedSensitivity} with empty parameter values.
	 */
	public FalsePositiveRateForFixedSensitivity() {
		super();
		try{
			parameters.add( new SimpleParameter( DataType.DOUBLE, "Sensitivity", "The fixed sensitivity for the false positive rate.", true, new NumberValidator<Double>(0d,1d),0.95 ) );
		}catch(ParameterException doesnothappen){}
	}
	
	/**
	 * Constructs a new instance of the performance measure {@link FalsePositiveRateForFixedSensitivity} with given <code>sensitivity</code>.
	 * 
	 * @param sensitivity the sensitivity for which the false positive rate should be computed
	 * 
	 * @throws Exception if the internal parameters can not be created or the value can not be set
	 */
	public FalsePositiveRateForFixedSensitivity(double sensitivity) throws Exception {
		this();
		getParameterAt( 0 ).setValue( sensitivity );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link FalsePositiveRateForFixedSensitivity} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link FalsePositiveRateForFixedSensitivity} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public FalsePositiveRateForFixedSensitivity( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "False positive rate for fixed sensitivity";
	}

	@Override
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] weightsClass0, double[] sortedScoresClass1, double[] weightsClass1 ) {
		double sensitivity = (Double)getParameterAt( 0 ).getValue();
		double threshold = findThreshold( sortedScoresClass0, sortedScoresClass1, weightsClass0, 1-sensitivity, false );
		double fpr;
		int i = findSplitIndex( sortedScoresClass1, threshold );
		int d = sortedScoresClass1.length;
		if( weightsClass1 == null ) {
			fpr = (double)( d - i ) / (double)d;
		} else {
			fpr = ToolBox.sum( i, d, weightsClass1 );
			fpr = fpr / (ToolBox.sum( 0, i, weightsClass1 ) + fpr); 
		}

		return new NumericalResultSet( new NumericalResult[]{
				new NumericalResult( getName() + " of "+sensitivity, "", fpr ),
				new NumericalResult( "Threshold for the " + getName().toLowerCase() + " of "+sensitivity, "", threshold )
		});
	}
}
