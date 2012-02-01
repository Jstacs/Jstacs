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
 * This class implements the sensitivity for a fixed specificity.
 * The sensitivity is defined as {@latex.inline $\\frac{TP}{TP+FN}$} and the specificity is defined as {@latex.inline $\\frac{TN}{TN+FP}$}. 
 * The classification threshold for computing the sensitivity is chosen such that the classifier yields at least the specified specificity.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class SensitivityForFixedSpecificity extends TwoClassAbstractPerformanceMeasure implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link SensitivityForFixedSpecificity} with empty parameter values.
	 */
	public SensitivityForFixedSpecificity() {
		super();
		try {
			parameters.add( new SimpleParameter( DataType.DOUBLE, "Specificity", "The fixed specificity for the sensitivity.", true, new NumberValidator<Double>(0d,1d),0.999 ) );
		} catch ( ParameterException doesnothappen ) { }
	}
	
	/**
	 * Constructs a new instance of the performance measure {@link SensitivityForFixedSpecificity} with given <code>specificity</code>.
	 * 
	 * @param specificity the specificity for which the sensitivity should be computed
	 * 
	 * @throws Exception if the internal parameters can not be created or the value can not be set
	 */
	public SensitivityForFixedSpecificity(double specificity) throws Exception {
		this();
		getParameterAt( 0 ).setValue( specificity );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link SensitivityForFixedSpecificity} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SensitivityForFixedSpecificity} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public SensitivityForFixedSpecificity( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Sensitivity for a fixed specificity";
	}

	@Override
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
		double specificity = (Double)getParameterAt( 0 ).getValue();
		int i = 0, m = sortedScoresClass0.length;
		double threshold = sortedScoresClass1[(int)Math.ceil( specificity * ( sortedScoresClass1.length - 1 ) )];
		while( i < m && sortedScoresClass0[i] <= threshold ) {
			i++;
		}
		
		return new NumericalResultSet(new NumericalResult[]{
				new NumericalResult( getName() +" of "+specificity, "", (double)( m - i ) / (double)m  ),
				new NumericalResult( "Threshold for the "+getName().toLowerCase() +" of "+specificity, "", threshold )
		});
	}

	public NumericalResultSet compute( double[][][] classSpecificScores ) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}
}
