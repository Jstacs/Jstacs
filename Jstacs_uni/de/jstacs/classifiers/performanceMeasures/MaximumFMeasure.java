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

//http://en.wikipedia.org/wiki/F1_score

/**
 * Computes the maximum of the general F-measure given a positive real parameter {@latex.inline $\\beta$}.
 * The F-measure is defined as {@latex.inline $F_\\beta = (1+\\beta^2) \\frac{PPV*Sn}{(\\beta^2 PPV) + Sn}$}, where PPV denotes the
 * positive predictive value and Sn denotes sensitivity. For {@latex.inline $\\beta = 1$} it is equal to the traditional F1-measure.
 * 
 * @author Jens Keilwagen
 * @see PositivePredictiveValueForFixedSensitivity
 * @see SensitivityForFixedSpecificity
 */
public class MaximumFMeasure extends MaximumNumericalTwoClassMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link MaximumFMeasure} with empty parameters.
	 * 
	 */
	public MaximumFMeasure() {
		super();
		try{
			parameters.add( new SimpleParameter( DataType.DOUBLE, "beta", "the beta defining the F measure", true, new NumberValidator<Double>(0d,Double.POSITIVE_INFINITY), 1d ) );
		}catch(ParameterException doesnothappen){ }	
	}
	
	/**
	 * Constructs a new instance of the performance measure {@link MaximumFMeasure} with given <code>beta</code>.
	 * 
	 * @param beta the beta for which the maximum F-measure should be computed
	 * 
	 * @throws Exception if the internal parameters can not be created or the value can not be set
	 */
	public MaximumFMeasure( double beta ) throws Exception {
		this();
		parameters.get(0).setValue( beta );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link MaximumFMeasure} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MaximumFMeasure} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public MaximumFMeasure( StringBuffer xml ) throws NonParsableException {
		super(xml);
	}

	@Override
	protected String getMeasureName() {
		return "F-Measure";
	}
	
	@Override
	protected String getSpecificName() {
		return getMeasureName() + " with beta=" + parameters.get("beta").getValue();
	}

	@Override
	protected double getMeasure(double tp, double fp, double fn, double tn) {
		double b = (Double) parameters.get( "beta").getValue();
		b *= b;
		double precision = tp / (tp+fp);
		double recall = tp / (tp+fn);
		return (1+b) * precision*recall / (b*precision + recall);
	}
}