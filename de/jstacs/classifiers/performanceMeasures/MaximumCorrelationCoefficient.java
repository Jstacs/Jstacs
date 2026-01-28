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

import de.jstacs.io.NonParsableException;

/**
 * This class implements the maximum of the correlation coefficient \( \frac{ TP\cdot TN - FN \cdot FP }{ \sqrt{ (TP+FN)\cdot(TN+FP)\cdot(TP+FP)\cdot(TN+FN) } } \).
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class MaximumCorrelationCoefficient extends MaximumNumericalTwoClassMeasure implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link MaximumCorrelationCoefficient}.
	 */
	public MaximumCorrelationCoefficient() {}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link MaximumCorrelationCoefficient} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MaximumCorrelationCoefficient} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public MaximumCorrelationCoefficient( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	protected String getMeasureName() {
		return "correlation coefficient";
	}
	
	@Override
	protected String getSpecificName() {
		return getMeasureName();
	}

	// double since otherwise an overflow is easily possible
	protected double getMeasure( double tp, double fp, double fn, double tn ) {
		return ( tp * tn - fn * fp ) / Math.sqrt( ( tp + fn ) * ( tn + fp ) * ( tp + fp ) * ( tn + fn ) );
	}
}
