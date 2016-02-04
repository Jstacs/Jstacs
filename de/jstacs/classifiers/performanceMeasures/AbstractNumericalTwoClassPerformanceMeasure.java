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
import de.jstacs.results.NumericalResultSet;

/**
 *  This class is the abstract super class of any performance measure that can only be computed for binary classifiers.
 * 
 * @author Jens Keilwagen
 */
public abstract class AbstractNumericalTwoClassPerformanceMeasure extends AbstractTwoClassPerformanceMeasure implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new {@link AbstractNumericalTwoClassPerformanceMeasure} with empty parameter values.
	 */
	protected AbstractNumericalTwoClassPerformanceMeasure() {
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link AbstractNumericalTwoClassPerformanceMeasure} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AbstractNumericalTwoClassPerformanceMeasure} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	protected AbstractNumericalTwoClassPerformanceMeasure(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	public final NumericalResultSet compute( double[][][] classSpecificScores ) {
		return compute( classSpecificScores, null );
	}
	
	public final NumericalResultSet compute( double[][][] classSpecificScores, double[][] weights ) {
		return (NumericalResultSet) super.compute( classSpecificScores, weights );
	}

	@Override
	public final NumericalResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
		return (NumericalResultSet) compute( sortedScoresClass0, null, sortedScoresClass1, null );
	}
	
	public abstract NumericalResultSet compute(double[] sortedScoresClass0, double[] weightClass0, double[] sortedScoresClass1, double[] weightsClass1 );
}
