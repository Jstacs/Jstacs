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
 * This class implements the area under curve of the Receiver Operating Characteristics curve.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see ROCCurve
 */
public class AucROC extends ROCCurve implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link AucROC}.
	 */
	public AucROC() {
		super();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link AucROC} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AucROC} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public AucROC(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
		return (NumericalResultSet) super.compute( sortedScoresClass0, sortedScoresClass1 );
	}

	public NumericalResultSet compute( double[][][] classSpecificScores ) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}
	
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] weightsClass0, double[] sortedScoresClass1, double[] weightsClass1 ) {
		return (NumericalResultSet) super.compute( sortedScoresClass0, weightsClass0, sortedScoresClass1, weightsClass1 );
	}

	public NumericalResultSet compute( double[][][] classSpecificScores, double[][] weights ) {
		return (NumericalResultSet) super.compute( classSpecificScores, weights );
	}
	
	@Override
	public String getName() {
		return "Area under curve (" + NAME + ")";
	}
}