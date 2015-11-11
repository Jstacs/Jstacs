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
package de.jstacs.algorithms.alignment.cost;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * Class for matrix costs, i.e., the cost of any match/mismatch is stored in
 * a matrix allowing a huge degree of freedom.
 * 
 * @author Jens Keilwagen
 */
public class MatrixCosts implements Costs {

	private double[][] matrix;
	private double gap;
	
	/**
	 * Creates a new instance of {@link MatrixCosts} where the costs
	 * for mismatch and match are given in <code>matrix</code>.
	 * Additionally, the costs for a gap can be specified.
	 *
	 * @param matrix
	 *            the match and mismatch costs
	 * @param gap
	 *            the cost for a gap 
	 *            
	 * @throws CloneNotSupportedException if <code>matrix</code> could not be cloned
	 */
	public MatrixCosts( double[][] matrix, double gap ) throws CloneNotSupportedException {
		this.matrix = ArrayHandler.clone( matrix );
		this.gap = gap;
	}
	
	/**
	 * Restores {@link MatrixCosts} object from its XML representation.
	 * @param xml the XML representation
	 * @throws NonParsableException if the XML could not be parsed
	 */
	public MatrixCosts( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, "MatrixCosts" );
		matrix = (double[][])XMLParser.extractObjectForTags( xml, "matrix" );
		gap = (Double)XMLParser.extractObjectForTags( xml, "gap" );
	}

	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, matrix, "matrix" );
		XMLParser.appendObjectWithTags( xml, gap, "gap" );
		XMLParser.addTags( xml, "MatrixCosts" );
		return xml;
	}
	
	public double getCostFor( Sequence s1, Sequence s2, int i, int j ) {
		return matrix[s1.discreteVal( i - 1 )][s2.discreteVal( j - 1 )];
	}

	@Override
	public double getGapCosts() {
		return gap;
	}
}