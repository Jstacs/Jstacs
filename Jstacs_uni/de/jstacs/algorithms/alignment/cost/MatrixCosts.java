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

import de.jstacs.data.Sequence;
import de.jstacs.io.ArrayHandler;

/**
 * Class for matrix costs, i.e., the cost of any match/mismatch is stored in
 * a matrix allowing a huge degree of freedom. Additionally, the class used
 * affine gap costs. The costs <code>start</code> to start a new gap, and
 * the costs <code>elong</code> to elongate a gap by one position.
 * 
 * @author Jens Keilwagen
 */
public class MatrixCosts extends AffineCosts {

	private double[][] matrix;
	
	/**
	 * Creates a new instance of matrix costs where the costs
	 * for mismatch and match are given in <code>matrix</code>.
	 * Additionally, costs <code>start</code> to start a new gap,
	 * and costs <code>elong</code> to elongate a gap by one
	 * position.
	 *
	 * @param matrix
	 *            the match and mismatch costs
	 * @param start
	 *            the costs to start a gap
	 * @param elong
	 *            the costs to elongate a gap
	 *            
	 * @throws CloneNotSupportedException if <code>matrix</code> could not be cloned
	 */
	public MatrixCosts( double[][] matrix, double start, double elong ) throws CloneNotSupportedException {
		super( start, elong );
		this.matrix = ArrayHandler.clone( matrix );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.Alignment.Costs#getCostFor(de.jstacs.data.Sequence, de.jstacs.data.Sequence, int, int, de.jstacs.algorithms.Alignment.Costs.Direction)
	 */
	public double getCostFor( Sequence s1, Sequence s2, int i, int j, Direction from ) {
		if( from == Direction.TOP || from == Direction.LEFT ) {
			return getGapCostsFor(1);
		} else {
			return matrix[s1.discreteVal( i - 1 )][s2.discreteVal( j - 1 )];
		}
	}
}