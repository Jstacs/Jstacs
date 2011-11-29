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

/**
 * The general interface for the costs of an alignment.
 * 
 * @author Jan Grau
 */
public interface Costs {

	/**
	 * The direction of the predecessor in the DP (<b>d</b>ynamic
	 * <b>p</b>rogramming)-matrix.
	 * 
	 * @author Jan Grau
	 */
	public enum Direction {

		/**
		 * The predecessor is above the current element.
		 */
		TOP,
		/**
		 * The predecessor is left of the current element.
		 */
		LEFT,
		/**
		 * The predecessor is top left of the current element.
		 */
		DIAGONAL,
		/**
		 * The predecessor is at the same position as the current element.
		 */
		SELF;

	}

	/**
	 * Returns the costs for the alignment if <code>s1(i)</code> and
	 * <code>s2(j)</code> coming from <code>from</code>.
	 * 
	 * @param s1
	 *            the first sequence
	 * @param s2
	 *            the second sequence
	 * @param i
	 *            the index in the first sequence
	 * @param j
	 *            the index in the second sequence
	 * @param from
	 *            the direction within the DP-matrix
	 * 
	 * @return the costs
	 */
	public double getCostFor( Sequence s1, Sequence s2, int i, int j, Direction from );

	/**
	 * Returns the costs for a gap of length <code>length</code>.
	 * 
	 * @param length
	 *            the length of the gap
	 * 
	 * @return the corresponding costs
	 */
	public double getGapCostsFor( int length );

	/**
	 * Returns the costs to elongate a gap by one position.
	 * 
	 * @return the corresponding costs
	 */
	public double getElongateCosts();
}