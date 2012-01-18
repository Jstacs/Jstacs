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

/**
 * Class for simple costs with costs <code>mismatch</code> for a mismatch,
 * costs <code>start</code> to start a new gap, costs <code>elong</code> to
 * elongate a gap by one position and costs of <code>match</code> for a match.
 */
public class SimpleCosts extends AffineCosts {

	private double match;
	
	private double mismatch;

	/**
	 * Creates a new instance of simple costs with costs
	 * <code>mismatch</code> for a mismatch, costs <code>start</code> to
	 * start a new gap, costs <code>elong</code> to elongate a gap by one
	 * position and costs of <code>match</code> for a match.
	 *
	 * @param match
	 *            the match costs
	 * @param mismatch
	 *            the mismatch costs
	 * @param start
	 *            the costs to start a gap
	 * @param elong
	 *            the costs to elongate a gap
	 */
	public SimpleCosts( double match, double mismatch, double start, double elong ) {
		super( start, elong );
		if( mismatch <= 0 ) {
			throw new IllegalArgumentException( "Problem: mismatch <= 0" );
		}
		this.mismatch = mismatch;
		if( match > 0 ) {
			throw new IllegalArgumentException( "Problem: match > 0" );
		}
		this.match = match;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.Alignment.Costs#getCostFor(de.jstacs.data.Sequence, de.jstacs.data.Sequence, int, int, de.jstacs.algorithms.Alignment.Costs.Direction)
	 */
	public double getCostFor( Sequence s1, Sequence s2, int i, int j, Direction from ) {
		if( from == Direction.TOP || from == Direction.LEFT ) {
			return getGapCostsFor(1);
		} else {
			if( s1.discreteVal( i - 1 ) != s2.discreteVal( j - 1 ) ) {
				return mismatch;
			} else {
				return match;
			}
		}
	}
}
