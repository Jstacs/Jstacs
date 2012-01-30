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
 * This class implements affine gap costs, i.e., the costs for starting a new gap are given by <code>start</code>, and
 * the costs for elongating a gap by one position are given by <code>elong</code>.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class AffineCosts implements Costs {
	private double start;
	private double elong;
	private Costs c;
	
	/**
	 * This constructor creates a new instance of cost using affine gap costs.
	 * The costs for match, mismatch and gap elongation are defined by <code>c</code>,
	 * while the costs for a start of a gap are given by <code>start</code>.
	 * 
	 * @param start the costs for starting a gap
	 * @param c the cost for match, mismatch and gap elongation
	 */
	public AffineCosts( double start, Costs c ) {
		this.c = c;
		this.elong = c.getGapCosts();
		if( start < 0 && -start > elong ) {
			throw new IllegalArgumentException( "Problem: start < 0 && -start > elong" );
		}
		this.start = start;
	}
	
	/**
	 * Returns the costs to elongate a gap by one position.
	 * 
	 * @return the corresponding costs
	 */
	public double getElongateCosts() {
		return elong;
	}

	/**
	 * Returns the costs for a gap of length <code>length</code>.
	 * 
	 * @param length
	 *            the length of the gap
	 * 
	 * @return the corresponding costs
	 */
	public double getGapCostsFor( int length ) {
		return start + ( length * elong );
	}
	
	@Override
	public double getCostFor(Sequence s1, Sequence s2, int i, int j) {
		return c.getCostFor(s1, s2, i, j);
	}

	@Override
	public double getGapCosts() {
		return c.getGapCosts();
	}
}