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


/**
 * This class implements affine gap costs, i.e., the costs for starting a new gap are given by <code>start</code>, and
 * the costs for elongating a gap by one position are given by <code>elong</code>.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public abstract class AffineCosts implements Costs {
	private double start;
	private double elong;
	
	/**
	 * This constructor must be used by any constructor of any subclass to define the affine gap costs of an alignment.
	 * 
	 * @param start the costs for starting a gap
	 * @param elong the cost for elongating a gap
	 */
	protected AffineCosts( double start, double elong ) {
		if( start < 0 && -start > elong ) {
			throw new IllegalArgumentException( "Problem: start < 0 && -start > elong" );
		}
		this.start = start;
		this.elong = elong;
	}
	

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.Alignment.Costs#getElongateCosts()
	 */
	public double getElongateCosts() {
		return elong;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.Alignment.Costs#getGapCostsFor(int)
	 */
	public double getGapCostsFor( int length ) {
		return start + ( length * elong );
	}
}
