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

import de.jstacs.Storable;
import de.jstacs.data.sequences.Sequence;

/**
 * The general interface for the costs of an alignment.
 * 
 * @author Jan Grau
 */
public interface Costs extends Storable {
	/**
	 * Returns the costs for the alignment of <code>s1(i)</code> and
	 * <code>s2(j)</code>.
	 * 
	 * @param s1
	 *            the first sequence
	 * @param s2
	 *            the second sequence
	 * @param i
	 *            the index in the first sequence
	 * @param j
	 *            the index in the second sequence
	 * 
	 * @return the costs
	 * 
	 * @see Sequence#discreteVal(int)
	 */
	public double getCostFor( Sequence s1, Sequence s2, int i, int j );
	
	/**
	 * Returns the costs for an insert gap, i.e., a gap in the first string.
	 *  
	 * @return the costs for an insert gap
	 */
	public double getInsertCosts();
	
	/**
	 * Returns the costs for a delete gap, i.e., a gap in the second string.
	 *  
	 * @return the costs for a delete gap
	 */
	public double getDeleteCosts();
	
	
}