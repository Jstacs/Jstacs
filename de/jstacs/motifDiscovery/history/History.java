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

package de.jstacs.motifDiscovery.history;

import de.jstacs.Storable;

/**
 * This interface is used to manage the history of some process.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public interface History extends Cloneable, Storable {
	
	/**
	 * Returns <code>true</code> if the specified operation is allowed, i.e. it does not conflict with any operation from the history.
	 * 
	 * @param op the operation to be tested
	 * 
	 * @return <code>true</code> if the specified operation is allowed
	 */
	public boolean operationAllowed(int... op);
	
	/**
	 * This method puts an operation to the history.
	 *  
	 * @param op the performed operation
	 */
	public void operationPerfomed(int... op);
	
	/**
	 * This method clears the history, i.e. it removes all operations from the history.
	 */
	public void clear();
	
	/**
	 * This method returns a deep copy of the instance
	 * 
	 * @return a deep copy
	 * 
	 * @throws CloneNotSupportedException if not possible
	 * 
	 * @see Cloneable
	 * @see Object#clone()
	 */
	public History clone() throws CloneNotSupportedException;
}