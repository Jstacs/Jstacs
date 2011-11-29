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

package de.jstacs.sampling;

import de.jstacs.Storable;

/**
 * This is the abstract super class for any test of the length of the burn-in
 * phase.
 * 
 * @author Jens Keilwagen
 */
public interface BurnInTest extends Cloneable, Storable {

	/**
     * Return a deep copy of this object.
     * 
     * @return a deep copy of this object
     * 
	 * @throws CloneNotSupportedException if the object can not be cloned 
     * 
     * @see Object#clone()
     */
	public BurnInTest clone() throws CloneNotSupportedException;

	/**
	 * This method sets the value of the current sampling. This allows to assign
	 * the values from {@link BurnInTest#setValue(double)} to a sampling.
	 * 
	 * @param index
	 *            the index of the sampling
	 */
	public void setCurrentSamplingIndex( int index );

	/**
	 * This method can be used to fill the internal memory with the values that
	 * will be used to determine the length of the burn-in phase.
	 * 
	 * @param val
	 *            the value
	 */
	public void setValue( double val );

	/**
	 * This method can be used to remove all values from the internal memory.
	 */
	public void resetAllValues();

	/**
	 * Computes and returns the length of the burn-in phase using the values
	 * from {@link BurnInTest#setValue(double)}.
	 * 
	 * @return the length of the burn-in phase
	 */
	public int getLengthOfBurnIn();

	/**
	 * Returns a short description of the burn-in test.
	 * 
	 * @return a short description of the burn-in test.
	 */
	public String getInstanceName();
}
