/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.parameters.validation;

import de.jstacs.Storable;

/**
 * Interface for a parameter validator, i.e. a class that can validate some
 * possible parameter values.
 * 
 * @author Jan Grau
 */
public interface ParameterValidator extends Storable, Cloneable {

	/**
	 * Returns <code>true</code> if the value is valid and <code>false</code>
	 * otherwise.
	 * 
	 * @param value
	 *            the {@link Object} to be checked
	 * 
	 * @return if <code>value</code> is valid
	 */
	public boolean checkValue(Object value);

	/**
	 * Returns the error message if {@link #checkValue(Object)} returned false.
	 * 
	 * @return the error message
	 */
	public String getErrorMessage();

	/**
	 * This method returns a deep copy of the current instance.
	 * 
	 * @return a deep copy of the current index
	 * 
	 * @throws CloneNotSupportedException
	 *             if the {@link ParameterValidator} could not be cloned
	 * 
	 * @see Cloneable
	 */
	public ParameterValidator clone() throws CloneNotSupportedException;

}
