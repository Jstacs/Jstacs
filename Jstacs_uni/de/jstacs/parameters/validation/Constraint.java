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
 * Interface for a constraint that must be fulfilled in a
 * {@link ConstraintValidator}.
 * 
 * @author Jan Grau
 */
public interface Constraint extends Storable, Cloneable {

	/**
	 * The condition is equality
	 */
	public static final int EQUALS = 1;
	/**
	 * The condition is less than
	 */
	public static final int LT = 2;
	/**
	 * The condition is greater than
	 */
	public static final int GT = 3;
	/**
	 * The condition is less or equal
	 */
	public static final int LEQ = 4;
	/**
	 * The condition is greater or equal
	 */
	public static final int GEQ = 5;

	/**
	 * Checks <code>value</code> for the constraint defined in the
	 * {@link Constraint}. If the constraint is fulfilled, <code>true</code> is
	 * returned.
	 * 
	 * @param value
	 *            the {@link Object} to be checked
	 * 
	 * @return <code>true</code> if the constraint is fulfilled,
	 *         <code>false</code> otherwise
	 */
	public boolean check(Object value);

	/**
	 * Returns the message of the last error (missed constraint) or
	 * <code>null</code> if the constraint was fulfilled by the last checked
	 * value.
	 * 
	 * @return the message of the last error or <code>null</code>
	 */
	public String getErrorMessage();

	/**
	 * This method returns a deep copy of the current instance.
	 * 
	 * @return a deep copy of the current instance
	 * 
	 * @throws CloneNotSupportedException
	 *             if the {@link Constraint} could not be cloned
	 * 
	 * @see Cloneable
	 */
	public Constraint clone() throws CloneNotSupportedException;

}
