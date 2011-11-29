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

package de.jstacs.parameters;

/**
 * Interface for {@link Parameter}s that potentially can be varied over a
 * certain interval or a list of values. If {@link #isRangeable()} returns
 * <code>true</code> an instance of a {@link RangeIterator} that is defined
 * equivalently (data type, name, comment, etc.) to the current instance but can
 * have multiple values can be obtained using the method
 * {@link #getRangedInstance()}.
 * 
 * @author Jan Grau
 * 
 * @see RangeIterator
 */
public interface Rangeable {

	/**
	 * Returns <code>true</code> if the parameters can be varied over a range of
	 * values.
	 * 
	 * @return <code>true</code> if the parameter can be varied,
	 *         <code>false</code> otherwise
	 */
	public boolean isRangeable();

	/**
	 * Returns an instance of {@link RangeIterator} that has the same properties
	 * as the current instance, but accepts a range or list of values. These
	 * values can be obtained by the methods of {@link RangeIterator}.
	 * 
	 * @return an instance with the same properties as the current instance
	 * 
	 * @throws Exception
	 *             is thrown if {@link #isRangeable()} returns
	 *             <code>false</code> or the ranged instance could not be
	 *             created for some other reason
	 */
	public Parameter getRangedInstance() throws Exception;
}
