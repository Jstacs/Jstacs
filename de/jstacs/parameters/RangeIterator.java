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
 * Interface for a {@link Parameter} or {@link ParameterSet} that can have
 * multiple values at the same time. This includes methods to loop over all
 * values while preserving the expected behaviour of the methods defined in
 * {@link Parameter} respectively {@link ParameterSet}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public interface RangeIterator {
	/**
	 * Switches to the next value in the collection of values in the specified
	 * range.
	 * 
	 * @return <code>true</code> if the next element exists, <code>false</code>
	 *         otherwise
	 * 
	 * @throws ParameterException
	 *             if the next value could not be fetched
	 */
	public boolean next() throws ParameterException;

	/**
	 * Resets the current value in the collection to the first value.
	 */
	public void resetToFirst();

	/**
	 * Returns the number of values in the collection.
	 * 
	 * @return the number of values
	 */
	public int getNumberOfValues();

	/**
	 * Returns a {@link String} representation of the set of values.
	 * 
	 * @return the {@link String} representation
	 */
	public String valuesToString();

	/**
	 * Returns <code>true</code> if this {@link RangeIterator} is ranging over a
	 * set of values.
	 * 
	 * @return <code>true</code> if this {@link RangeIterator} is ranging over a
	 *         set of values, <code>false</code> otherwise
	 */
	public boolean isRanged();
}
