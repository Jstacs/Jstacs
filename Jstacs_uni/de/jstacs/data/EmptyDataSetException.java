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

package de.jstacs.data;


/**
 * An {@link EmptyDataSetException} will be thrown if no {@link de.jstacs.data.sequences.Sequence} is in a
 * {@link DataSet} (i.e. the {@link DataSet} is empty).
 * 
 * @author Jens Keilwagen
 */
public class EmptyDataSetException extends Exception {

	private static final long serialVersionUID = 1L;

	/**
	 * This constructor creates an instance with default error message
	 * (&quot;The created DataSet is empty.&quot;).
	 */
	public EmptyDataSetException() {
		super( "The created DataSet is empty." );
	}
}
