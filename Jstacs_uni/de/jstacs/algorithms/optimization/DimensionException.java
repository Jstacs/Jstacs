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

package de.jstacs.algorithms.optimization;

/**
 * This class is for {@link Exception}s depending on wrong dimensions of vectors
 * for a given function.
 * 
 * @author Jens Keilwagen
 */
public class DimensionException extends Exception {

	private static final long serialVersionUID = 604134826262890722L;

	/**
	 * Creates a new {@link DimensionException} with standard error message
	 * (&quot;The vector has wrong dimension for this function.&quot;).
	 */
	public DimensionException() {
		super( "The vector has wrong dimension for this function." );
	}

	/**
	 * Creates a new {@link DimensionException} with a more detailed error
	 * message.
	 * 
	 * @param dimV
	 *            the dimension of the vector
	 * @param dimF
	 *            the dimension of the function
	 */
	public DimensionException( int dimV, int dimF ) {
		super( "The vector (dim = " + dimV + ") has wrong dimension for this function (dim = " + dimF + ")." );
	}
}
