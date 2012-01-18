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

package de.jstacs.classifiers;

/**
 * This class indicates that a classifier is defined for less than 2 classes or
 * is defined over a different number of classes than given (e.g. by an array of
 * {@link de.jstacs.data.DataSet}s).
 * 
 * @author Jens Keilwagen
 */
public class ClassDimensionException extends Exception {

	private static final long serialVersionUID = 1L;

	/**
	 * This constructor creates a {@link ClassDimensionException} with the
	 * default error message (&quot;The number of classes in the classfier
	 * differs from the given number.&quot;).
	 */
	public ClassDimensionException() {
		super( "The number of classes in the classfier differs from the given number." );
	}

	/**
	 * This constructor creates a {@link ClassDimensionException} with given
	 * error message.
	 * 
	 * @param message
	 *            the error message
	 */
	public ClassDimensionException( String message ) {
		super( message );
	}
}
