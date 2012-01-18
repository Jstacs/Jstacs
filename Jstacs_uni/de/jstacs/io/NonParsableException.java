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

package de.jstacs.io;

/**
 * A {@link NonParsableException} is thrown if some object could not be restored
 * (parsed) from a {@link StringBuffer}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class NonParsableException extends Exception {

	private static final long serialVersionUID = 3257284725474670903L;

	/**
	 * Creates a new {@link NonParsableException} with standard error message
	 * (&quot;StringBuffer not parsable.&quot;).
	 */
	public NonParsableException() {
		super( "StringBuffer not parsable." );
	}

	/**
	 * Creates a new {@link NonParsableException} with given error
	 * <code>message</code>.
	 * 
	 * @param message
	 *            the error message
	 */
	public NonParsableException( String message ) {
		super( "StringBuffer not parsable. (" + message + ")" );
	}
}