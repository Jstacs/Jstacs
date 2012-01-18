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
 * A {@link WrongAlphabetException} is thrown if the given {@link de.jstacs.data.alphabets.Alphabet} or {@link de.jstacs.data.AlphabetContainer}
 * does not support some data.
 * 
 * @author Jan Grau
 */
public class WrongAlphabetException extends Exception {

	private static final long serialVersionUID = 3257008769530474547L;

	/**
	 * Creates a new {@link WrongAlphabetException} with standard error message
	 * (&quot;The data of the selected file does not match the entered
	 * alphabet.&quot;).
	 */
	public WrongAlphabetException() {
		super( "The data of the selected file does not match the entered alphabet." );
	}

	/**
	 * Creates a new {@link WrongAlphabetException} with given error
	 * <code>message</code>.
	 * 
	 * @param message
	 *            the error message
	 */
	public WrongAlphabetException( String message ) {
		super( message );
	}
}