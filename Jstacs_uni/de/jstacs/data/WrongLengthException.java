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
 * A {@link WrongLengthException} is thrown if a given (sub)sequence length is
 * not correct.
 * 
 * @author Jan Grau
 */
public final class WrongLengthException extends Exception {

	private static final long serialVersionUID = 3257853172968927799L;

	/**
	 * Simple Constructor.
	 * 
	 * @param length
	 *            the length that was not supported
	 */
	public WrongLengthException( int length ) {
		super( "The subsequence length of " + length + " is not supported by this sample set." );
	}

	/**
	 * Creates a new {@link WrongLengthException} with your own
	 * <code>message</code>.
	 * 
	 * @param message
	 *            the message
	 */
	public WrongLengthException( String message ) {
		super( message );
	}
}