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

package de.jstacs.data.alphabets;

/**
 * A {@link DoubleSymbolException} is thrown if a symbol occurred more than once
 * in an alphabet.
 * 
 * @author Jan Grau
 */
public class DoubleSymbolException extends Exception {

	private static final long serialVersionUID = 3258132444627677750L;

	/**
	 * Constructor for a {@link DoubleSymbolException} that takes the symbol
	 * that occurs more than once in the error message.
	 * 
	 * @param symbol
	 *            the symbol that occurred more than once
	 * 
	 * @see Exception#Exception(String)
	 */
	public DoubleSymbolException( String symbol ) {

		super( "Alphabet contained symbol \"" + symbol + "\" at least twice!" );
	}
}
