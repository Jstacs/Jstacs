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
 * Class for an exception that is thrown if some error occurs while setting a
 * parameter's value or constructing a parameter.
 * 
 * @author Jan Grau
 */
public class ParameterException extends Exception {

	private static final long serialVersionUID = -5051644692781763386L;

	/**
	 * Constructor for a {@link ParameterException} with the specified error
	 * message.
	 * 
	 * @param message
	 *            the error message
	 */
	public ParameterException(String message) {
		super("Error in parameter: " + message);
	}

	/**
	 * Constructor for a {@link ParameterException} without a specific message.
	 */
	public ParameterException() {
		this("");
	}

}
