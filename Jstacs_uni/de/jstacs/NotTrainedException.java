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

package de.jstacs;

/**
 * A {@link NotTrainedException} is thrown if the user tries to use an untrained
 * model.
 * 
 * @author Stefan Haufe, Jens Keilwagen
 */
public class NotTrainedException extends Exception {

	private static final long serialVersionUID = -1782795697441206482L;

	/**
	 * Creates a new {@link NotTrainedException} with standard error message
	 * (&quot;The model is not trained yet. Please try to train or load before
	 * an invocation of this method.&quot;).
	 */
	public NotTrainedException() {
		super( "The model is not trained yet. Please try to train or load before an invocation of this method." );
	}

	/**
	 * Creates a new {@link NotTrainedException} with given error
	 * <code>message</code>.
	 * 
	 * @param message
	 *            the error message
	 */
	public NotTrainedException( String message ) {
		super( message );
	}
}