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
 * This class is for an {@link Exception} that is thrown if something with a
 * termination was not correct.
 * 
 * @author Jens Keilwagen
 */
public class TerminationException extends Exception {

	private static final long serialVersionUID = -2379596737222969758L;

	/**
	 * Creates a new {@link TerminationException} with standard error message
	 * (&quot;The termination mode was incorrect, please check your
	 * choice.&quot;).
	 */
	public TerminationException() {
		super( "The termination mode was incorrect, please check your choice." );
	}
}
