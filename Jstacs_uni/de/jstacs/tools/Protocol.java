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

package de.jstacs.tools;

import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.tools.ui.galaxy.Galaxy;

/**
 * Interface for the protocol of the run of a {@link JstacsTool}.
 * Within a {@link JstacsTool}, only this generic interface should be used.
 * The implementation of the methods highly depends on the interface, e.g., {@link CLI} or {@link Galaxy}.
 * 
 * @author Jan Grau
 *
 */
public interface Protocol {

	/**
	 * Appends a general message to the protocol
	 * @param str the message
	 */
	public void append(String str);
	
	/**
	 * Appends a heading to the protocol, which is generally highlighted in some appropriate way.
	 * @param heading the heading
	 */
	public void appendHeading(String heading);
	
	/**
	 * Appends a warning to the protocol, which is generally highlighted in some appropriate way.
	 * @param warning the warning
	 */
	public void appendWarning(String warning);
	
	/**
	 * Appends a {@link Throwable} to the protocol. The message of the protocol is generally highlighted in some appropriate way.
	 * Depending on the implementation, the stack trace of the {@link Throwable} may also be displayed.
	 * @param th the {@link Throwable}
	 * @see Throwable#printStackTrace()
	 */
	public void appendThrowable(Throwable th);
	
	/**
	 * Appends some verbatim text (i.e., text that is displayed "as is" regardless of default formatting) to the protocol.
	 * @param verbatim the verbatim text
	 */
	public void appendVerbatim(String verbatim);
}
