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

package de.jstacs.utils;

/**
 * Interface for supervising the progress of long time processes like cross
 * validation.
 * 
 * @author Andre Gohr, Jan Grau
 */
public interface ProgressUpdater {

	/**
	 * Sets the maximal value that will be set by {@link #setValue(int)}, so a
	 * value of max indicates the end of the supervised method call.
	 * 
	 * @param max
	 *            the maximal value
	 */
	public void setMax( int max );

	/**
	 * Sets the current value the supervised process has reached.
	 * 
	 * @param value
	 *            the current value
	 */
	public void setValue( int value );

	/**
	 * Specifies if the process is cancelled by the user.
	 * 
	 * @return <code>true</code> if process was cancelled, <code>false</code>
	 *         otherwise
	 */
	public boolean isCancelled();
}