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
 * Simple class that implements {@link ProgressUpdater} and prints the
 * percentage of iterations that is already done on the screen.
 * 
 * @author Andre Gohr, Jan Grau, Jens Keilwagen
 */
public class DefaultProgressUpdater implements ProgressUpdater {

	/**
	 * The maximal number of steps.
	 */
	protected int max;

	private String old;

	/**
	 * Creates a {@link DefaultProgressUpdater}.
	 */
	public DefaultProgressUpdater() {
		this.old = "";
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.ProgressUpdater#setMax(int)
	 */
	public void setMax( int max ) {
		this.max = max;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.ProgressUpdater#setValue(int)
	 */
	public void setValue( int value ) {
		String str = Math.round( value * 100d / max ) + "%";

		if( !str.equals( old ) ) {
			for( int i = 0; i < old.length(); i++ )
				System.err.print( "\b" );
			System.err.print( str );
			old = str;
		}
		if( value == max ) System.err.println();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.ProgressUpdater#isCancelled()
	 */
	public boolean isCancelled() {
		return false;
	}
}