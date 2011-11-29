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

package de.jstacs.motifDiscovery;


/**
 * This interface allows to modify a motif model.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public interface Mutable {
	/**
	 * Manually modifies the model. The two offsets <code>offsetLeft</code>
	 * and <code>offsetRight</code> define how many positions the left or
	 * right border positions shall be moved. Negative numbers indicate moves to
	 * the left while positive numbers correspond to moves to the right.
	 * 
	 * @param offsetLeft
	 *            the offset on the left side
	 * @param offsetRight
	 *            the offset on the right side
	 * 
	 * @return <code>true</code> if the motif model was modified otherwise
	 *         <code>false</code>
	 */
	public boolean modify( int offsetLeft, int offsetRight	);
}
