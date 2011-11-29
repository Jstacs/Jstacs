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
 * This interface is used to determine the next start distance that will be used
 * in a line search.
 * 
 * @author Jens Keilwagen
 * 
 * @see Optimizer
 * @see OneDimensionalSubFunction
 */
public interface StartDistanceForecaster {

	/**
	 * This method returns the new positive start distance. It enables the
	 * implementor to compute any interpolation using the last values.
	 * 
	 * @return the new positive start distance
	 * 
	 * @see StartDistanceForecaster#setLastDistance(double)
	 */
	public double getNewStartDistance();

	/**
	 * Sets the last used distance.
	 * 
	 * @param last
	 *            the last used distance
	 */
	public void setLastDistance( double last );

	/**
	 * Resets the object to the initial state.
	 */
	public void reset();
}
