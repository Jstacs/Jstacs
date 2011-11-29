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
 * The most simple {@link StartDistanceForecaster} that returns always the same
 * value.
 * 
 * @author Jens Keilwagen
 */
public class ConstantStartDistance implements StartDistanceForecaster {

	private double startDistance;

	/**
	 * This constructor creates an instance of {@link ConstantStartDistance}
	 * that returns always the given <code>value</code>.
	 * 
	 * @param value
	 *            the value that will always be returned
	 */
	public ConstantStartDistance( double value ) {
		if( value <= 0 ) {
			throw new IllegalArgumentException( "The start distance has to be positive." );
		}
		startDistance = value;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.StartDistanceForecaster#getNewStartDistance()
	 */
	public double getNewStartDistance() {
		return startDistance;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.StartDistanceForecaster#setLastDistance(double)
	 */
	public void setLastDistance( double last ) {}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.StartDistanceForecaster#reset()
	 */
	public void reset() {}
}
