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

import java.util.Arrays;

/**
 * This class implements a {@link StartDistanceForecaster} that uses the
 * median of a limited memory over the last values. The method
 * {@link LimitedMedianStartDistance#getNewStartDistance()} returns
 * <code>0.667*median</code>.
 * 
 * @author Jens Keilwagen
 */
public class LimitedMedianStartDistance implements StartDistanceForecaster {

	private double initial;

	private double[] memory, sorted;

	private int index;

	/**
	 * This constructor creates an instance with <code>slots</code> memory slots
	 * that will initially be filled with <code>value</code>.
	 * 
	 * @param slots
	 *            the number of slots in the memory
	 * @param value
	 *            the initial value for the slots
	 */
	public LimitedMedianStartDistance( int slots, double value ) {
		if( value <= 0 ) {
			throw new IllegalArgumentException( "The initial start distance has to be positive." );
		}
		initial = value;
		if( slots < 1 ) {
			throw new IllegalArgumentException( "The number of slots has to be positive." );
		}
		memory = new double[slots];
		sorted = new double[slots];
		reset();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.StartDistanceForecaster#getNewStartDistance()
	 */
	public double getNewStartDistance() {
		// naive implementation
		System.arraycopy( memory, 0, sorted, 0, memory.length );
		Arrays.sort( sorted );
		double median;
		if( memory.length % 2 == 1 ) {
			median = sorted[memory.length / 2];
		} else {
			median = 0.5d * ( sorted[memory.length / 2 - 1] + sorted[memory.length / 2] );
		}
		return 0.667d*median;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.StartDistanceForecaster#setLastDistance(double)
	 */
	public void setLastDistance( double last ) {
		memory[index] = last;
		index++;
		index %= memory.length;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.StartDistanceForecaster#reset()
	 */
	public void reset() {
		index = 0;
		Arrays.fill( memory, initial );
	}

}
