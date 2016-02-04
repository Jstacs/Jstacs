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

import java.util.Comparator;

/**
 * This class implements a {@link Comparator} of double arrays.
 * It just uses the double value at a specific index of the array to compare two arrays.
 * This index can be set in the constructor.
 * 
 * @author Jens Keilwagen
 */
public class DoubleArrayComparator implements Comparator<double[]> {
	
	private int index;

	/**
	 * Creates a new instance with a specific index for comparing double arrays.
	 * 
	 * @param idx the index of the element which should be used to compare two arrays
	 */
	public DoubleArrayComparator( int idx ) {
		index = idx;
	}
	
	@Override
	public int compare(double[] o1, double[] o2) {
		return (int) Math.signum( o1[index] - o2[index] );
	}		
}