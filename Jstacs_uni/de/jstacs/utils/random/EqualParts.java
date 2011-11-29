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

package de.jstacs.utils.random;

import java.util.Arrays;

/**
 * This class is no real random generator it just returns 1/n for all values.
 * The content of container {@link MRGParams} is never used in this class.
 * 
 * @author Jens Keilwagen
 */
public class EqualParts extends MultivariateRandomGenerator {

	/* (non-Javadoc)
	 * @see de.jstacs.utils.random.MultivariateRandomGenerator#generate(double[], int, int, de.jstacs.utils.random.MRGParams)
	 * @see EqualParts#generate(double[], int, int)
	 */
	@Override
	public void generate( double[] d, int start, int number, MRGParams p ) {
		generate( d, start, number );
	}

	/**
	 * Returns an array of length <code>n</code> with entries <code>1/n</code>.
	 * 
	 * @param n
	 *            the length of the array and the divisor for the calculation of
	 *            the equal parts as its entries
	 * 
	 * @return an array of length <code>n</code> with entries 1/n
	 * 
	 * @see EqualParts#generate(double[], int, int)
	 */
	public double[] generate( int n ) {
		double[] erg = new double[n];
		generate( erg, 0, n );
		return erg;
	}

	/**
	 * Modifies a given <code>number</code> of array entries to
	 * <code>1/number</code>.
	 * 
	 * @param d
	 *            the array which will be modified
	 * @param start
	 *            the start index of the modification of the array
	 * @param number
	 *            the number of entries <code>1/number</code>
	 */
	public void generate( double[] d, int start, int number ) {
		Arrays.fill( d, start, start + number, 1d / (double)number );
	}
}
