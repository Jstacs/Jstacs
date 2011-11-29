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

/**
 * This class is the abstract super class for any multivariate random generator
 * (MRG). Some random generators need parameters. These will be given by an
 * object {@link MRGParams}.
 * 
 * @author Jens Keilwagen
 * 
 * @see MRGParams
 */
public abstract class MultivariateRandomGenerator {

	/**
	 * Generates a <code>n</code>-dimensional random array.
	 * 
	 * @param n
	 *            the dimension of the random array
	 * @param p
	 *            the parameter of the underlying distribution
	 * 
	 * @return a <code>n</code>-dimensional random array
	 * 
	 * @throws ClassCastException
	 *             if the object <code>p</code> could not be parsed to the
	 *             correct subclass
	 * @throws IllegalArgumentException
	 *             if an argument is not correct
	 * 
	 * @see MultivariateRandomGenerator#generate(double[], int, int, MRGParams)
	 */
	public double[] generate( int n, MRGParams p ) throws ClassCastException, IllegalArgumentException {
		double[] erg = new double[n];
		generate( erg, 0, n, p );
		return erg;
	}

	/**
	 * Generates a <code>n</code>-dimensional random array as part of the array
	 * <code>d</code> beginning at <code>start</code>.
	 * 
	 * @param d
	 *            the array
	 * @param start
	 *            the start index for generated values
	 * @param n
	 *            the dimension of the random array
	 * @param p
	 *            the parameter of the underlying distribution
	 * 
	 * @throws ClassCastException
	 *             if the object <code>p</code> could not be parsed to the
	 *             correct subclass
	 * @throws IllegalArgumentException
	 *             if an argument is not correct
	 */
	public abstract void generate( double[] d, int start, int n, MRGParams p ) throws ClassCastException, IllegalArgumentException;
}
