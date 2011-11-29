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
import java.util.Random;

/**
 * This random generator returns <code>1-epsilon</code> for one and equal parts
 * for the rest of a random vector. The content of container {@link MRGParams}
 * is never used in this class.
 * 
 * @author Jens Keilwagen
 */
public class SoftOneOfN extends MultivariateRandomGenerator {

	private Random r;

	private double epsilon;

	private double p;

	/**
	 * This constructor can be used for (soft) sampling one of n. One item will
	 * get <code>1-epsilon</code>, all the others will get an equal part.
	 * 
	 * @param epsilon
	 *            the value that should be subtracted from 1 for one item
	 * 
	 * @throws IllegalArgumentException
	 *             if the value of <code>epsilon</code> is not in [0,1]
	 */
	public SoftOneOfN( double epsilon ) throws IllegalArgumentException {
		r = new Random();
		if( epsilon < 0 || epsilon > 1 ) {
			throw new IllegalArgumentException( "The value of epsilon has to be in [0,1]" );
		}
		this.epsilon = epsilon;
		p = 1d - epsilon;
	}

	/**
	 * This constructor can be used for (hard) sampling one of n. One item will
	 * get 1, all the others 0.
	 */
	public SoftOneOfN() {
		this( 0 );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.random.MultivariateRandomGenerator#generate(double[], int, int, de.jstacs.utils.random.MRGParams)
	 */
	@Override
	public void generate( double[] d, int start, int number, MRGParams p ) {
		generate( d, start, number );
	}

	/**
	 * Generates an array of length <code>number</code> with one entry getting
	 * the value <code>1-epsilon</code> and all the others equal parts of
	 * <code>epsilon</code>.
	 * 
	 * @param number
	 *            the length of the generated array
	 * 
	 * @return array of length <code>number</code> with one entry getting the
	 *         value <code>1-epsilon</code> and all the others equal parts of
	 *         <code>epsilon</code>
	 * 
	 * @see SoftOneOfN#generate(double[], int, int)
	 */
	public double[] generate( int number ) {
		double[] erg = new double[number];
		generate( erg, 0, number );
		return erg;
	}

	/**
	 * Generates an array of length <code>number</code> as part of the array
	 * <code>d</code> beginning at index <code>start</code> with one entry
	 * getting the value <code>1-epsilon</code> and all the others equal parts
	 * of <code>epsilon</code>.
	 * 
	 * @param d
	 *            the array
	 * @param start
	 *            the start index for generated values
	 * @param number
	 *            the dimension of the generated array
	 */
	public void generate( double[] d, int start, int number ) {
		Arrays.fill( d, start, start + number, epsilon / (double)( number - 1 ) );
		d[start + r.nextInt( number )] = p;
	}
}
