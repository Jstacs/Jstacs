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
 * The container for parameters of an Erlang multivariate random generator.
 * 
 * @author Jens Keilwagen
 * 
 * @see MRGParams
 */
public class ErlangMRGParams extends MRGParams {

	private int[] alpha;

	private int sum;

	/**
	 * Constructor which creates a new hyperparameter vector for an Erlang
	 * random generator. The hyperparameter vector of the underlying Erlang
	 * distribution consists of the same value at every position.
	 * 
	 * @param alpha
	 *            the value for the hyperparameter vector
	 * @param n
	 *            the dimension of the hyperparameter vector
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>n</code> is chosen incorrectly (has to be greater
	 *             than 2) or if <code>alpha</code> is chosen incorrectly (has
	 *             to be positive)
	 */
	public ErlangMRGParams( int alpha, int n ) throws IllegalArgumentException {
		if( n < 2 ) {
			throw new IllegalArgumentException( "The parameter n has to be at least 2." );
		}
		if( alpha <= 0d ) {
			throw new IllegalArgumentException( "The parameter alpha has to be positive." );
		}
		this.alpha = new int[n];
		Arrays.fill( this.alpha, alpha );
		sum = n * alpha;
	}

	/**
	 * Constructor which creates a new hyperparameter vector for an Erlang
	 * random generator.
	 * 
	 * @param alpha
	 *            the hyperparameter vector
	 * 
	 * @throws IllegalArgumentException
	 *             if at least one of the hyperparameters is not positive
	 */
	public ErlangMRGParams( int[] alpha ) throws IllegalArgumentException {
		this.alpha = new int[alpha.length];
		sum = 0;
		for( int i = 0; i < alpha.length; i++ ) {
			if( alpha[i] <= 0d ) {
				throw new IllegalArgumentException( "Each parameter alpha[i] has to be positive." );
			}
			this.alpha[i] = alpha[i];
			sum += alpha[i];
		}
	}

	/**
	 * Returns the dimension of the hyperparameter vector of the underlying
	 * Erlang distribution.
	 * 
	 * @return the dimension of the hyperparameter vector of the underlying
	 *         Erlang distribution
	 */
	public int getDimension() {
		return alpha.length;
	}

	/**
	 * Returns the value at position <code>i</code> of the hyperparameter vector
	 * of the underlying Erlang distribution.
	 * 
	 * @param i
	 *            the position of the hyperparameter vector
	 * 
	 * @return the value at position <code>i</code> of the hyperparameter vector
	 *         of the underlying Erlang distribution
	 */
	public int getHyperparameter( int i ) {
		return alpha[i];
	}

	/**
	 * Returns the sum of the hyperparameters (entries of the hyperparameter
	 * vector) of the underlying Erlang distribution.
	 * 
	 * @return the sum of the hyperparameters (entries of the hyperparameter
	 *         vector) of the underlying Erlang distribution
	 */
	public int getSumOfHyperparameter() {
		return sum;
	}
}
