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
 * The container for parameters of a Dirichlet random generator.
 * 
 * @author Jens Keilwagen
 * 
 * @see DiMRGParams
 * @see FastDirichletMRGParams
 */
public class DirichletMRGParams extends DiMRGParams {

	private double[] alpha;

	private double sum;

	/**
	 * Constructor which creates a new hyperparameter vector for a Dirichlet
	 * random generator. The hyperparameter vector of the underlying Dirichlet
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
	 * 
	 * @see FastDirichletMRGParams#FastDirichletMRGParams(double)
	 */
	public DirichletMRGParams( double alpha, int n ) throws IllegalArgumentException {
		if( n < 2 ) {
			throw new IllegalArgumentException( "The parameter n has to be at least 2." );
		}
		if( alpha <= 0d ) {
			throw new IllegalArgumentException( "The parameter alpha has to be positive." );
		}
		this.alpha = new double[n];
		Arrays.fill( this.alpha, alpha );
		sum = n * alpha;
	}

	/**
	 * Constructor which creates a new hyperparameter vector for a Dirichlet
	 * random generator.
	 * 
	 * @param alpha
	 *            the hyperparameter vector
	 *  
	 * @throws IllegalArgumentException
	 *             if at least one of the hyperparameters is not positive
	 */
	public DirichletMRGParams( double... alpha ) throws IllegalArgumentException {
		this( 0, alpha.length, alpha );
	}
	
	/**
	 * Constructor which creates a new hyperparameter vector for a Dirichlet
	 * random generator.
	 * 
	 * @param start start index (inclusive)
	 * @param end end index (exclusive)
	 * @param alpha
	 *            the hyperparameter vector
	 *  
	 * @throws IllegalArgumentException
	 *             if at least one of the hyperparameters is not positive
	 */
	public DirichletMRGParams( int start, int end, double... alpha ) throws IllegalArgumentException {
		this.alpha = new double[end-start];
		sum = 0;
		for( int i = 0; start < end; start++, i++ ) {
			if( alpha[start] <= 0d ) {
				throw new IllegalArgumentException( "Each parameter alpha["+start+"] has to be positive." );
			}
			this.alpha[i] = alpha[start];
			sum += this.alpha[i];
		}
	}
	
	

	/* (non-Javadoc)
	 * @see de.jstacs.utils.random.DiMRGParams#getDimension()
	 */
	@Override
	public int getDimension() {
		return alpha.length;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.random.DiMRGParams#getHyperparameter(int)
	 */
	@Override
	public double getHyperparameter( int i ) {
		return alpha[i];
	}

	/**
	 * Returns the sum of the hyperparameters (entries of the hyperparameter
	 * vector) of the underlying Dirichlet distribution.
	 * 
	 * @return the sum of the hyperparameters (entries of the hyperparameter
	 *         vector) of the underlying Dirichlet distribution
	 */
	public double getSumOfHyperparameter() {
		return sum;
	}
}
