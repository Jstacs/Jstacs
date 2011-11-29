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
 * The super container for parameters of Dirichlet multivariate random
 * generators.
 * 
 * @author Jens Keilwagen
 * 
 * @see DirichletMRG
 */
public abstract class DiMRGParams extends MRGParams {

	/**
	 * Returns the dimension of the hyperparameter vector of the underlying
	 * Dirichlet distribution and therefore the dimension of the generated
	 * random array.
	 * 
	 * @return the dimension of the hyperparameter vector of the underlying
	 *         Dirichlet distribution and therefore the dimension of the
	 *         generated random array
	 */
	public abstract int getDimension();

	/**
	 * Returns the value at position <code>i</code> of the hyperparameter vector
	 * of the underlying Dirichlet distribution.
	 * 
	 * @param i
	 *            the position of the hyperparameter vector
	 * 
	 * @return the value at position <code>i</code> of the hyperparameter vector
	 *         of the underlying Dirichlet distribution
	 */
	public abstract double getHyperparameter( int i );
}
