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
 * The container for parameters of a Dirichlet random generator that uses the
 * same hyperparameter at all positions.
 * 
 * @author Jens Keilwagen
 * 
 * @see MRGParams
 * @see DirichletMRGParams
 */
public class FastDirichletMRGParams extends DiMRGParams {

	private double alpha;

	/**
	 * Creates the hyperparameter for a Dirichlet random generator which is used
	 * at all positions for the hyperparameter vector of the underlying
	 * Dirichlet distribution.
	 * 
	 * @param alpha
	 *            the hyperparameter for all positions of the hyperparameter
	 *            vector
	 * @throws IllegalArgumentException
	 *             if the value for <code>alpha</code> is chosen incorrectly
	 *             (negative or 0), the parameter has to be positive
	 */
	public FastDirichletMRGParams( double alpha ) throws IllegalArgumentException {
		if( alpha <= 0d ) {
			throw new IllegalArgumentException( "The parameter alpha has to be positive." );
		}
		this.alpha = alpha;
	}

	/* (non-Javadoc)
	 * Since there is no need to save the hyperparameters for this case in an array this method returns -1 as dimension.
	 * @see de.jstacs.utils.random.DiMRGParams#getDimension()
	 */
	@Override
	public int getDimension() {
		return -1;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.random.DiMRGParams#getHyperparameter(int)
	 */
	@Override
	public double getHyperparameter( int i ) {
		return alpha;
	}
}
