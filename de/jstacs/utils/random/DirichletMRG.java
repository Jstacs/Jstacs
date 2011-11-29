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
 * This class is a multivariate random generator based on a Dirichlet
 * distribution.
 * 
 * @author Jens Keilwagen
 */
public class DirichletMRG extends MultivariateRandomGenerator {

	/**
	 * This instance shall be used, since quite often two instance of this class
	 * return the same values. With this a new Dirichlet random generator is
	 * created.
	 */
	public static final DirichletMRG DEFAULT_INSTANCE = new DirichletMRG();

	private RandomNumberGenerator r;

	private DirichletMRG() {
		r = new RandomNumberGenerator();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.random.MultivariateRandomGenerator#generate(double[], int, int, de.jstacs.utils.random.MRGParams)
	 */
	@Override
	public void generate( double[] d, int start, int n, MRGParams p ) {
		DiMRGParams param = (DiMRGParams)p;
		int i = param.getDimension();
		if( i > 0 && i != n ) {
			throw new IllegalArgumentException( "Hyperparameter doesnot have a correct dimension." );
		}
		double sum = 0d;
		for( i = 0; i < n; i++ ) {
			d[start + i] = r.nextGamma( param.getHyperparameter( i ), 1 );
			sum += d[start + i];
		}
		// we obtain: i == n;
		i += start;
		for( ; start < i; start++ ) {
			d[start] /= sum;
		}
	}
}
