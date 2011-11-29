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

import java.util.Random;

/**
 * This class is a multivariate random generator based on a Dirichlet
 * distribution for <code>alpha_i \in N</code>. It is a kind of
 * &quot;multivariate Erlang distribution&quot;.
 * 
 * @author Jens Keilwagen
 */
public class ErlangMRG extends MultivariateRandomGenerator {

	private Random r;

	/**
	 * Constructor that creates a new multivariate random generator with
	 * underlying Erlang distribution.
	 */
	public ErlangMRG() {
		r = new Random();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.random.MultivariateRandomGenerator#generate(double[], int, int, de.jstacs.utils.random.MRGParams)
	 */
	@Override
	public void generate( double[] d, int start, int n, MRGParams p ) {
		ErlangMRGParams param = (ErlangMRGParams)p;
		if( param.getDimension() != n ) {
			throw new IllegalArgumentException( "Hyperparameter doesnot have a correct dimension." );
		}
		double sum = 0d;
		int i;
		for( i = 0; i < n; i++ ) {
			d[start + i] = erlangDistributed( param.getHyperparameter( i ), r );
			sum += d[start + i];
		}
		//we obtain: i == n;
		i += start;
		for( ; start < i; start++ ) {
			d[start] /= sum;
		}
	}

	private static double erlangDistributed( int alpha, Random r ) {
		double erg = 0;
		for( int i = 0; i < alpha; i++ ) {
			erg -= Math.log( 1 - r.nextDouble() );
		}
		return erg;
	}
}
