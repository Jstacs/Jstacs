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

import Jama.EigenvalueDecomposition;
import Jama.Matrix;


/**
 * This class can be used to determine the stationary distribution.
 * 
 * @author Jens Keilwagen
 */
public class StationaryDistribution {

	/**
	 * This method return the stationary distribution. This distribution can be
	 * defined over a single symbol, pairs of symbols, ...
	 * 
	 * @param condProbs
	 *            an array containing all conditional probabilities;
	 *            <ul>
	 *            <li><code>condProbs[0] = P(a_1|b_1)</code>
	 *            <li><code>condProbs[1] = P(a_2|b_1)</code>
	 *            <li>...
	 *            </ul>
	 * @param alphabetSize
	 *            the number of symbols in the alphabet
	 * 
	 * @return the stationary distribution
	 */
	public static double[] getStationaryDistribution(double[] condProbs,
			int alphabetSize) {
		if (condProbs.length % alphabetSize != 0) {
			throw new IllegalArgumentException("wrong dimension");
		}
		int dim = condProbs.length / alphabetSize;
		// show( condProbs, alphabetSize );

		// create matrix
		double[][] matrix = new double[dim][dim];
		int help = dim / alphabetSize, h;
		for (int offset1, offset2, j, i = 0; i < dim; i++) {
			offset1 = i / alphabetSize;
			offset2 = i % alphabetSize;
			for (j = 0; j < alphabetSize; j++) {
				h = offset1 + j * help;
				matrix[i][h] = condProbs[h * alphabetSize + offset2];
			}
		}

		Matrix m = new Matrix(matrix, dim, dim);
		// show( m.getArray() );
		EigenvalueDecomposition eigen = new EigenvalueDecomposition(m);
		int ind = getIndex(eigen.getRealEigenvalues(), 1);
		Matrix m2 = eigen.getV().getMatrix(0, dim - 1, ind, ind).transpose();
		double[] res = m2.getArray()[0];
		// System.out.println();
		// System.out.println( Arrays.toString( eigen.getRealEigenvalues() ) );
		// System.out.println( "index " + getIndex( eigen.getRealEigenvalues(),
		// 1 ) );

		double sum = 0;
		for (int i = 0; i < dim; i++) {
			sum += res[i];
		}
		for (int i = 0; i < dim; i++) {
			res[i] /= sum;
		}
		/*
		 * //test System.out.println( "\ntest" ); show( m2.times( m.transpose()
		 * ).getArray() );
		 */

		return res;
	}

	private static int getIndex(double[] eigenValues, double wanted) {
		int i = 1, best = 0;
		double diff = Math.abs(1 - eigenValues[best]), current;
		for (; i < eigenValues.length; i++) {
			current = Math.abs(1 - eigenValues[i]);
			if (current < diff) {
				best = i;
				diff = current;
			}
		}
		return best;
	}

	/**
	 * This method returns the conditional stationary distributions, i.e. for a
	 * homogeneous Markov model of order 2 it returns an array containing the
	 * stationary symbol distribution as first entry, the conditional stationary
	 * distribution of order 1 as second entry and the conditional distribution
	 * of order 2 as third entry.
	 * 
	 * @param probs
	 *            an array containing all (unconditional) probabilities
	 * @param alphabetSize
	 *            the number of symbols in the alphabet
	 * 
	 * @return the conditional stationary distribution
	 * 
	 * @see StationaryDistribution#getStationaryDistribution(double[], int)
	 */
	public static double[][][] getAllConditionalStationaryDistributions(
			double[] probs, int alphabetSize) {
		double[] help = probs.clone();
		int dim = (int) (Math.log(probs.length) / Math.log(alphabetSize));
		double[][][] res = new double[dim][][];
		for (int i = dim - 1; i > 0; i--) {
			makeCondProb(help, alphabetSize);
			res[i] = fold(help, alphabetSize);
			help = getStationaryDistribution(help, alphabetSize);

			// System.out.println();
			// show( help, alphabetSize );
		}
		res[0] = new double[][] { help };
		return res;
	}

	private static void makeCondProb(double[] jointProb, int alphabetSize) {
		double sum;
		for (int j, i = 0; i < jointProb.length;) {
			sum = 0;
			for (j = 0; j < alphabetSize; j++) {
				sum += jointProb[i + j];
			}
			for (j = 0; j < alphabetSize; j++, i++) {
				// System.out.print( jointProb[i] );
				jointProb[i] /= sum;
				// System.out.println( "\t" + (sum*jointProb[i]) );
			}
		}
	}

	/**
	 * This method folds an one-dimensional array containing conditional
	 * probability distributions to a two-dimensional array <code>res</code>
	 * where each sub-array <code>res[i]</code> contains a conditional
	 * distribution.
	 * 
	 * @param prob
	 *            a one-dimensional array containing conditional probabilities
	 * @param alphabetSize
	 *            the number of symbols in the alphabet
	 * 
	 * @return a two-dimensional array <code>res</code> where each sub-array
	 *         <code>res[i]</code> contains a conditional distribution
	 */
	private static double[][] fold(double[] prob, int alphabetSize) {
		int l = prob.length / alphabetSize;
		double[][] res = new double[l][alphabetSize];
		for (int i = 0; i < l; i++) {
			System.arraycopy(prob, i * alphabetSize, res[i], 0, alphabetSize);
		}
		return res;
	}

	private static void show(double[] condProb, int alphabetSize) {
		for (int i = 0, j; i < condProb.length;) {
			for (j = 0; j < alphabetSize; j++, i++) {
				System.out.print(condProb[i] + "\t");
			}
			System.out.println();
		}
	}

	private static void show(double[][] matrix) {
		for (int i = 0, j; i < matrix.length; i++) {
			for (j = 0; j < matrix[i].length; j++) {
				System.out.print(matrix[i][j] + "\t");
			}
			System.out.println();
		}
	}
}