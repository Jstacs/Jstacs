/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.utils;

/**
 * This class enables the user to compute some divergences.
 * 
 * In Latex notation the divergence is
 * <p>
 * <code>
 * Y_{\epsilon}(p,q) = \frac{\left[\sum_{i,j} p_{i,j} \cdot \left(\frac{p_{i,j}}{q_{i,j}}\right)^{\epsilon-1}\right] - 1}{\frac{\epsilon(\epsilon-1)}{2}}
 * </code>
 * </p>
 * 
 * <b>NOTE:</b> There is no checking if the user given matrices and vectors are
 * stochastic.
 * 
 * @author Jens Keilwagen
 */
public class StatisticalTest {
	/**
	 * Computes the generalized divergence for two given stochastic matrices
	 * over the same domain, i.e. the matrices have to have the same
	 * dimensionality.
	 * 
	 * <br>
	 * <br>
	 * 
	 * For <code>epsilon=1</code> this method returns 2*(mutual information).
	 * 
	 * @param p
	 *            one stochastic matrix
	 * @param q
	 *            another stochastic matrix
	 * @param epsilon
	 *            the positive divergence parameter
	 * 
	 * @return the value of the generalized divergence for the two matrices
	 * 
	 * @throws IllegalArgumentException
	 *             if some arguments are not correct
	 */
	public static double getGeneralizedDivergence(double[][] p, double[][] q,
			double epsilon) throws IllegalArgumentException {
		if (epsilon <= 0) {
			throw new IllegalArgumentException("epsilon has to be positive");
		}
		int i = 0, j, length1 = p.length, length2 = p[0].length;
		if (q.length != length1) {
			throw new IllegalArgumentException(
					"The matrices have to have same dimension (" + length1
							+ " x " + length2 + ").");
		}
		double res = 0;
		if (epsilon == 1d) {
			for (; i < length1; i++) {
				if (p[i].length != length2 || q[i].length != length2) {
					throw new IllegalArgumentException(
							"The matrices have to have same dimension ("
									+ length1 + " x " + length2 + ").");
				}
				for (j = 0; j < length2; j++) {
					if (p[i][j] != 0) // lim_{p_{i,j} \to 0} p_{i,j} \ln p_{i,j}
										// = 0
					{
						res += p[i][j] * Math.log(p[i][j] / q[i][j]);
					}
				}
			}
			res *= 2d;
		} else {
			double e = epsilon - 1d;
			for (; i < length1; i++) {
				if (p[i].length != length2 || q[i].length != length2) {
					throw new IllegalArgumentException(
							"The matrices have to have same dimension ("
									+ length1 + " x " + length2 + ").");
				}
				for (j = 0; j < length2; j++) {
					res += p[i][j] * Math.pow(p[i][j] / q[i][j], e);
				}
			}
			res = (res - 1d) / (epsilon * e / 2d);
		}
		return res;
	}

	/**
	 * Computes the generalized divergence for two stochastic matrices over the
	 * same domain, i.e. the matrices have to have the same dimensionality. The
	 * second matrix is computed as joint probability of the given stochastic
	 * vectors.
	 * 
	 * <br>
	 * <br>
	 * 
	 * For <code>epsilon=1</code> this method returns 2*(mutual information).
	 * 
	 * @param p
	 *            a stochastic matrix
	 * @param r
	 *            one stochastic vector
	 * @param s
	 *            another stochastic vector
	 * @param epsilon
	 *            the divergence parameter
	 * 
	 * @return the value of the generalized divergence for two stochastic
	 *         matrices
	 */
	public static double getGeneralizedDivergence(double[][] p, double[] r,
			double[] s, double epsilon) {
		return getGeneralizedDivergence(p, getJointDistribution(r, s), epsilon);
	}

	/**
	 * Computes the generalized divergence for two stochastic matrices over the
	 * same domain, i.e. the matrices have to have the same dimensionality. The
	 * second matrix is computed as joint probability of the marginal
	 * distributions.
	 * 
	 * <br>
	 * <br>
	 * 
	 * For <code>epsilon=1</code> this method returns 2*(mutual information).
	 * 
	 * @param p
	 *            a stochastic matrix
	 * @param epsilon
	 *            the divergence parameter
	 * 
	 * @return the value of the generalized divergence for two stochastic
	 *         matrices
	 */
	public static double getGeneralizedDivergence(double[][] p, double epsilon) {
		double[][] marginal = getMarginalDistributions(p);
		return getGeneralizedDivergence(p, marginal[0], marginal[1], epsilon);
	}

	/**
	 * Returns the marginal distributions for a given stochastic matrix.
	 * 
	 * @param p
	 *            a stochastic matrix
	 * 
	 * @return the marginals
	 *         <ol>
	 *         <li><code>marginal[0][i] = \sum_j p[i][j]</code></li>
	 *         <li><code>marginal[1][j] = \sum_i p[i][j]</code></li>
	 *         </ol>
	 */
	private static double[][] getMarginalDistributions(double[][] p) {
		double[][] marginal = new double[2][];
		marginal[0] = new double[p.length];
		marginal[1] = new double[p[0].length];
		for (int i = 0, j; i < marginal[0].length; i++) {
			for (j = 0; j < marginal[1].length; j++) {
				marginal[0][i] += p[i][j];
				marginal[1][j] += p[i][j];
			}
		}
		return marginal;
	}

	/**
	 * Returns a matrix containing the joint distribution.
	 * 
	 * @param r
	 *            a stochastic vector
	 * @param s
	 *            a stochastic vector
	 * 
	 * @return a matrix <code>q</code> containing
	 *         <code>q[i][j] = r[i]*s[j]</code>
	 */
	private static double[][] getJointDistribution(double[] r, double[] s) {
		double[][] q = new double[r.length][s.length];
		for (int i = 0, j; i < r.length; i++) {
			for (j = 0; j < s.length; j++) {
				q[i][j] = r[i] * s[j];
			}
		}
		return q;
	}
}
