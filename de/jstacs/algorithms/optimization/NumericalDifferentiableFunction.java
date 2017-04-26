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

package de.jstacs.algorithms.optimization;

/**
 * This class is the framework for any numerical differentiable function
 *{@latex.inline $f: \\mathbb{R}^n \\to \\mathbb{R}$}. The gradient is computed numerically. Each partial
 * differentiation is computed by the formula<br>
 * 
 * {@latex.inline $\\partial_k f(x) = \\frac{f(x)-f(x+\\varepsilon*e_k)}{\\varepsilon}$}.
 * 
 * @author Jens Keilwagen
 */
public class NumericalDifferentiableFunction extends DifferentiableFunction {

	/**
	 * The constant used in the computation of the gradient. Should be close to
	 * 0 but not exactly 0.
	 */
	protected double eps;
	
	protected Function f;

	/**
	 * Sets the function and value for epsilon for this
	 * {@link NumericalDifferentiableFunction}.
	 * 
	 * @param f
	 *            the function to be used 			
	 * @param epsilon
	 *            the epsilon used for the numerical differentiation
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>epsilon = 0</code>
	 */
	public NumericalDifferentiableFunction( Function f, double epsilon ) throws IllegalArgumentException {
		if( epsilon == 0 ) {
			throw new IllegalArgumentException( "Epsilon can not be 0." );
		}
		this.f = f;
		eps = epsilon;
	}
	
	/**
	 * Evaluates the gradient of a function at a certain vector (in mathematical
	 * sense) <code>x</code> numerically.
	 * 
	 * @see DifferentiableFunction#evaluateGradientOfFunction(double[])
	 */
	@Override
	public double[] evaluateGradientOfFunction( double[] x ) throws DimensionException, EvaluationException {
		int n = getDimensionOfScope();
		if( x == null || x.length != getDimensionOfScope() ) {
			if( x != null ) {
				throw new DimensionException( x.length, n );
			} else {
				throw new DimensionException( 0, n );
			}
		}
		double[] gradient = new double[n];
		double current = evaluateFunction( x ), h;
		for( int i = 0; i < n; i++ ) {
			h = x[i];
			x[i] += eps;
			gradient[i] = ( evaluateFunction( x ) - current ) / eps;
			x[i] = h;
		}

		return gradient;
	}

	@Override
	public double evaluateFunction(double[] x) throws DimensionException, EvaluationException {
		return f.evaluateFunction(x);
	}

	@Override
	public int getDimensionOfScope() {
		return f.getDimensionOfScope();
	}
}