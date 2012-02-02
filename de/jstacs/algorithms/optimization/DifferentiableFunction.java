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
 * This class is the framework for any (at least) one time differentiable
 * function {@latex.inline $f: \\mathbb{R}^n \\to \\mathbb{R}$}.
 * 
 * @author Jens Keilwagen
 */
public abstract class DifferentiableFunction implements Function {

	/**
	 * Evaluates the gradient of a function at a certain vector (in mathematical
	 * sense) <code>x</code>, i.e.,
	 * {@latex.inline $\\nabla f(\\underline{x}) = \\left(\\frac{\\partial f(\\underline{x})}{\\partial x_1},\\ldots,\\frac{\\partial f(\\underline{x})}{\\partial x_n}\\right)$}.
	 * 
	 * @param x
	 *            the current vector
	 * 
	 * @return the evaluation of the gradient of a function, has dimension
	 *         {@link Function#getDimensionOfScope()}
	 * 
	 * @throws DimensionException
	 *             if <code>dim(x) != n</code>, with {@latex.inline $f: \\mathbb{R}^n \\to \\mathbb{R}$}
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             gradient
	 * 
	 * @see Function#getDimensionOfScope()
	 */
	public abstract double[] evaluateGradientOfFunction( double[] x ) throws DimensionException, EvaluationException;

	/**
	 * This method is used to find an approximation of an one-dimensional
	 * subfunction. That means it will find an approximation of the minimum
	 * starting at point <code>x</code> and search in direction
	 * <code>d</code>, {@latex.ilb %preamble{\\usepackage{amsmath}} \\[\\operatorname{argmin}_{\\alpha
	 * \\ge 0} f(\\underline{x} + \\alpha\\underline{d})\\]}.
	 * 
	 * This method is a standard implementation. You are enabled to overwrite
	 * this method to be faster if you know anything about the problem or if you
	 * just want to test other line search methods.
	 * 
	 * @param x
	 *            the start point
	 * @param d
	 *            the search direction
	 * @param alpha_0
	 *            the initial alpha
	 * @param fAlpha_0
	 *            the initial function value (this value is known in most cases
	 *            and does not have to be computed again)
	 * @param linEps
	 *            the tolerance for stopping this method
	 * @param startDistance
	 *            the initial distance for bracketing the minimum
	 * 
	 * @return <code>double[2] res = { alpha*, f(alpha*) }</code>
	 * 
	 * @throws DimensionException
	 *             if there is something wrong with the dimension
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * 
	 * @see OneDimensionalFunction#findMin(double, double, double, double)
	 */
	protected double[] findOneDimensionalMin( double[] x, double[] d, double alpha_0, double fAlpha_0, double linEps,
			double startDistance ) throws DimensionException, EvaluationException {
		int i=0;
		while(i<d.length && d[i] == 0){
			i++;
		}
		OneDimensionalSubFunction fun = new OneDimensionalSubFunction( this, x, d );
		if(i==d.length){
			return new double[]{0,fun.evaluateFunction( 0 )};
		}
		return fun.findMin( alpha_0, fAlpha_0, linEps, startDistance );
	}
}
