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
 * This class implements the interface {@link Function} for an one-dimensional
 * function.
 * 
 * @author Jens Keilwagen
 */
public abstract class OneDimensionalFunction implements Function {

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#evaluateFunction(double[])
	 */
	public double evaluateFunction( double[] x ) throws DimensionException, EvaluationException {
		if( x == null || x.length != 1 ) {
			throw new DimensionException();
		}
		return evaluateFunction( x[0] );
	}

	/**
	 * Evaluates the function at position <code>x</code>.
	 * 
	 * @param x
	 *            the current position
	 * 
	 * @return <code>f(x)</code>
	 * 
	 * @throws EvaluationException
	 *             if there was a mistake during evaluation of the function
	 */
	public abstract double evaluateFunction( double x ) throws EvaluationException;

	/**
	 * This method returns a minimum <code>x</code> and the value
	 * <code>f(x)</code>, starting the search at <code>lower</code>. First the
	 * minimum is tried to be bracketed and then the interval is tried to be
	 * reduced.<br>
	 * 
	 * This method is a standard implementation and can be overwritten any time.
	 * It uses methods from the class {@link Optimizer}.
	 * 
	 * @param lower
	 *            the initial value of <code>x</code>
	 * @param fLower
	 *            the value <code>f(x)</code>
	 * @param eps
	 *            the threshold to stop the search of the minimum
	 * @param startDistance
	 *            the initial distance for bracketing the minimum
	 * 
	 * @return an array containing <code>x</code> at position 0 and
	 *         <code>f(x)</code> at position 1
	 * 
	 * @throws EvaluationException
	 *             if there was a mistake during the evaluation of the function
	 * 
	 * @see Optimizer
	 */
	public double[] findMin( double lower, double fLower, double eps, double startDistance ) throws EvaluationException {
		//necessary for heuristic methods, additional runtime should be negligible
		fLower = this.evaluateFunction(lower);//TODO
		double[] bracket = Optimizer.findBracket( this, lower, fLower, startDistance );
		return Optimizer.brentsMethod( this, bracket[0], bracket[2], bracket[3], bracket[4], eps );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#getDimensionOfScope()
	 */
	public final int getDimensionOfScope() {
		return 1;
	}
}
