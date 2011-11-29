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
 * The negative function <code>-f</code> for a given
 * {@link DifferentiableFunction} <code>f</code>.
 * 
 * @author Jens Keilwagen
 */
public class NegativeDifferentiableFunction extends DifferentiableFunction {

	private DifferentiableFunction f;

	/**
	 * Creates the {@link DifferentiableFunction} <code>f</code> for which
	 * <code>-f</code> should be calculated.
	 * 
	 * @param f
	 *            the {@link DifferentiableFunction}
	 */
	public NegativeDifferentiableFunction( DifferentiableFunction f ) {
		this.f = f;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#evaluateFunction(double[])
	 */
	public double evaluateFunction( double[] x ) throws DimensionException, EvaluationException {
		return -f.evaluateFunction( x );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#getDimensionOfScope()
	 */
	public int getDimensionOfScope() {
		return f.getDimensionOfScope();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.DifferentiableFunction#evaluateGradientOfFunction(double[])
	 */
	@Override
	public double[] evaluateGradientOfFunction( double[] x ) throws DimensionException, EvaluationException {
		double[] erg = f.evaluateGradientOfFunction( x );
		for( int counter = 0; counter < erg.length; counter++ ) {
			erg[counter] *= -1d;
		}
		return erg;
	}

}
