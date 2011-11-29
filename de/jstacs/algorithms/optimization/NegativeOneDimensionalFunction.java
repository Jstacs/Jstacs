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
 * This class extends the class {@link OneDimensionalFunction}. That class is
 * useful when you want to find a maximum of any one-dimensional function. Just
 * use this class and find the minimum of your new function. This minimum has
 * the same abscissa and just a negative ordinate.
 * 
 * @author Jens Keilwagen
 */
public class NegativeOneDimensionalFunction extends OneDimensionalFunction {

	private OneDimensionalFunction f;

	/**
	 * Creates the {@link OneDimensionalFunction} <code>f</code> for which
	 * <code>-f</code> should be calculated.
	 * 
	 * @param f
	 *            the {@link OneDimensionalFunction}
	 */
	public NegativeOneDimensionalFunction( OneDimensionalFunction f ) {
		this.f = f;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.OneDimensionalFunction#evaluateFunction(double[])
	 */
	@Override
	public double evaluateFunction( double[] x ) throws DimensionException, EvaluationException {
		if( x == null || x.length != 1 ) {
			throw new DimensionException();
		}
		return evaluateFunction( x[0] );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.OneDimensionalFunction#evaluateFunction(double)
	 */
	@Override
	public double evaluateFunction( double x ) throws EvaluationException {
		return -f.evaluateFunction( x );
	}
}
