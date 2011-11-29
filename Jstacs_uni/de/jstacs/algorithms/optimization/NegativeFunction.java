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
 * The negative function <i>-f</i> for a given {@link Function} <code>f</code>.
 * 
 * @author Jens Keilwagen
 */
public class NegativeFunction implements Function {

	private Function f;

	/**
	 * Creates the {@link Function} <code>f</code> for which <code>-f</code>
	 * should be calculated.
	 * 
	 * @param f
	 *            the {@link Function}
	 */
	public NegativeFunction( Function f ) {
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
}
