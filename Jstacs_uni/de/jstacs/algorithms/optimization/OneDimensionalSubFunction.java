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
 * This class is used to do the line search.
 * 
 * @author Jens Keilwagen
 */
public class OneDimensionalSubFunction extends OneDimensionalFunction {

	private double[] d;

	private double[] current;

	private Function f;

	/**
	 * Creates a new {@link OneDimensionalSubFunction} from a {@link Function}
	 * <code>f</code> for the line search.
	 * 
	 * @param f
	 *            the high dimensional function
	 * @see OneDimensionalSubFunction#set(double[], double[])
	 */
	public OneDimensionalSubFunction( Function f ) {
		this.f = f;
		d = null;
		current = null;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.OneDimensionalFunction#evaluateFunction(double)
	 */
	@Override
	public double evaluateFunction( double x ) throws EvaluationException {
		double[] p = new double[d.length];
		for( int counter = 0; counter < p.length; counter++ ) {
			p[counter] = current[counter] + x * d[counter];
		}
		try {
			return f.evaluateFunction( p );
		} catch ( DimensionException impossible ) {
			EvaluationException ee = new EvaluationException( impossible.getMessage() );
			ee.setStackTrace( impossible.getStackTrace() );
			throw ee;
		}
	}
	
	/**
	 * Sets the current values and direction.
	 * 
	 * @param current
	 *            the current vector
	 * @param d
	 *            the direction along which the line search will be performed
	 * 
	 * @throws DimensionException
	 *             if there is something wrong with the dimension of a vector
	 */
	public void set( double[] current, double[] d ) throws DimensionException {
		int n = f.getDimensionOfScope();
		if( n != d.length && n != current.length ) {
			throw new DimensionException();
		}
		if( this.d == null || this.current.length != current.length ) {
			this.d = d.clone();
			this.current = current.clone();
		} else {
			System.arraycopy( current, 0, this.current, 0, this.current.length );
			System.arraycopy( d, 0, this.d, 0, this.d.length );
		}
	}
}