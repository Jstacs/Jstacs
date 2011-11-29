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

import java.io.IOException;

/**
 * This class implements a quadratic function.
 * 
 * @author Jens Keilwagen
 */
public class QuadraticFunction extends OneDimensionalFunction {

	private double a, b, c;

	/**
	 * This constructor creates a quadratic function with:
	 * <code>a*x<sup>2</sup> + b*x + c</code>.
	 * 
	 * @param a
	 *            the coefficient for <code>x<sup>2</sup></code>
	 * @param b
	 *            the coefficient for <code>x</code>
	 * @param c
	 *            the constant
	 */
	public QuadraticFunction( double a, double b, double c ) {
		this.a = a;
		this.b = b;
		this.c = c;
	}

	/**
	 * This constructor creates an instance from 3 points <code>(x,f(x))</code>.
	 * 
	 * @param x1
	 *            the first <code>x</code>-coordinate
	 * @param fx1
	 *            <code>f(x1)</code>
	 * @param x2
	 *            the second <code>x</code>-coordinate
	 * @param fx2
	 *            <code>f(x2)</code>
	 * @param x3
	 *            the third <code>x</code>-coordinate
	 * @param fx3
	 *            <code>f(x3)</code>
	 * 
	 * @throws IOException
	 *             if the quadratic function could not be determined
	 */
	public QuadraticFunction( double x1, double fx1, double x2, double fx2, double x3, double fx3 ) throws IOException {
		if( x1 == x2 || x1 == x3 || x2 == x3 ) {
			throw new IOException( "The values of x1, x2 and x3 have to be different." );
		}
		// solve the equation system
		// fx1 = a*x1^2 + b*x1 + c
		// fx2 = a*x2^2 + b*x2 + c
		// fx3 = a*x3^2 + b*x3 + c

		// to save some time we use the variable
		// - a sum/difference has to be used at least 3 times
		// - a product has to be used at least twice

		double x1x1 = x1 * x1, x2x2 = x2 * x2, x2_x1 = x2 - x1;

		this.a = ( ( fx3 - fx1 ) * ( x2_x1 ) - ( fx2 - fx1 ) * ( x3 - x1 ) ) / ( ( x3 * x3 - x1x1 ) * ( x2_x1 ) - ( x2x2 - x1x1 ) * ( x3 - x1 ) );
		this.b = ( fx2 - fx1 - a * ( x2x2 - x1x1 ) ) / ( x2_x1 );
		this.c = fx1 - a * x1x1 - b * x1;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.OneDimensionalFunction#evaluateFunction(double)
	 */
	@Override
	public double evaluateFunction( double x ) {
		return a * x * x + b * x + c;
	}

	/**
	 * This method returns the maximum of the {@link QuadraticFunction}.
	 * 
	 * @return the maximum
	 * 
	 * @throws IOException
	 *             if no maximum exists
	 */
	public double findMax() throws IOException {
		if( a >= 0 ) {
			throw new IOException( "There is no maximum." );
		}
		return getExtremum();
	}

	/**
	 * This method returns the minimum of the {@link QuadraticFunction}.
	 * 
	 * @return the minimum
	 * 
	 * @throws IOException
	 *             if no minimum exists
	 */
	public double findMin() throws IOException {
		if( a <= 0 ) {
			throw new IOException( "There is no minimum." );
		}
		return getExtremum();
	}

	/**
	 * This method returns the extremum of the {@link QuadraticFunction}.
	 * 
	 * @return the extremum
	 * 
	 * @throws IOException
	 *             if no extremum exists
	 */
	public double getExtremum() throws IOException {
		if( a == 0 ) {
			throw new IOException( "There's no extremum." );
		}
		return -b / ( 2 * a );
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String erg = "f(x) = ";
		if( a != 0 ) {
			erg += a + " * x^2 ";
		}
		if( b != 0 ) {
			if( b > 0 ) {
				erg += "+ " + b + " * x ";
			} else {
				erg += b + " * x ";
			}
		}
		if( c != 0 ) {
			if( c > 0 ) {
				erg += "+ " + c;
			} else {
				erg += c;
			}
		}
		if( a == 0 && b == 0 && c == 0 ) {
			erg += 0;
		}
		return erg;
	}
}
