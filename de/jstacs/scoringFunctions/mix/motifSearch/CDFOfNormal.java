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

package de.jstacs.scoringFunctions.mix.motifSearch;

/**
 * This class enables to compute the 
 * 
 * <br><br>
 * 
 * The code is a simplified version of pnorm.c from R (version 2.8.0 downloaded at 07.12.2008 about 19:30).
 * 
 * @author Jens Keilwagen
 */
public class CDFOfNormal {

    private final static double[] a = { 2.2352520354606839287,
										161.02823106855587881,
										1067.6894854603709582,
										18154.981253343561249,
										0.065682337918207449113 };
    
    private final static double[] b = { 47.20258190468824187,
                                        976.09855173777669322,
                                        10260.932208618978205,
                                        45507.789335026729956 };

	private final static double[] c = { 0.39894151208813466764,
										8.8831497943883759412,
										93.506656132177855979,
										597.27027639480026226,
										2494.5375852903726711,
										6848.1904505362823326,
										11602.651437647350124,
										9842.7148383839780218,
										1.0765576773720192317e-8 };

	private final static double[] d = { 22.266688044328115691,
										235.38790178262499861,
										1519.377599407554805,
										6485.558298266760755,
										18615.571640885098091,
										34900.952721145977266,
										38912.003286093271411,
										19685.429676859990727 };

	private final static double[] p = { 0.21589853405795699,
										0.1274011611602473639,
										0.022235277870649807,
										0.001421619193227893466,
										2.9112874951168792e-5,
										0.02307344176494017303 };

	private final static double[] q = { 1.28426009614491121,
										0.468238212480865118,
										0.0659881378689285515,
										0.00378239633202758244,
										7.29751555083966205e-5 };
	
    private static final double M_SQRT_32 = Math.sqrt(32);
    private static final double M_1_SQRT_2PI = 1/Math.sqrt(2*Math.PI);
	
    /**
     * This method computes the logarithm of the cumulative density function of a standard normal distribution.
     * 
     * @param x the value
     * 
     * @return logarithm of the cumulative density function of a standard normal distribution for the value <code>x</code>
     */
	public static double getLogCDF( double x ) {
		double xden, xnum, temp, eps = 1E-5 * 0.5, xsq, xAbs = Math.abs( x ), cum;
		int i;
		boolean negative = x < 0;

		if( xAbs <= 0.67448975 ) {
			/* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
			if( xAbs > eps ) {
				xsq = x * x;
				xnum = a[4] * xsq;
				xden = xsq;
				for( i = 0; i < 3; ++i ) {
					xnum = ( xnum + a[i] ) * xsq;
					xden = ( xden + b[i] ) * xsq;
				}
			} else {
				xnum = xden = 0.0;
			}

			temp = x * ( xnum + a[3] ) / ( xden + b[3] );
			cum = Math.log( 0.5 + temp );
		} else if( xAbs <= M_SQRT_32 ) {
			/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */
			xnum = c[8] * xAbs;
			xden = xAbs;
			for( i = 0; i < 7; ++i ) {
				xnum = ( xnum + c[i] ) * xAbs;
				xden = ( xden + d[i] ) * xAbs;
			}
			temp = ( xnum + c[7] ) / ( xden + d[7] );

			cum = compute( xAbs, temp, negative );
		} else {
			/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
			xsq = 1.0 / ( x * x );
			xnum = p[5] * xsq;
			xden = xsq;
			for( i = 0; i < 4; ++i ) {
				xnum = ( xnum + p[i] ) * xsq;
				xden = ( xden + q[i] ) * xsq;
			}
			temp = xsq * ( xnum + p[4] ) / ( xden + q[4] );
			temp = ( M_1_SQRT_2PI - temp ) / xAbs;

			cum = compute( x, temp, negative );
		}

		return cum;
	}

	private static double compute( double xAbs, double temp, boolean negative ) {
		double xSq = Math.floor( xAbs * 16 ) / 16, t = -0.5 * ( xSq * xSq + ( xAbs - xSq ) * ( xAbs + xSq ) );
		if( negative ) {
			return t + Math.log( temp );
		} else {
			return Math.log1p( -Math.exp( t ) * temp );
		}
	}
}
