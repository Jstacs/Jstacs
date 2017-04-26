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
import java.io.OutputStream;
import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.Time;

/**
 * This class can be used for optimization purposes. The function will be
 * minimized. To maximize a function use the classes {@link NegativeFunction}
 * and {@link NegativeDifferentiableFunction} as input.<br>
 * 
 * All methods for multidimensional optimization need an implementation of
 * {@link DifferentiableFunction}, that will be minimized. After running an
 * multidimensional optimization it is not guaranteed that the optimal parameters
 * are set to the function (due to internal line search, ...). You like to use
 * the optimal parameters in your function you have to set the parameters on
 * your own using the the vector of parameters that has been used during
 * optimization.<br>
 * 
 * The optimizing
 * methods try to do as less as possible function evaluation, which is important
 * for functions that are costly to evaluate. So if you have a costly function
 * and you will fasten the computation, try to do as much effort as you can to
 * make the methods {@link DifferentiableFunction#evaluateFunction(double[])} and
 * {@link DifferentiableFunction#evaluateGradientOfFunction(double[])}
 * in your implementation as fast as possible.<br>
 * 
 * You are enabled to choose different termination modes using different
 * instances of {@link TerminationCondition}.
 * 
 * @author Jens Keilwagen
 * 
 * @see OneDimensionalFunction
 * @see DifferentiableFunction
 * @see NegativeOneDimensionalFunction
 * @see NegativeDifferentiableFunction
 * @see NumericalDifferentiableFunction
 */
public class Optimizer {

	private static final double limFib = ( Math.sqrt( 5 ) - 1 ) / 2d;

	private static final double limFibPlus1 = limFib + 1;

	private static final double limFibPlus2 = limFib + 2;

	/**
	 * This class is used for the limited memory BFGS (<i><b>B</b>royden</i>-
	 * <i>
	 * <b>F</b>letcher</i>-<i><b>G</b>oldfarb</i>-<i><b>S</b>hanno</i>)-method.
	 * 
	 * @author Jens Keilwagen
	 */
	private static class VectorPair {

		private double[] v;

		private double[] s;

		private double alpha;

		private double rho;

		private VectorPair( int n ) {
			v = new double[n];
			s = new double[n];
		}
	}

	/**
	 * This method returns a bracket containing a minimum. This method can be
	 * used to find an initial bracket for a minimization method.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param lower
	 *            the initial value of <code>x</code>, lower bound
	 * @param fLower
	 *            the value <code>f(x)</code>
	 * @param startDistance
	 *            the initial distance for bracketing the minimum
	 * 
	 * @return an array containing <code>x_lower</code>, <code>f(x_lower)</code>
	 *         , <code>x_middle</code>, <code>f(x_middle)</code>,
	 *         <code>x_upper</code> and <code>f(x_upper)</code> at positions 0
	 *         to 5
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * 
	 * @see Optimizer#brentsMethod(OneDimensionalFunction f, double a, double x,
	 *      double fx, double b, double tol )
	 * @see Optimizer#goldenRatio(OneDimensionalFunction f, double lower, double
	 *      p1, double fP1, double upper, double eps )
	 */
	public static double[] findBracket( OneDimensionalFunction f, double lower, double fLower, double startDistance ) throws EvaluationException {
		double[] erg = { lower, fLower, lower + startDistance, 0, 0, 0 };
		erg[3] = f.evaluateFunction( erg[2] );
		//System.out.println("find bracket start" + Arrays.toString(erg) );
		if( Double.isNaN(erg[3]) || erg[1] <= erg[3] ) {
			//shrink upper bound
			//System.out.println("if");
			do {
				erg[4] = erg[2];
				erg[5] = erg[3];
				erg[2] = erg[0] + ( 1 - limFib ) * ( erg[2] - erg[0] );
				erg[3] = f.evaluateFunction( erg[2] );
				//System.out.println(Arrays.toString(erg) );
			} while( Double.isNaN(erg[3]) || (erg[0] != erg[2] && erg[1] < erg[3]) );
		} else {
			//expand upper bound
			//System.out.println("else");
			boolean b;
			do {
				erg[4] = limFibPlus2 * erg[2] - limFibPlus1 * erg[0];
				erg[5] = f.evaluateFunction( erg[4] );
				b = erg[3]>erg[5];
				if( b ) {
					erg[0] = erg[2];
					erg[1] = erg[3];
					erg[2] = erg[4];
					erg[3] = erg[5];
				}
				//System.out.println(Arrays.toString(erg) );
			} while( !Double.isNaN(erg[5]) && b );
		}
		if( Double.isNaN(erg[5]) ) {
			//nicht optimal
			//System.out.println("PROBLEM");
			//System.out.println(Arrays.toString(erg) );
			do{
				double x = erg[2] + ( 1 - limFib ) * ( erg[4] - erg[2] );
				if( Double.compare( erg[2], x) == 0 ) {
					break;
				}
				double fx = f.evaluateFunction( x );
				if( Double.isNaN(fx) || fx > erg[3] ) {
					erg[4]=x;
					erg[5]=fx;
				} else {					
					if( fx != erg[3] ) {
						erg[0] = erg[2];
						erg[1] = erg[3];
					}
					erg[2] = x;
					erg[3] = fx;
				}
				//System.out.println(Arrays.toString(erg) );
			} while( Double.isNaN(erg[5]) );
			
			if( Double.isNaN(erg[3]) ) {
				erg[5] = erg[3];
				erg[4] = erg[2];
			}
		}
		//System.out.println("find bracket end" + Arrays.toString(erg) );
		return erg;
	}

	/**
	 * This method returns a bracket containing a minimum.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param lower
	 *            the initial value of <code>x</code>, lower bound
	 * @param startDistance
	 *            the initial distance for bracketing the minimum
	 * 
	 * @return an array containing <code>x_lower</code>, <code>f(x_lower)</code>
	 *         , <code>x_inside</code>, <code>f(x_inside)</code>,
	 *         <code>x_upper</code> and <code>f(x_upper)</code> at positions 0
	 *         to 5
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * 
	 * @see Optimizer#findBracket(OneDimensionalFunction, double, double,
	 *      double)
	 */
	public static double[] findBracket( OneDimensionalFunction f, double lower, double startDistance ) throws EvaluationException {
		return findBracket( f, lower, f.evaluateFunction( lower ), startDistance );
	}

	/**
	 * Approximates a minimum (not necessary the global) in the interval
	 * <code>[lower,upper]</code>. The method is using the &quot;Brent's
	 * method&quot;. <br>
	 * <br>
	 * 
	 * This method is a modified copy of a method downloaded from the internet.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param a
	 *            the initial lower boundary
	 * @param x
	 *            an abscissa inside the interval <code>[a,b]</code>
	 * @param fx
	 *            the ordinate corresponding to <code>x</code>
	 * @param b
	 *            the initial upper bound
	 * @param tol
	 *            the termination parameter (interval size)
	 * 
	 * @return an array containing the minimum <code>x</code> (position 0) and
	 *         the value of the minimum <code>f(x)</code> (position 1)
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 */
	public static double[] brentsMethod( OneDimensionalFunction f, double a, double x, double fx, double b, double tol ) throws EvaluationException {
		int i = 0;
		double c = 1 - limFib, d = 0, e, eps, xm, p, q, r, tol1, t2, u, v, w, fu, fv, fw, tol3;

		eps = 1.2e-16;
		tol1 = eps + 1.0;
		eps = Math.sqrt( eps );
		e = 0.0;
		v = w = x;
		fv = fw = fx;

		tol3 = tol / 3.0;

		xm = .5 * ( a + b );
		tol1 = eps * Math.abs( x ) + tol3;
		t2 = 2.0 * tol1;

		// main loop
		while( Math.abs( x - xm ) > ( t2 - .5 * ( b - a ) ) ) {
			p = q = r = 0.0;
			if( Math.abs( e ) > tol1 ) {
				// fit the parabola
				r = ( x - w ) * ( fx - fv );
				q = ( x - v ) * ( fx - fw );
				p = ( x - v ) * q - ( x - w ) * r;
				q = 2.0 * ( q - r );
				if( q > 0.0 ) {
					p = -p;
				} else {
					q = -q;
				}
				r = e;
				e = d;
			}
			if( ( Math.abs( p ) < Math.abs( .5 * q * r ) ) && ( p > q * ( a - x ) ) && ( p < q * ( b - x ) ) ) {
				// a parabolic interpolation step
				d = p / q;
				u = x + d;
				// f must not be evaluated too close to a or b
				if( ( ( u - a ) < t2 ) || ( ( b - u ) < t2 ) ) {
					d = tol1;
					if( x >= xm ) d = -d;
				}
			} else {
				// a golden-section step
				if( x < xm ) {
					e = b - x;
				} else {
					e = a - x;
				}
				d = c * e;
			}
			// f must not be evaluated too close to x
			if( Math.abs( d ) >= tol1 ) {
				u = x + d;
			} else {
				if( d > 0.0 ) {
					u = x + tol1;
				} else {
					u = x - tol1;
				}
			}
			fu = f.evaluateFunction( u );
			// Update a, b, v, w, and x
			if( fx <= fu ) {
				if( u < x ) {
					a = u;
				} else {
					b = u;
				}
			}
			if( fu <= fx ) {
				if( u < x ) {
					b = x;
				} else {
					a = x;
				}
				v = w;
				fv = fw;
				w = x;
				fw = fx;
				x = u;
				fx = fu;
			} else {
				if( ( fu <= fw ) || ( w == x ) ) {
					v = w;
					fv = fw;
					w = u;
					fw = fu;
				} else if( !( ( fu > fv ) && ( v != x ) && ( v != w ) ) ) {
					v = u;
					fv = fu;
				}
			}
			xm = .5 * ( a + b );
			tol1 = eps * Math.abs( x ) + tol3;
			t2 = 2.0 * tol1;
			i++;
		}
		double[] erg = { x, fx };
		return erg;
	}

	/**
	 * Approximates a minimum (not necessary the global) in the interval
	 * <code>[lower,upper]</code>. The method is using the &quot;Brent's
	 * method&quot;. <br>
	 * <br>
	 * 
	 * This method is an modified copy of a method downloaded from the internet.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param a
	 *            the initial lower boundary
	 * @param x
	 *            an abscissa inside the interval <code>[a,b]</code>
	 * @param b
	 *            the initial upper bound
	 * @param tol
	 *            the termination parameter (interval size)
	 * 
	 * @return an array containing the minimum <code>x</code> (position 0) and
	 *         the value of the minimum <code>f(x)</code> (position 1)
	 * 
	 * @throws EvaluationException
	 *             if there was a mistake in evaluating the function
	 * 
	 * @see Optimizer#brentsMethod(OneDimensionalFunction, double, double,
	 *      double, double, double)
	 */
	public static double[] brentsMethod( OneDimensionalFunction f, double a, double x, double b, double tol ) throws EvaluationException {
		return brentsMethod( f, a, x, f.evaluateFunction( x ), b, tol );
	}

	/**
	 * Approximates a minimum (not necessary the global) in the interval
	 * <code>[lower,upper]</code>. The method is using the &quot;golden
	 * section&quot;.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param lower
	 *            the lower interval boundary
	 * @param p1
	 *            an abscissa in the interval <code>[lower,upper]</code>
	 * @param fP1
	 *            the ordinate corresponding to <code>p1</code>
	 * @param upper
	 *            the upper interval boundary
	 * @param eps
	 *            the minimal interval size
	 * 
	 * @return an array containing the minimum <code>x</code> (position 0) and
	 *         the value of the minimum <code>f(x)</code> (position 1)
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 */
	public static double[] goldenRatio( OneDimensionalFunction f, double lower, double p1, double fP1, double upper, double eps ) throws EvaluationException {
		double[] erg = { lower + limFib * ( upper - lower ), 0 };
		boolean fP1Computed = true;

		int i = 0;
		while( upper - lower > eps ) {
			i++;
			if( !fP1Computed ) {
				fP1 = f.evaluateFunction( p1 );
			} else {
				erg[1] = f.evaluateFunction( erg[0] );
			}
			if( fP1 < erg[1] ) {
				// lower and fLower don't change
				upper = erg[0];
				erg[0] = p1;
				erg[1] = fP1;
				p1 = lower + ( 1 - limFib ) * ( upper - lower );
				fP1Computed = false;
			} else {
				// upper and fUpper don't change
				lower = p1;
				p1 = erg[0];
				fP1 = erg[1];
				erg[0] = lower + limFib * ( upper - lower );
				fP1Computed = true;
			}
		}
		if( fP1Computed ) {
			erg[0] = p1;
			erg[1] = fP1;
		}
		return erg;
	}

	/**
	 * Approximates a minimum (not necessary the global) in the interval
	 * <code>[lower,upper]</code>. The method is using the &quot;golden
	 * section&quot;.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param lower
	 *            the lower interval boundary
	 * @param upper
	 *            the upper interval boundary
	 * @param eps
	 *            the minimal interval size
	 * 
	 * @return an array containing the minimum <code>x</code> (position 0) and
	 *         the value of the minimum <code>f(x)</code> (position 1)
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * 
	 * @see Optimizer#goldenRatio(OneDimensionalFunction, double, double,
	 *      double, double, double)
	 */
	public static double[] goldenRatio( OneDimensionalFunction f, double lower, double upper, double eps ) throws EvaluationException {
		return goldenRatio( f, lower, lower + ( 1 - limFib ) * ( upper - lower ), f.evaluateFunction( lower + ( 1 - limFib )
																										* ( upper - lower ) ), upper, eps );
	}

	/**
	 * Approximates a minimum (not necessary the global) in the interval
	 * <code>[lower,upper]</code>. The method is using the &quot;parabolic
	 * interpolation&quot;.<br>
	 * <br>
	 * 
	 * This method should work well for convex functions.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param lower
	 *            the lower interval boundary
	 * @param fLower
	 *            the ordinate corresponding to lower
	 * @param p1
	 *            an abscissa in the interval <code>[lower,upper]</code>
	 * @param fP1
	 *            the ordinate corresponding to <code>p1</code>
	 * @param upper
	 *            the upper interval boundary
	 * @param fUpper
	 *            the ordinate corresponding to upper
	 * @param eps
	 *            the minimal interval size
	 * 
	 * @return an array containing the minimum <code>x</code> (position 0) and
	 *         the value of the minimum <code>f(x)</code> (position 1)
	 * 
	 * @throws Exception
	 *             if there was something wrong during the evaluation of the
	 *             function
	 */
	public static double[] parabolicInterpolation( OneDimensionalFunction f, double lower, double fLower, double p1, double fP1,
			double upper, double fUpper, double eps ) throws Exception {
		double[] erg = { p1, fP1 };
		double eps2 = 2 * eps;
		int i = 0;
		while( upper - lower > eps ) {
			if( fLower == fP1 && fP1 == fUpper ) {
				// the parabolic interpolation will not find any minimum, since
				// the points are colinear
				return erg;
			}
			i++;
			erg[0] = ( new QuadraticFunction( lower, fLower, p1, fP1, upper, fUpper ) ).findMin();
			// the "if" is to ensure, that the function is not evaluated to
			// close
			if( Math.abs( erg[0] - p1 ) <= eps ) {
				// erg[0] is inside [lower,p1] and close to p1
				if( erg[0] < p1 ) {
					if( p1 - lower <= eps2 ) {
						erg[0] = lower + ( p1 - lower ) / 2d;
					} else {
						erg[0] = p1 - eps;
					}
				}
				// erg[0] is inside [p1, upper] and close to p1
				else {
					if( upper - p1 <= eps2 ) {
						erg[0] = p1 + ( p1 - lower ) / 2d;
					} else {
						erg[0] = p1 + eps;
					}
				}
			}
			erg[1] = f.evaluateFunction( erg[0] );
			if( erg[0] < p1 ) {
				if( erg[1] < fP1 ) {
					// lower and fLower don't change
					upper = p1;
					fUpper = fP1;
					p1 = erg[0];
					fP1 = erg[1];
				} else {
					// upper and fUpper don't change
					lower = erg[0];
					fLower = erg[1];
				}
			} else {
				if( erg[1] < fP1 ) {
					// upper and fUpper don't change
					lower = p1;
					fLower = fP1;
					p1 = erg[0];
					fP1 = erg[1];
				} else {
					// lower and fLower don't change
					upper = erg[0];
					fUpper = erg[1];
				}
			}
		}
		return erg;
	}

	/**
	 * Approximates a minimum (not necessary the global) in the interval
	 * <code>[lower,upper]</code>. The method is using the &quot;parabolic
	 * interpolation&quot;.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param lower
	 *            the lower interval boundary
	 * @param p1
	 *            an abscissa in the interval <code>[lower,upper]</code>
	 * @param upper
	 *            the upper interval boundary
	 * @param eps
	 *            the minimal interval size
	 * 
	 * @return an array containing the minimum <code>x</code> (position 0) and
	 *         the value of minimum <code>f(x)</code> (position 1)
	 * 
	 * @throws Exception
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * 
	 * @see Optimizer#parabolicInterpolation(OneDimensionalFunction, double,
	 *      double, double, double, double, double, double)
	 */
	public static double[] parabolicInterpolation( OneDimensionalFunction f, double lower, double p1, double upper, double eps ) throws Exception {
		return parabolicInterpolation( f,
				lower,
				f.evaluateFunction( lower ),
				p1,
				f.evaluateFunction( p1 ),
				upper,
				f.evaluateFunction( upper ),
				eps );
	}

	private static void checkStartDistance( double startDistance ) throws IllegalArgumentException {
		if( startDistance <= 0 ) {
			throw new IllegalArgumentException( "The startDistance has to be greater than 0." );
		}
	}

	/**
	 * The steepest descent.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param currentValues
	 *            at the begin the start vector and at the end the minimum
	 *            vector
	 * @param terminationMode
	 *            the choice how to terminate the algorithm
	 * @param linEps
	 *            the bound for stopping the linear search
	 * @param startDistance
	 *            the initial step for the linear search
	 * @param out
	 *            the {@link OutputStream} for writing some information
	 *            about the iterations, if <code>null</code> nothing will be
	 *            written
	 * @param t
	 *            a {@link Time} instance used to measure the elapsed time
	 * 
	 * @return the number of iterations
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * @throws TerminationException
	 *             if the termination condition is unknown
	 * @throws DimensionException
	 *             if the dimension of the scope of the function and the array
	 *             of arguments do not match
	 * @throws IOException
	 *             if there was something wrong with the IO
	 */
	public static int steepestDescent( DifferentiableFunction f, double[] currentValues, TerminationCondition terminationMode,
			double linEps, StartDistanceForecaster startDistance, OutputStream out, Time t ) throws DimensionException,
			TerminationException,
			EvaluationException,
			IOException {
		SafeOutputStream myOut = SafeOutputStream.getSafeOutputStream( out );
		int i = 0, n = f.getDimensionOfScope();
		double[] d = new double[n];
		Arrays.fill( d, Double.MAX_VALUE );
		double current = Double.POSITIVE_INFINITY;
		double[] help = { Double.POSITIVE_INFINITY, f.evaluateFunction( currentValues ) }, gradient = f.evaluateGradientOfFunction( currentValues );
		
		int counter;
		double sd;
		myOut.writeln( "iteration\ttime\tf(x)\tdelta\tstart distance\tlinesearch" );
		myOut.writeln( "0\t0\t" + help[1] + "\t0" );
		while( terminationMode.doNextIteration( i, current, help[1], gradient, d, help[0], t ) ) {
			current = help[1];
			d = gradient;
			for( counter = 0; counter < n; counter++ ) {
				d[counter] *= -1d;
			}
			sd = startDistance.getNewStartDistance();
			checkStartDistance( sd );
			help = f.findOneDimensionalMin( currentValues, d, 0, current, linEps, sd );
			myOut.writeln( ( ++i ) + " \t"
							+ t.getElapsedTime()
							+ " \t"
							+ help[1]
							+ " \t"
							+ ( current - help[1] )
							+ "\t"
							+ sd
							+ "\t"
							+ help[0] );
			startDistance.setLastDistance( help[0] );
			for( counter = 0; counter < n; counter++ ) {
				currentValues[counter] += d[counter] * help[0];
			}
			gradient = f.evaluateGradientOfFunction( currentValues );
		}
		return i;
	}

	/**
	 * The conjugate gradient algorithm by <i>Fletcher</i> and <i>Reeves</i>.
	 * 
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param currentValues
	 *            at the begin the start vector and at the end the minimum
	 *            vector
	 * @param terminationMode
	 *            the choice how to terminate the algorithm
	 * @param linEps
	 *            the bound for stopping the linear search
	 * @param startDistance
	 *            the initial step for the linear search
	 * @param out
	 *            the {@link OutputStream} for writing some information
	 *            about the iterations, if <code>null</code> nothing will be
	 *            written
	 * @param t
	 *            a {@link Time} instance used to measure the elapsed time
	 * 
	 * @return the number of iterations
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * @throws TerminationException
	 *             if the termination condition is unknown
	 * @throws DimensionException
	 *             if the dimension of the scope of the function and the array
	 *             of arguments do not match
	 * @throws IOException
	 *             if there was something wrong with the IO
	 */
	public static int conjugateGradientsFR( DifferentiableFunction f, double[] currentValues, TerminationCondition terminationMode,
			double linEps, StartDistanceForecaster startDistance, OutputStream out, Time t ) throws DimensionException,
			TerminationException,
			IOException,
			EvaluationException {
		SafeOutputStream myOut = SafeOutputStream.getSafeOutputStream( out );
		int i = 0, n = f.getDimensionOfScope();
		double[] d = new double[n];
		Arrays.fill( d, Double.MAX_VALUE );
		double current = Double.POSITIVE_INFINITY;
		double[] help = { Double.POSITIVE_INFINITY, f.evaluateFunction( currentValues ) }, gradient = f.evaluateGradientOfFunction( currentValues );
		
		int counter;
		double mu1 = 0, mu2 = 1, sd;
		for( counter = 0; counter < n; counter++ ) {
			mu1 += gradient[counter] * gradient[counter];
		}
		boolean next = terminationMode.doNextIteration( i, current, help[1], gradient, d, help[0], t );
		Arrays.fill( d, 0 );
		
		myOut.writeln( "iteration\ttime\tf(x)\tdelta\tstart distance\tlinesearch" );
		myOut.writeln( "0\t0\t" + help[1] + "\t0" );
		while( next ) {
			current = help[1];
			for( counter = 0; counter < n; counter++ ) {
				d[counter] = -gradient[counter] + mu1 / mu2 * d[counter];
			}
			sd = startDistance.getNewStartDistance();
			checkStartDistance( sd );
			help = f.findOneDimensionalMin( currentValues, d, 0, current, linEps, sd );
			myOut.writeln( ( ++i ) + "\t" + t.getElapsedTime() + "\t" + help[1] + "\t" + ( current - help[1] ) + "\t" + sd + "\t" + help[0] );
			startDistance.setLastDistance( help[0] );
			for( counter = 0; counter < n; counter++ ) {
				currentValues[counter] += d[counter] * help[0];
			}
			mu2 = mu1;
			mu1 = 0;
			gradient = f.evaluateGradientOfFunction( currentValues );
			for( counter = 0; counter < n; counter++ ) {
				mu1 += gradient[counter] * gradient[counter];
			}
			next = terminationMode.doNextIteration( i, current, help[1], gradient, d, help[0], t );
		}
		return i;
	}

	/**
	 * The conjugate gradient algorithm by <i>Polak</i> and <i>Ribi&egrave;re</i>.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param currentValues
	 *            at the begin the start vector and at the end the minimum
	 *            vector
	 * @param terminationMode
	 *            the choice how to terminate the algorithm
	 * @param linEps
	 *            the bound for stopping the linear search
	 * @param startDistance
	 *            the initial step for the linear search
	 * @param out
	 *            the {@link OutputStream} for writing some information
	 *            about the iterations, if <code>null</code> nothing will be
	 *            written
	 * @param t
	 *            a {@link Time} instance used to measure the elapsed time
	 * 
	 * @return the number of iterations
	 * 
	 * @throws EvaluationException
	 *             if there something wrong during the evaluation of the
	 *             function
	 * @throws TerminationException
	 *             if the termination condition is unknown
	 * @throws DimensionException
	 *             if the dimension of the scope of the function and the array
	 *             of arguments do not match
	 * @throws IOException
	 *             if there was something wrong with the IO
	 * 
	 * @see Optimizer#conjugateGradientsPRP(DifferentiableFunction, double[], TerminationCondition, double, StartDistanceForecaster, OutputStream, boolean, Time)
	 */
	public static int conjugateGradientsPR( DifferentiableFunction f, double[] currentValues, TerminationCondition terminationMode,
			double linEps, StartDistanceForecaster startDistance, OutputStream out, Time t ) throws DimensionException,
			TerminationException,
			IOException,
			EvaluationException {
		return conjugateGradientsPRP( f, currentValues, terminationMode, linEps, startDistance, out, false, t );
	}

	/**
	 * The conjugate gradient algorithm by <i>Polak</i> and <i>Ribi&egrave;re</i>
	 * called &quot;Polak-Ribi&egrave;re-Positive&quot;.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param currentValues
	 *            at the begin the start vector and at the end the minimum
	 *            vector
	 * @param terminationMode
	 *            the choice how to terminate the algorithm
	 * @param linEps
	 *            the bound for stopping the linear search
	 * @param startDistance
	 *            the initial step for the linear search
	 * @param out
	 *            the {@link OutputStream} for writing some information
	 *            about the iterations, if <code>null</code> nothing will be
	 *            written.
	 * @param t
	 *            a {@link Time} instance used to measure the elapsed time
	 * 
	 * @return the number of iterations
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong with the evaluation of the
	 *             function
	 * @throws TerminationException
	 *             if the termination condition is unknown
	 * @throws DimensionException
	 *             if the dimension of the scope of the function and the array
	 *             of arguments do not match
	 * @throws IOException
	 *             if there was something wrong with the IO
	 * 
	 * @see Optimizer#conjugateGradientsPRP(DifferentiableFunction, double[], TerminationCondition, double, StartDistanceForecaster, OutputStream, boolean, Time)
	 */
	public static int conjugateGradientsPRP( DifferentiableFunction f, double[] currentValues, TerminationCondition terminationMode,
			double linEps, StartDistanceForecaster startDistance, OutputStream out, Time t ) throws DimensionException,
			TerminationException,
			IOException,
			EvaluationException {
		return conjugateGradientsPRP( f, currentValues, terminationMode, linEps, startDistance, out, true, t );
	}

	private static int conjugateGradientsPRP( DifferentiableFunction f, double[] currentValues, TerminationCondition terminationMode,
			double linEps, StartDistanceForecaster startDistance, OutputStream out, boolean positive, Time t ) throws DimensionException,
			TerminationException,
			IOException,
			EvaluationException,
			IllegalArgumentException {
		SafeOutputStream myOut = SafeOutputStream.getSafeOutputStream( out );
		int i = 0, n = f.getDimensionOfScope();
		double[] d = new double[n];
		Arrays.fill( d, Double.MAX_VALUE );
		double current = Double.POSITIVE_INFINITY;
		double[] help = { Double.POSITIVE_INFINITY, f.evaluateFunction( currentValues ) };

		int counter;
		double[] gradient_old, gradient_new = f.evaluateGradientOfFunction( currentValues );
		gradient_old = gradient_new;
		double mu1, mu2, sd;

		myOut.writeln( "iteration\ttime\tf(x)\tdelta\tstart distance\tlinesearch" );
		myOut.writeln( "0\t0\t" + help[1] + "\t0" );
		while( terminationMode.doNextIteration( i, current, help[1], gradient_new, d, help[0], t ) ) {
			current = help[1];
			mu2 = 0;
			mu1 = 0;
			for( counter = 0; counter < n; counter++ ) {
				mu1 += ( gradient_new[counter] - gradient_old[counter] ) * gradient_new[counter];
				mu2 += gradient_old[counter] * gradient_old[counter];
			}
			if( positive ) {
				mu1 = Math.max( 0, mu1 / mu2 );
			} else {
				mu1 /= mu2;
			}
			for( counter = 0; counter < n; counter++ ) {
				d[counter] = -gradient_new[counter] + mu1 * d[counter];
			}
			sd = startDistance.getNewStartDistance();
			checkStartDistance( sd );
			help = f.findOneDimensionalMin( currentValues, d, 0, current, linEps, sd );
			myOut.writeln( ( ++i ) + "\t" + t.getElapsedTime() + "\t" + help[1] + "\t" + ( current - help[1] ) + "\t" + sd + "\t" + help[0] );
			startDistance.setLastDistance( help[0] );
			for( counter = 0; counter < n; counter++ ) {
				currentValues[counter] += d[counter] * help[0];
			}
			gradient_old = gradient_new;
			gradient_new = f.evaluateGradientOfFunction( currentValues );
		}
		return i;
	}

	/**
	 * The <i>Davidon</i>-<i>Fletcher</i>-<i>Powell</i> version of the
	 * quasi-Newton method.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param currentValues
	 *            at the begin the start vector and at the end the minimum
	 *            vector
	 * @param terminationMode
	 *            the choice how to terminate the algorithm
	 * @param linEps
	 *            the bound for stopping the linear search
	 * @param startDistance
	 *            the initial step for the linear search
	 * @param out
	 *            the {@link OutputStream} for writing some information
	 *            about the iterations, if <code>null</code> nothing will be
	 *            written.
	 * @param t
	 *            a {@link Time} instance used to measure the elapsed time
	 * 
	 * @return the number of iterations
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * @throws TerminationException
	 *             if the termination condition is unknown
	 * @throws DimensionException
	 *             if the dimension of the scope of the function and the array
	 *             of arguments do not match
	 * @throws IOException
	 *             if there was something wrong with the IO
	 */
	public static int quasiNewtonDFP( DifferentiableFunction f, double[] currentValues, TerminationCondition terminationMode,
			double linEps, StartDistanceForecaster startDistance, OutputStream out, Time t ) throws DimensionException,
			TerminationException,
			IOException,
			EvaluationException {
		SafeOutputStream myOut = SafeOutputStream.getSafeOutputStream( out );
		int i = 0, n = f.getDimensionOfScope();
		double[] d = new double[n];
		double current = Double.POSITIVE_INFINITY;
		double[] help = { Double.POSITIVE_INFINITY, f.evaluateFunction( currentValues ) };
		
		int counter1, counter2;
		double[] s = new double[n], v = new double[n], matrixv = new double[n];
		double[] gradient_old, gradient_new = f.evaluateGradientOfFunction( currentValues );
		double[][] matrix = new double[n][n];
		double vmatrixv, sv, vv, sd;
		for( counter1 = 0; counter1 < n; counter1++ ) {
			d[counter1] = -gradient_new[counter1];
			matrix[counter1][counter1] = 1;
		}
		boolean next = terminationMode.doNextIteration( i, current, help[1], gradient_new, d, help[0], t );

		myOut.writeln( "iteration\ttime\tf(x)\tdelta\tstart distance\tlinesearch" );
		myOut.writeln( "0\t0\t" + help[1] + "\t0" );
		while( next ) {
			current = help[1];
			sd = startDistance.getNewStartDistance();
			checkStartDistance( sd );
			help = f.findOneDimensionalMin( currentValues, d, 0, current, linEps, sd );
			myOut.writeln( ( ++i ) + "\t" + t.getElapsedTime() + "\t" + help[1] + "\t" + ( current - help[1] ) + "\t" + sd + "\t" + help[0] );
			startDistance.setLastDistance( help[0] );
			for( counter1 = 0; counter1 < n; counter1++ ) {
				currentValues[counter1] += d[counter1] * help[0];
			}
			gradient_old = gradient_new;
			gradient_new = f.evaluateGradientOfFunction( currentValues );
			next = terminationMode.doNextIteration( i, current, help[1], gradient_new, d, help[0], t );
			if( next ) {
				for( counter1 = 0; counter1 < n; counter1++ ) {
					s[counter1] = d[counter1] * help[0];
					v[counter1] = gradient_new[counter1] - gradient_old[counter1];
				}

				vmatrixv = sv = vv = 0;

				for( counter1 = 0; counter1 < n; counter1++ ) {
					sv += s[counter1] * v[counter1];
					vv += v[counter1] * v[counter1];
					matrixv[counter1] = 0;
					for( counter2 = 0; counter2 < n; counter2++ ) {
						matrixv[counter1] += matrix[counter1][counter2] * v[counter2];
					}
					vmatrixv += v[counter1] * matrixv[counter1];
				}
				for( counter1 = 0; counter1 < n; counter1++ ) {
					for( counter2 = 0; counter2 < n; counter2++ ) {
						matrix[counter1][counter2] += s[counter1] * s[counter2] / sv - matrixv[counter1] * matrixv[counter2] / vmatrixv;
					}
					d[counter1] = 0;
					for( counter2 = 0; counter2 < n; counter2++ ) {
						d[counter1] += -matrix[counter1][counter2] * gradient_new[counter2];
					}
				}
			}
		}
		return i;
	}

	/**
	 * The <i>Broyden</i>-<i>Fletcher</i>-<i>Goldfarb</i>-<i>Shanno</i> version
	 * of the quasi-Newton method.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param currentValues
	 *            at the begin the start vector and at the end the minimum
	 *            vector
	 * @param terminationMode
	 *            the choice how to terminate the algorithm
	 * @param linEps
	 *            the bound for stopping the linear search
	 * @param startDistance
	 *            the initial step for the linear search
	 * @param out
	 *            the {@link OutputStream} for writing some information
	 *            about the iterations, if <code>null</code> nothing will be
	 *            written.
	 * @param t
	 *            a {@link Time} instance used to measure the elapsed time
	 * 
	 * @return the number of iterations
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * @throws TerminationException
	 *             if the termination condition is unknown
	 * @throws DimensionException
	 *             if the dimension of the scope of the function and the array
	 *             of arguments do not match
	 * @throws IOException
	 *             if there was something wrong with the IO
	 */
	public static int quasiNewtonBFGS( DifferentiableFunction f, double[] currentValues, TerminationCondition terminationMode,
			double linEps, StartDistanceForecaster startDistance, OutputStream out, Time t ) throws DimensionException,
			TerminationException,
			IOException,
			EvaluationException {
		SafeOutputStream myOut = SafeOutputStream.getSafeOutputStream( out );
		int i = 0, n = f.getDimensionOfScope();
		double[] d = new double[n];
		double current = Double.POSITIVE_INFINITY;
		double[] help = { Double.POSITIVE_INFINITY, f.evaluateFunction( currentValues ) };
		
		int counter1, counter2;
		double[] s = new double[n], v = new double[n], u = new double[n], matrixv = new double[n];
		double[] gradient_old, gradient_new = f.evaluateGradientOfFunction( currentValues );
		double[][] matrix = new double[n][n];
		double vmatrixv, sv, vv, sd;
		for( counter1 = 0; counter1 < n; counter1++ ) {
			d[counter1] = -gradient_new[counter1];
			matrix[counter1][counter1] = 1;
		}
		boolean next = terminationMode.doNextIteration( i, current, help[1], gradient_new, d, help[0], t );
		
		myOut.writeln( "iteration\ttime\tf(x)\tdelta\tstart distance\tlinesearch" );
		myOut.writeln( "0\t0\t" + help[1] + "\t0" );
		while( next ) {
			current = help[1];
			sd = startDistance.getNewStartDistance();
			checkStartDistance( sd );
			help = f.findOneDimensionalMin( currentValues, d, 0, current, linEps, sd );
			myOut.writeln( ( ++i ) + "\t" + t.getElapsedTime() + "\t" + help[1] + "\t" + ( current - help[1] ) + "\t" + sd + "\t" + help[0] );
			startDistance.setLastDistance( help[0] );
			for( counter1 = 0; counter1 < n; counter1++ ) {
				currentValues[counter1] += d[counter1] * help[0];
			}
			gradient_old = gradient_new;
			gradient_new = f.evaluateGradientOfFunction( currentValues );
			next = terminationMode.doNextIteration( i, current, help[1], gradient_new, d, help[0], t ); 
			if( next ) {
				for( counter1 = 0; counter1 < n; counter1++ ) {
					s[counter1] = d[counter1] * help[0];
					v[counter1] = gradient_new[counter1] - gradient_old[counter1];
				}

				vmatrixv = sv = vv = 0;

				for( counter1 = 0; counter1 < n; counter1++ ) {
					sv += s[counter1] * v[counter1];
					vv += v[counter1] * v[counter1];
					matrixv[counter1] = 0;
					for( counter2 = 0; counter2 < n; counter2++ ) {
						matrixv[counter1] += matrix[counter1][counter2] * v[counter2];
					}
					vmatrixv += v[counter1] * matrixv[counter1];
				}
				for( counter1 = 0; counter1 < n; counter1++ ) {
					u[counter1] = s[counter1] / sv - matrixv[counter1] / vmatrixv;
				}
				for( counter1 = 0; counter1 < n; counter1++ ) {
					for( counter2 = 0; counter2 < n; counter2++ ) {
						matrix[counter1][counter2] += s[counter1] * s[counter2]
														/ sv
														- matrixv[counter1]
														* matrixv[counter2]
														/ vmatrixv
														+ vmatrixv
														* u[counter1]
														* u[counter2];
					}
					d[counter1] = 0;
					for( counter2 = 0; counter2 < n; counter2++ ) {
						d[counter1] += -matrix[counter1][counter2] * gradient_new[counter2];
					}
				}
			}
		}
		return i;
	}

	/**
	 * The <i>Broyden</i>-<i>Fletcher</i>-<i>Goldfarb</i>-<i>Shanno</i> version
	 * of limited memory quasi-Newton methods.
	 * 
	 * @param f
	 *            the function to be minimized
	 * @param currentValues
	 *            at the begin the start vector and at the end the minimum
	 *            vector
	 * @param m
	 *            the current number of vectors used to approximate the Hessian,
	 *            typically between 3 and 10
	 * @param terminationMode
	 *            the choice how to terminate the algorithm
	 * @param linEps
	 *            the bound for stopping the linear search
	 * @param startDistance
	 *            the initial step for the linear search
	 * @param out
	 *            the {@link OutputStream} for writing some information
	 *            about the iterations, if <code>null</code> nothing will be
	 *            written.
	 * @param t
	 *            a {@link Time} instance used to measure the elapsed time
	 * 
	 * @return the number of iterations
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * @throws TerminationException
	 *             if the termination condition is unknown
	 * @throws DimensionException
	 *             if the dimension of the scope of the function and the array
	 *             of arguments do not match
	 * @throws IOException
	 *             if there was something wrong with the IO
	 */
	public static int limitedMemoryBFGS( DifferentiableFunction f, double[] currentValues, byte m, TerminationCondition terminationMode,
			double linEps, StartDistanceForecaster startDistance, OutputStream out, Time t ) throws DimensionException,
			TerminationException,
			IOException,
			EvaluationException {
		if( m <= 2 ) {
			throw new IllegalArgumentException( "This choice of m is not allowed." );
		}
		SafeOutputStream myOut = SafeOutputStream.getSafeOutputStream( out );
		int i = 0, n = f.getDimensionOfScope();
		double[] d = new double[n];
		double current = Double.POSITIVE_INFINITY;
		double[] help = { Double.POSITIVE_INFINITY, f.evaluateFunction( currentValues ) };
		
		int s, counter1, counter2;
		LinkedList<VectorPair> vp = new LinkedList<VectorPair>();
		VectorPair now, newVp;
		double[] gradient_old, gradient_new = f.evaluateGradientOfFunction( currentValues );
		double beta, yy, sd;
		for( counter1 = 0; counter1 < n; counter1++ ) {
			d[counter1] = -gradient_new[counter1];
		}
		boolean next = terminationMode.doNextIteration( i, current, help[1], gradient_new, d, help[0], t );
		
		myOut.writeln( "iteration\ttime\tf(x)\tdelta\tstart distance\tlinesearch" );
		myOut.writeln( "0\t0\t" + help[1] + "\t0" );
		while( next ) {
			current = help[1];
			sd = startDistance.getNewStartDistance();
			checkStartDistance( sd );
			help = f.findOneDimensionalMin( currentValues, d, 0, current, linEps, sd );
			myOut.writeln( ( ++i ) + "\t" + t.getElapsedTime() + "\t" + help[1] + "\t" + ( current - help[1] ) + "\t" + sd + "\t" + help[0] );
			startDistance.setLastDistance( help[0] );
			newVp = new VectorPair( n );
			newVp.rho = 0;
			for( counter1 = 0; counter1 < n; counter1++ ) {
				newVp.s[counter1] = d[counter1] * help[0];
				currentValues[counter1] += newVp.s[counter1];
			}
			gradient_old = gradient_new;
			gradient_new = f.evaluateGradientOfFunction( currentValues );
			next = terminationMode.doNextIteration( i, current, help[1], gradient_new, d, help[0], t );
			if( next ) {
				for( counter1 = 0; counter1 < n; counter1++ ) {
					newVp.v[counter1] = gradient_new[counter1] - gradient_old[counter1];
					d[counter1] = -gradient_new[counter1];
					newVp.rho += newVp.s[counter1] * newVp.v[counter1];
				}
				newVp.rho = 1d / newVp.rho;

				vp.addFirst( newVp );
				s = vp.size();
				if( s > m ) {
					vp.removeLast();
					s--;
				}
				for( counter1 = 0; counter1 < s; counter1++ ) {
					now = (VectorPair)vp.get( counter1 );
					now.alpha = 0;
					// compute alpha
					for( counter2 = 0; counter2 < n; counter2++ ) {
						now.alpha += now.s[counter2] * d[counter2];
					}
					now.alpha *= now.rho;
					// compute q (=d)
					for( counter2 = 0; counter2 < n; counter2++ ) {
						d[counter2] -= now.alpha * now.v[counter2];
					}
				}

				beta = 0;
				yy = 0;
				for( counter2 = 0; counter2 < n; counter2++ ) {
					beta += newVp.v[counter2] * newVp.s[counter2];
					yy += newVp.v[counter2] * newVp.v[counter2];
				}
				beta /= yy;
				for( counter2 = 0; counter2 < n; counter2++ ) {
					d[counter2] *= beta;
				}

				for( counter1 = s - 1; counter1 >= 0; counter1-- ) {
					now = (VectorPair)vp.get( counter1 );
					beta = 0;
					// compute beta
					for( counter2 = 0; counter2 < n; counter2++ ) {
						beta += now.v[counter2] * d[counter2];
					}
					beta *= now.rho;
					// compute d
					for( counter2 = 0; counter2 < n; counter2++ ) {
						d[counter2] += ( now.alpha - beta ) * now.s[counter2];
					}
				}
			}
		}
		return i;
	}

	/**
	 * This constant can be used to specify that the steepest descent should be
	 * used in the <code>optimize</code>-method.
	 */
	public static final byte STEEPEST_DESCENT = 16;

	/**
	 * This constant can be used to specify that conjugate gradients (by
	 * <i>Fletcher</i> and <i>Reeves</i>) should be used in the
	 * <code>optimize</code>-method.
	 */
	public static final byte CONJUGATE_GRADIENTS_FR = 17;

	/**
	 * This constant can be used to specify that conjugate gradients (by
	 * <i>Polak</i> and <i>Ribi&egrave;re</i> should be used in the
	 * <code>optimize</code>-method.
	 */
	public static final byte CONJUGATE_GRADIENTS_PRP = 18;

	/**
	 * This constant can be used to specify that the quasi-Newton method of
	 * <i>Davidon</i>-<i>Fletcher</i>-<i>Powell</i> should be used in the
	 * <code>optimize</code> -method.
	 */
	public static final byte QUASI_NEWTON_DFP = 19;

	/**
	 * This constant can be used to specify that the quasi-Newton method of
	 * <i>Broyden</i>-<i>Fletcher</i>-<i>Goldfarb</i>-<i>Shanno</i> should be
	 * used in the <code>optimize</code>-method.
	 */
	public static final byte QUASI_NEWTON_BFGS = 20;

	/**
	 * This method enables you to use all different implemented optimization
	 * algorithms by only one method. You just have to change the parameter
	 * <code>algorithm</code>.
	 * 
	 * @param algorithm
	 *            the algorithm that should be used, either you use 2 &lt; x &lt; 11 for
	 *            {@link Optimizer#limitedMemoryBFGS(DifferentiableFunction, double[], byte, TerminationCondition, double, StartDistanceForecaster, OutputStream, Time)}
	 *            with <code>m=x</code> or the constants of the class:
	 *            {@link Optimizer#STEEPEST_DESCENT}, {@link Optimizer#CONJUGATE_GRADIENTS_FR},
	 *            {@link Optimizer#CONJUGATE_GRADIENTS_PRP}, {@link Optimizer#QUASI_NEWTON_DFP},
	 *            {@link Optimizer#QUASI_NEWTON_BFGS}
	 * @param f
	 *            the function to be minimized
	 * @param currentValues
	 *            at the begin the start vector and at the end the minimum
	 *            vector
	 * @param terminationMode
	 *            the choice how to terminate the algorithm
	 * @param linEps
	 *            the bound for stopping the linear search
	 * @param startDistance
	 *            the initial step for the linear search
	 * @param out
	 *            the {@link OutputStream} for writing some information
	 *            about the iterations, if <code>null</code> nothing will be
	 *            written.
	 * 
	 * @return the number of iterations
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * @throws TerminationException
	 *             if the termination condition is unknown
	 * @throws DimensionException
	 *             if the dimension of the scope of the function and the array
	 *             of arguments do not match
	 * @throws IOException
	 *             if there was something wrong with the IO

	 * @see Optimizer#optimize(byte, DifferentiableFunction, double[],
	 *      TerminationCondition, double, StartDistanceForecaster,
	 *      OutputStream, Time)
	 */
	public static int optimize( byte algorithm, DifferentiableFunction f, double[] currentValues, TerminationCondition terminationMode,
			double linEps, StartDistanceForecaster startDistance, OutputStream out ) throws DimensionException,
			TerminationException,
			IOException,
			EvaluationException {
		return optimize( algorithm, f, currentValues, terminationMode, linEps, startDistance, SafeOutputStream.getSafeOutputStream( out ), Time.getTimeInstance( out ) );
	}

	/**
	 * This method enables you to use all different implemented optimization
	 * algorithms by only one method. You just have to change the parameter
	 * <code>algorithm</code>.
	 * 
	 * @param algorithm
	 *            the algorithm that should be used, either you use 2 &lt; x &lt; 11 for
	 *            {@link Optimizer#limitedMemoryBFGS(DifferentiableFunction, double[], byte, TerminationCondition, double, StartDistanceForecaster, OutputStream, Time)}
	 *            with <code>m=x</code> or the constants of the class:
	 *            {@link Optimizer#STEEPEST_DESCENT}, {@link Optimizer#CONJUGATE_GRADIENTS_FR},
	 *            {@link Optimizer#CONJUGATE_GRADIENTS_PRP}, {@link Optimizer#QUASI_NEWTON_DFP},
	 *            {@link Optimizer#QUASI_NEWTON_BFGS}
	 * @param f
	 *            the function to be minimized
	 * @param currentValues
	 *            at the begin the start vector and at the end the minimum
	 *            vector
	 * @param terminationMode
	 *            the choice how to terminate the algorithm
	 * @param linEps
	 *            the bound for stopping the linear search
	 * @param startDistance
	 *            the initial step for the linear search
	 * @param out
	 *            the {@link OutputStream} for writing some information
	 *            about the iterations, if <code>null</code> nothing will be
	 *            written.
	 * @param t
	 *            a {@link Time} instance used to measure the elapsed time
	 * 
	 * @return the number of iterations
	 * 
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 * @throws TerminationException
	 *             if the termination condition is unknown
	 * @throws DimensionException
	 *             if the dimension of the scope of the function and the array
	 *             of arguments do not match
	 * @throws IOException
	 *             if there was something wrong with the IO
	 */
	public static int optimize( byte algorithm, DifferentiableFunction f, double[] currentValues, TerminationCondition terminationMode,
			double linEps, StartDistanceForecaster startDistance, OutputStream out, Time t ) throws DimensionException,
			TerminationException,
			IOException,
			EvaluationException {
		SafeOutputStream myOut = SafeOutputStream.getSafeOutputStream( out );
		int counter1;
		switch( algorithm ) {
			case STEEPEST_DESCENT:
				myOut.writeln( "STEEPEST_DESCENT" );
				counter1 = Optimizer.steepestDescent( f, currentValues, terminationMode, linEps, startDistance, out, t );
				break;
			case CONJUGATE_GRADIENTS_FR:
				myOut.writeln( "CONJUGATE_GRADIENTS_FR" );
				counter1 = Optimizer.conjugateGradientsFR( f, currentValues, terminationMode, linEps, startDistance, out, t );
				break;
			case CONJUGATE_GRADIENTS_PRP:
				myOut.writeln( "CONJUGATE_GRADIENTS_PRP" );
				counter1 = Optimizer.conjugateGradientsPRP( f, currentValues, terminationMode, linEps, startDistance, out, t );
				break;
			case QUASI_NEWTON_DFP:
				myOut.writeln( "QUASI_NEWTON_DFP" );
				counter1 = Optimizer.quasiNewtonDFP( f, currentValues, terminationMode, linEps, startDistance, out, t );
				break;
			case QUASI_NEWTON_BFGS:
				myOut.writeln( "QUASI_NEWTON_BFGS" );
				counter1 = Optimizer.quasiNewtonBFGS( f, currentValues, terminationMode, linEps, startDistance, out, t );
				break;
			default:
				if( algorithm >= 3 && algorithm <= 10 ) {
					myOut.writeln( "lm-BFGS (n=" + algorithm + ")" );
					counter1 = Optimizer.limitedMemoryBFGS( f,
							currentValues,
							algorithm,
							terminationMode,
							linEps,
							startDistance,
							out,
							t );
				} else {
					throw new IllegalArgumentException( "The algorithm choice is impossible." );
				}
		}
		return counter1;
	}
}
