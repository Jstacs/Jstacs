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
 * \( f: \mathbb{R}^n \to \mathbb{R} \). The gradient is computed numerically. Each partial
 * differentiation is computed by the formula<br>
 * 
 * \( \partial_k f(x) = \frac{f(x)-f(x+\varepsilon*e_k)}{\varepsilon} \).
 * 
 * Be careful with the choice of epsilon, too small values might result in precision problems
 * 
 * @author Jens Keilwagen
 */
public class NumericalDifferentiableFunction extends DifferentiableFunction {

	/**
	 * The constant used in the computation of the gradient. Should be close to
	 * 0 but not exactly 0.
	 */
	protected double eps;
	
	/**
	 * The function to be differentiated numerically
	 */
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
		/*
		//(too) simple approximation
		double current = evaluateFunction( x ), h;
		//System.out.println(current);
		for( int i = 0; i < n; i++ ) {
			h = x[i];
			x[i] += eps;
			double fx = evaluateFunction( x );
			gradient[i] = ( fx - current ) / eps;
			//System.out.println( i + "\t" + fx + "\t" + (fx-current) + "\t" + gradient[i] + "\t" + Arrays.toString(x));
			x[i] = h;
		}/**/
		
		//better approximation
		for( int i = 0; i < n; i++ ) {
			double h = x[i];
			
			x[i] = h + eps;
			double fx1 = evaluateFunction( x );
			x[i] = h - eps;
			double fx2 = evaluateFunction( x );
			gradient[i] = ( fx1 - fx2 ) / (2*eps);
			//System.out.println( i + "\t" + fx1 + "\t" + fx2 + "\t" + gradient[i] );
			
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
		
	/**
	 * Compares the analytical and numerical gradient.
	 *  
	 * @param df the differentiable function
	 * @param x the parameters
	 * @param eps the epsilon for the numerical gradient
	 * 
	 * @throws DimensionException if the dimension of parameters does not match the function
	 * @throws EvaluationException if the evaluation of the function does not work
	 */
	public static void compare( DifferentiableFunction df, double[] x, double eps ) throws DimensionException, EvaluationException {
		NumericalDifferentiableFunction nFun = new NumericalDifferentiableFunction(df, eps);

		double[] anaGrad = df.evaluateGradientOfFunction(x);
		double[] numGrad = nFun.evaluateGradientOfFunction(x);
		
		double max = Double.NEGATIVE_INFINITY, maxPerc=0;
		System.out.println( "#\tparameter\tanalytical gradient\tnumerical gradient\tabs difference\t|%|" );
		for( int i = 0; i < anaGrad.length; i++ ) {
			double m = Math.abs(anaGrad[i]-numGrad[i]);
			double p = m/Math.abs(anaGrad[i])*100;
			System.out.println( i + "\t" + x[i] + "\t" + anaGrad[i] + "\t" + numGrad[i] + "\t" + m + "\t" + p );
			if( m > max ) {
				max = m;
			}
			if( p > maxPerc ) {
				maxPerc = p;
			}
		}
		System.out.println("<max>\t\t\t\t"+max+"\t"+maxPerc);
	}
}
