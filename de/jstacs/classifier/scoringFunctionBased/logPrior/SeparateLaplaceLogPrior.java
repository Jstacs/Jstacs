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

package de.jstacs.classifier.scoringFunctionBased.logPrior;

import de.jstacs.NonParsableException;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;

/**
 * Class for a {@link LogPrior} that defines a Laplace prior on the parameters
 * of a set of {@link de.jstacs.scoringFunctions.NormalizableScoringFunction}s
 * and a set of class parameters. The scale hyperparameters of the Laplace prior
 * are determined from a set of class-specific base variances and the class
 * variances, respectively. The formula to determine the shape hyperparameter
 * <code>b[i]</code> for a parameter <code>i</code> from a variance
 * <code>v[i]</code> is<br />
 * <code>b[i] = Math.sqrt(v[i]/2)</code>.<br />
 * The variances <code>v[i]</code> for a parameter <code>i</code> of a
 * {@link de.jstacs.scoringFunctions.NormalizableScoringFunction}
 * <code>fun[j]</code> are determined from the base variance <code>v[j]</code>
 * as<br />
 * <code>v[i] = v[j]*funs[j].getSizeOfEventSpaceForRandomVariablesOfParameter(j)</code>
 * .<br />
 * The mean parameters are set to 0 for the parameters of the
 * {@link de.jstacs.scoringFunctions.NormalizableScoringFunction}s and to the
 * user-specified means for the class parameters.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class SeparateLaplaceLogPrior extends SeparateLogPrior {

	private double[] bs;

	private double[] mus;

	/**
	 * Creates a new {@link SeparateLaplaceLogPrior} from a set of base
	 * variances <code>vars</code>, a set of class variances
	 * <code>classVars</code> and a set of class means <code>classMus</code>.
	 * 
	 * @param vars
	 *            the base variances for each class
	 * @param classVars
	 *            the class variances
	 * @param classMus
	 *            the class means
	 * 
	 * @see SeparateLogPrior#SeparateLogPrior(double[], double[], double[])
	 */
	public SeparateLaplaceLogPrior( double[] vars, double[] classVars, double[] classMus ) {
		super( vars, classVars, classMus );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SeparateLaplaceLogPrior} out of its XML
	 * representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SeparateLaplaceLogPrior} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see SeparateLogPrior#SeparateLogPrior(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public SeparateLaplaceLogPrior( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.logPrior.SeparateLogPrior#unset()
	 */
	@Override
	protected void unset() {
		this.bs = null;
		this.mus = null;
	}

	private void computeBs() {
		int num = funs.length - ( freeParameters ? 1 : 0 );
		for( int i = 0; i < funs.length; i++ ) {
			num += funs[i].getNumberOfParameters();
		}
		bs = new double[num];
		mus = new double[num];
		num = 0;
		for( ; num < funs.length - ( freeParameters ? 1 : 0 ); num++ ) {
			bs[num] = Math.sqrt( classVars[num] / 2d );
			mus[num] = classMus[num];
		}
		for( int i = 0; i < funs.length; i++ ) {
			for( int j = 0; j < funs[i].getNumberOfParameters(); j++ ) {
				bs[num] = Math.sqrt( vars[i] * funs[i].getSizeOfEventSpaceForRandomVariablesOfParameter( j ) / 2d );
				mus[num] = 0;
				num++;
			}
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.logPrior.LogPrior#addGradientFor(double[], double[])
	 */
	@Override
	public void addGradientFor( double[] params, double[] vector ) {
		if( bs == null ) {
			computeBs();
		}
		for( int j = 0; j < params.length; j++ ) {
			vector[j] -= Math.signum( params[j] - mus[j] ) / bs[j];
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#evaluateFunction(double[])
	 */
	public double evaluateFunction( double[] x ) throws DimensionException, EvaluationException {
		if( bs == null ) {
			computeBs();
		}
		double prior = 0;
		for( int i = 0; i < x.length; i++ ) {
			prior -= Math.abs( x[i] - mus[i] ) / ( bs[i] );
		}
		return prior;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#getDimensionOfScope()
	 */
	public int getDimensionOfScope() {
		if( bs == null ) {
			return UNKNOWN;
		} else {
			return bs.length;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.logPrior.LogPrior#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "Separate Laplace prior";
	}
}
