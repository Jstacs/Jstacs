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

package de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.io.NonParsableException;

/**
 * Class for a {@link LogPrior} that defines a Gaussian prior on the parameters
 * of a set of {@link de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel}s
 * and a set of class parameters. The variances <code>v[i]</code> for a
 * parameter <code>i</code> of a
 * {@link de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel}
 * <code>fun[j]</code> are determined from the base variance <code>v[j]</code>
 * as<br>
 * <code>v[i] = v[j]*funs[j].getSizeOfEventSpaceForRandomVariablesOfParameter(j)</code>
 * . <br>
 * The variances for the class parameters are used as defined by the user. The
 * mean parameters are set to 0 for the parameters of the
 * {@link de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel}s and to the
 * user-specified means for the class parameters.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class SeparateGaussianLogPrior extends SeparateLogPrior {

	/**
	 * The variances
	 */
	protected double[] vars2;

	/**
	 * The means
	 */
	protected double[] mus;

	/**
	 * Creates a new {@link SeparateGaussianLogPrior} from a set of base
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
	public SeparateGaussianLogPrior( double[] vars, double[] classVars, double[] classMus ) {
		super( vars, classVars, classMus );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SeparateGaussianLogPrior} out of its XML
	 * representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SeparateGaussianLogPrior} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see SeparateLogPrior#SeparateLogPrior(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public SeparateGaussianLogPrior( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.SeparateLogPrior#unset()
	 */
	@Override
	public void unset() {
		vars2 = null;
		mus = null;
	}

	private void computeVars() {
		int num = funs.length - ( freeParameters ? 1 : 0 );
		for( int i = 0; i < funs.length; i++ ) {
			num += funs[i].getNumberOfParameters();
		}
		vars2 = new double[num];
		mus = new double[num];
		num = 0;
		for( ; num < funs.length - ( freeParameters ? 1 : 0 ); num++ ) {
			vars2[num] = classVars[num];
			mus[num] = classMus[num];
		}
		for( int i = 0; i < funs.length; i++ ) {
			for( int j = 0; j < funs[i].getNumberOfParameters(); j++ ) {
				vars2[num] = vars[i] * funs[i].getSizeOfEventSpaceForRandomVariablesOfParameter( j );
				mus[num] = 0;
				num++;
			}
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior#addGradientFor(double[], double[])
	 */
	@Override
	public void addGradientFor( double[] params, double[] grad ) {
		if( vars2 == null ) {
			computeVars();
		}
		for( int j = 0; j < params.length; j++ ) {
			grad[j] -= ( params[j] - mus[j] ) / vars2[j];
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#evaluateFunction(double[])
	 */
	public double evaluateFunction( double[] x ) throws DimensionException, EvaluationException {
		if( vars2 == null ) {
			computeVars();
		}
		double prior = 0;
		for( int i = 0; i < x.length; i++ ) {
			prior -= ( ( x[i] - mus[i] ) * ( x[i] - mus[i] ) ) / vars2[i];
		}
		prior /= 2.0;
		return prior;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#getDimensionOfScope()
	 */
	public int getDimensionOfScope() {
		if( vars2 == null ) {
			return UNKNOWN;
		} else {
			return vars2.length;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "Separate Gaussian prior";
	}
}
