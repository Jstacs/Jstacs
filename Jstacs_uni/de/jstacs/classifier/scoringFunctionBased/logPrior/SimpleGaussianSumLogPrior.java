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
import de.jstacs.io.XMLParser;

/**
 * This class implements a prior that is a product of Gaussian distributions
 * with mean 0 and equal variance for each parameter.
 * 
 * @author Jens Keilwagen
 */
public class SimpleGaussianSumLogPrior extends LogPrior {

	private double sigmaSq;

	/**
	 * Creates a new {@link SimpleGaussianSumLogPrior} with mean 0 and variance
	 * <code>sigma</code> for all parameters, including the class parameters.
	 * 
	 * @param sigma
	 *            the variance
	 */
	public SimpleGaussianSumLogPrior( double sigma ) {
		sigmaSq = sigma * sigma;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SimpleGaussianSumLogPrior} out of its XML
	 * representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SimpleGaussianSumLogPrior} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public SimpleGaussianSumLogPrior( StringBuffer xml ) throws NonParsableException {
		StringBuffer content = XMLParser.extractForTag( xml, "SimpleGaussianProductPrior" );
		sigmaSq = XMLParser.extractObjectForTags( content, "sigmaSquared", double.class );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.logPrior.LogPrior#addGradientFor(double[], double[])
	 */
	@Override
	public void addGradientFor( double[] params, double[] grad ) {
		for( int i = 0; i < params.length; i++ ) {
			grad[i] -= params[i] / sigmaSq;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#evaluateFunction(double[])
	 */
	public double evaluateFunction( double[] params ) {
		// only relevant terms
		double erg = 0;
		for( int i = 0; i < params.length; i++ ) {
			erg -= params[i] * params[i];
		}
		return erg / ( 2d * sigmaSq );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#getDimensionOfScope()
	 */
	public int getDimensionOfScope() {
		return UNKNOWN;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.logPrior.LogPrior#getNewInstance()
	 */
	@Override
	public SimpleGaussianSumLogPrior getNewInstance() throws CloneNotSupportedException {
		SimpleGaussianSumLogPrior newI = new SimpleGaussianSumLogPrior( 0 );
		newI.sigmaSq = this.sigmaSq;
		return newI;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.logPrior.LogPrior#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 300 );
		XMLParser.appendObjectWithTags( xml, sigmaSq, "sigmaSquared" );
		XMLParser.addTags( xml, "SimpleGaussianProductPrior" );
		return xml;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.logPrior.LogPrior#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "Gaussian product prior";
	}
}
