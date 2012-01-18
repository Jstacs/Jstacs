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

/**
 * This class defines a {@link LogPrior} that does not penalize any parameter.
 * 
 * @author Jens Keilwagen
 */
public class DoesNothingLogPrior extends LogPrior {

	/**
	 * As this prior does not penalize parameters and does not have any
	 * parameters itself, this class does not have a constructor, but provides a
	 * default instance in order to reduce memory consumption.
	 */
	public static final DoesNothingLogPrior defaultInstance = new DoesNothingLogPrior();

	private DoesNothingLogPrior() {}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#evaluateFunction(double[])
	 */
	public double evaluateFunction( double[] params ) {
		return 0;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior#addGradientFor(double[], double[])
	 */
	@Override
	public void addGradientFor( double[] params, double[] grad ) {}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#getDimensionOfScope()
	 */
	public int getDimensionOfScope() {
		return UNKNOWN;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior#getNewInstance()
	 */
	@Override
	public LogPrior getNewInstance() throws CloneNotSupportedException {
		return defaultInstance;
	}

	@Override
	@Deprecated
	public StringBuffer toXML() throws RuntimeException {
		throw new RuntimeException( "Impossible to encode this instance as xml." );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "";
	}

}
