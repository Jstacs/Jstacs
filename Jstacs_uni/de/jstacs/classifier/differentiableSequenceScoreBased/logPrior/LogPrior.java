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

package de.jstacs.classifier.differentiableSequenceScoreBased.logPrior;

import de.jstacs.Storable;
import de.jstacs.algorithms.optimization.DifferentiableFunction;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;

/**
 * The abstract class for any log-prior used e.g. for maximum supervised
 * posterior optimization.
 * 
 * @author Jens Keilwagen
 */
public abstract class LogPrior extends DifferentiableFunction implements Storable {

	/**
	 * Sometimes the method
	 * {@link de.jstacs.algorithms.optimization.Function#getDimensionOfScope()}
	 * can not return a value. In these cases it is recommended to use this
	 * constant.
	 */
	public static final int UNKNOWN = -1;

	/**
	 * Adds the gradient of the log-prior using the current parameters to a
	 * given vector.
	 * 
	 * @param params
	 *            the parameters
	 * @param vector
	 *            the vector
	 * 
	 * @throws EvaluationException
	 *             if the gradient could not be evaluated
	 */
	public abstract void addGradientFor( double[] params, double[] vector ) throws EvaluationException;

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.DifferentiableFunction#evaluateGradientOfFunction(double[])
	 */
	@Override
	public final double[] evaluateGradientOfFunction( double[] params ) throws EvaluationException {
		double[] grad = new double[getDimensionOfScope()];
		addGradientFor( params, grad );
		return grad;
	}

	/**
	 * This method returns an empty new instance of the current prior. The
	 * method works similar to clone, but does not clone {@link DifferentiableSequenceScore}
	 * s that may be inside the instance. The {@link DifferentiableSequenceScore}s must be
	 * set by an invocation of the method
	 * {@link #set(boolean, DifferentiableSequenceScore...)}.
	 * 
	 * @return a new empty instance of the prior
	 * 
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * 
	 * @see LogPrior#set(boolean, DifferentiableSequenceScore...)
	 */
	public abstract LogPrior getNewInstance() throws CloneNotSupportedException;

	/**
	 * Resets all pre-computed values to their initial values using the
	 * {@link DifferentiableSequenceScore}s <code>funs</code>.
	 * 
	 * @param freeParameters
	 *            the switch for using only the free parameters or all
	 *            parameters in a {@link DifferentiableSequenceScore}
	 * @param funs
	 *            the {@link DifferentiableSequenceScore}s for the prior
	 * 
	 * @throws Exception
	 *             if the {@link DifferentiableSequenceScore}s could not be set
	 */
	public void set( boolean freeParameters, DifferentiableSequenceScore... funs ) throws Exception {

	}

	/**
	 * Encodes the prior as an XML representation. It does not encode all
	 * information, since the method {@link #set(boolean, DifferentiableSequenceScore...)}
	 * has to be invoked after decoding.
	 * 
	 * @return the XML representation;
	 */
	public abstract StringBuffer toXML();

	/**
	 * Returns a <b>short</b> instance name.
	 * 
	 * @return a <b>short</b> instance name
	 */
	public abstract String getInstanceName();
}
