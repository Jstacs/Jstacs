/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.scoringFunctions;

import de.jstacs.SequenceScoringFunction;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This interface is the main part of any {@link de.jstacs.classifier.scoringFunctionBased.ScoreClassifier}.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public interface ScoringFunction extends SequenceScoringFunction {
	/**
	 * Indicates that the number of parameters of this {@link ScoringFunction}
	 * is not known (yet).
	 */
	public static final int UNKNOWN = -1;

	/**
	 * Creates a clone (deep copy) of the current {@link ScoringFunction}
	 * instance.
	 * 
	 * @return the cloned instance of the current {@link ScoringFunction}
	 * 
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning the
	 *             {@link ScoringFunction}
	 */
	public ScoringFunction clone() throws CloneNotSupportedException;

	/**
	 * This method creates the underlying structure of the
	 * {@link ScoringFunction}.
	 * 
	 * @param index
	 *            the index of the class the {@link ScoringFunction} models
	 * @param freeParams
	 *            indicates whether the (reduced) parameterization is used
	 * @param data
	 *            the samples
	 * @param weights
	 *            the weights of the sequences in the samples
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public void initializeFunction(int index, boolean freeParams,
			DataSet[] data, double[][] weights) throws Exception;

	/**
	 * This method initializes the {@link ScoringFunction} randomly. It has to
	 * create the underlying structure of the {@link ScoringFunction}.
	 * 
	 * @param freeParams
	 *            indicates whether the (reduced) parameterization is used
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public void initializeFunctionRandomly(boolean freeParams) throws Exception;

	/**
	 * Returns the logarithmic score for a {@link Sequence} <code>seq</code> and
	 * fills lists with the indices and the partial derivations.
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * @param indices
	 *            an {@link IntList} of indices, after method invocation the
	 *            list should contain the indices <code>i</code> where
	 *            {@latex.inline $\\frac{\\partial \\log score(seq)}{\\partial \\lambda_i}$}
	 *            is not zero
	 * @param partialDer
	 *            a {@link DoubleList} of partial derivations, after method
	 *            invocation the list should contain the corresponding
	 *            {@latex.inline $\\frac{\\partial \\log score(seq)}{\\partial \\lambda_i}$}
	 *            that are not zero
	 * 
	 * @return the logarithmic score for the {@link Sequence}
	 */
	public double getLogScoreAndPartialDerivation(Sequence seq, IntList indices, DoubleList partialDer);

	/**
	 * Returns the logarithmic score for a {@link Sequence} beginning at
	 * position <code>start</code> in the {@link Sequence} and fills lists with
	 * the indices and the partial derivations.
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * @param start
	 *            the start position in the {@link Sequence}
	 * @param indices
	 *            an {@link IntList} of indices, after method invocation the
	 *            list should contain the indices i where
	 *            {@latex.inline $\\frac{\\partial \\log score(seq)}{\\partial \\lambda_i}$}
	 *            is not zero
	 * @param partialDer
	 *            a {@link DoubleList} of partial derivations, after method
	 *            invocation the list should contain the corresponding
	 *            {@latex.inline $\\frac{\\partial \\log score(seq)}{\\partial \\lambda_i}$}
	 *            that are not zero
	 * 
	 * @return the logarithmic score for the {@link Sequence}
	 */
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer);

	/**
	 * Returns the logarithmic score for a {@link Sequence} beginning at
	 * position <code>start</code> in the {@link Sequence} and fills lists with
	 * the indices and the partial derivations.
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * @param start
	 *            the start position in the {@link Sequence}
	 * @param indices
	 *            an {@link IntList} of indices, after method invocation the
	 *            list should contain the indices i where
	 *            {@latex.inline $\\frac{\\partial \\log score(seq)}{\\partial \\lambda_i}$}
	 *            is not zero
	 * @param partialDer
	 *            a {@link DoubleList} of partial derivations, after method
	 *            invocation the list should contain the corresponding
	 *            {@latex.inline $\\frac{\\partial \\log score(seq)}{\\partial \\lambda_i}$}
	 *            that are not zero
	 * 
	 * @return the logarithmic score for the {@link Sequence}
	 */
	public double getLogScoreAndPartialDerivation( Sequence seq, int start, int end, IntList indices, DoubleList partialDer ) throws Exception;

	/**
	 * Returns the number of parameters in this {@link ScoringFunction}. If the
	 * number of parameters is not known yet, the method returns
	 * {@link #UNKNOWN}.
	 * 
	 * @return the number of parameters in this {@link ScoringFunction}
	 * 
	 * @see ScoringFunction#UNKNOWN
	 */
	public int getNumberOfParameters();

	/**
	 * This method returns the number of recommended optimization starts. The
	 * standard implementation returns 1.
	 * 
	 * @return the number of recommended optimization starts
	 */
	public int getNumberOfRecommendedStarts();

	/**
	 * Returns a <code>double</code> array of dimension
	 * {@link #getNumberOfParameters()} containing the current parameter values.
	 * If one likes to use these parameters to start an optimization it is
	 * highly recommended to invoke
	 * {@link #initializeFunction(int, boolean, DataSet[], double[][])} before.
	 * After an optimization this method can be used to get the current
	 * parameter values.
	 * 
	 * @return the current parameter values
	 * 
	 * @throws Exception
	 *             if no parameters exist (yet)
	 */
	public double[] getCurrentParameterValues() throws Exception;

	/**
	 * This method sets the internal parameters to the values of
	 * <code>params</code> between <code>start</code> and
	 * <code>start + {@link #getNumberOfParameters()} - 1</code>
	 * 
	 * @param params
	 *            the new parameters
	 * @param start
	 *            the start index in <code>params</code>
	 */
	public void setParameters(double[] params, int start);

	/**
	 * Returns the initial class parameter for the class this
	 * {@link ScoringFunction} is responsible for, based on the class
	 * probability <code>classProb</code>.
	 * 
	 * @param classProb
	 *            the class probability
	 * 
	 * @return the initial class parameter
	 */
	public double getInitialClassParam(double classProb);
}
