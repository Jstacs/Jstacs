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

import de.jstacs.data.Sequence;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This is an interface for all {@link NormalizableScoringFunction}s that allow to score
 * subsequences of arbitrary length. This {@link NormalizableScoringFunction} should be the
 * super class for non-motif {@link NormalizableScoringFunction}s like homogeneous Markov
 * models, cyclic Markov models, ... etc.
 * 
 * @author Jens Keilwagen
 */
public interface VariableLengthScoringFunction extends NormalizableScoringFunction {

	/**
	 * This method returns the logarithm of the normalization constant for a given sequence
	 * length.
	 * 
	 * @param length
	 *            the sequence length
	 * 
	 * @return the logarithm of the normalization constant
	 * 
	 * @see NormalizableScoringFunction#getLogNormalizationConstant()
	 */
	public abstract double getLogNormalizationConstant(int length);

	/**
	 * This method returns the logarithm of the partial normalization constant for a given
	 * parameter index and a sequence length.
	 * 
	 * @param parameterIndex
	 *            the index of the parameter
	 * @param length
	 *            the sequence length
	 * 
	 * @return the logarithm of the partial normalization constant
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see NormalizableScoringFunction#getLogPartialNormalizationConstant(int)
	 */
	public abstract double getLogPartialNormalizationConstant(int parameterIndex,
			int length) throws Exception;

	/**
	 * This method computes the logarithm of the score for a given subsequence.
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * @param start
	 *            the start index in the {@link Sequence}
	 * @param length
	 *            the length of the {@link Sequence} beginning at <code>start</code>
	 * 
	 * @return the logarithm of the score for the subsequence
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getLogScoreFor(Sequence,
	 *      int)
	 */
	public abstract double getLogScore(Sequence seq, int start, int length);

	/**
	 * This method computes the logarithm of the score and the partial
	 * derivations for a given subsequence.
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * @param start
	 *            the start index in the {@link Sequence}
	 * @param length
	 *            the end index in the {@link Sequence}
	 * @param indices
	 *            an {@link IntList} of indices, after method invocation the
	 *            list should contain the indices i where
	 *            {@latex.inline $\\frac{\\partial \\log score(seq)}{\\partial \\lambda_i}$}
	 *            is not zero
	 * @param dList
	 *            a {@link DoubleList} of partial derivations, after method
	 *            invocation the list should contain the corresponding
	 *            {@latex.inline $\\frac{\\partial \\log score(seq)}{\\partial \\lambda_i}$}
	 *            that are not zero
	 * 
	 * @return the logarithm of the score
	 * 
	 * @see ScoringFunction#getLogScoreAndPartialDerivation(Sequence, int,
	 *      IntList, DoubleList)
	 */
	public abstract double getLogScoreAndPartialDerivation(Sequence seq,
			int start, int length, IntList indices, DoubleList dList);

	/**
	 * This method sets the hyperparameters for the model parameters by
	 * evaluating the given statistic. The statistic can be interpreted as
	 * follows: The model has seen a number of sequences. From these sequences
	 * it is only known how long (<code>length</code>) and how often (
	 * <code>weight</code>) they have been seen.
	 * 
	 * @param length
	 *            the non-negative lengths of the sequences
	 * @param weight
	 *            the non-negative weight for the corresponding sequence
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see de.jstacs.motifDiscovery.Mutable
	 */
	public abstract void setStatisticForHyperparameters(int[] length,
			double[] weight) throws Exception;
}