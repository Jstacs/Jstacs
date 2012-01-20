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

package de.jstacs.sequenceScores.statisticalModels.differentiable;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This is an interface for all {@link DifferentiableStatisticalModel}s that allow to score
 * subsequences of arbitrary length. This {@link DifferentiableStatisticalModel} should be the
 * super class for non-motif {@link DifferentiableStatisticalModel}s like homogeneous Markov
 * models, cyclic Markov models, ... etc.
 * 
 * @author Jens Keilwagen
 */
public interface VariableLengthDiffSM extends DifferentiableStatisticalModel {

	/**
	 * This method returns the logarithm of the normalization constant for a given sequence
	 * length.
	 * 
	 * @param length
	 *            the sequence length
	 * 
	 * @return the logarithm of the normalization constant
	 * 
	 * @see DifferentiableStatisticalModel#getLogNormalizationConstant()
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
	 * @see DifferentiableStatisticalModel#getLogPartialNormalizationConstant(int)
	 */
	public abstract double getLogPartialNormalizationConstant(int parameterIndex,
			int length) throws Exception;

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScore(de.jstacs.data.Sequence)
	 */
	@Override
	public abstract double getLogScoreFor( Sequence seq, int startpos, int endpos );
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScoreAndPartialDerivation(int, de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public abstract double getLogScoreAndPartialDerivation( Sequence seq, int startpos, int endpos, IntList indices, DoubleList partialDer );
	
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