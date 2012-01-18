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

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This abstract class implements some methods declared in {@link DifferentiableStatisticalModel} based on the declaration
 * of methods in {@link VariableLengthDiffSM}.
 * 
 * It allows to score subsequences of arbitrary length. This {@link DifferentiableStatisticalModel} should be the
 * super class for non-motif {@link DifferentiableStatisticalModel}s like homogeneous Markov
 * models, cyclic Markov models, ... etc.
 * 
 * @author Jens Keilwagen
 */
public abstract class AbstractVariableLengthDiffSM extends
		AbstractDifferentiableStatisticalModel implements VariableLengthDiffSM {

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link VariableLengthDiffSM} that models sequences of arbitrary
	 * length.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer} of this
	 *            {@link VariableLengthDiffSM}
	 * 
	 * @see AbstractVariableLengthDiffSM#AbstractVariableLengthDiffSM(AlphabetContainer,
	 *      int)
	 */
	protected AbstractVariableLengthDiffSM(AlphabetContainer alphabets) {
		this(alphabets, 0);
	}

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link VariableLengthDiffSM} that models sequences of a given
	 * length.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer} of this
	 *            {@link VariableLengthDiffSM}
	 * @param length
	 *            the length of the modeled sequences
	 */
	protected AbstractVariableLengthDiffSM(AlphabetContainer alphabets,
			int length) {
		super(alphabets, length);
		if (!alphabets.isSimple() || !alphabets.isDiscrete()) {
			throw new IllegalArgumentException(
					"The AlphabetContainer has to be simple and discrete.");
		}
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link VariableLengthDiffSM} out of its XML
	 * representation.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed.
	 */
	protected AbstractVariableLengthDiffSM(StringBuffer source)
			throws NonParsableException {
		super(source);
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogNormalizationConstant()
	 */
	@Override
	public double getLogNormalizationConstant() {
		return getLogNormalizationConstant(length);
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogPartialNormalizationConstant(int)
	 */
	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex)
			throws Exception {
		return getLogPartialNormalizationConstant(parameterIndex, length);
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getLogScoreFor(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		if (length != 0) {
			return getLogScoreFor( seq, start, start+length-1 );
		} else {
			return getLogScoreFor( seq, start, seq.getLength()-1 );
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScoreAndPartialDerivation(de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList dList) {
		if (length != 0) {
			return getLogScoreAndPartialDerivation( seq, start, start+length-1, indices, dList );
		} else {
			return getLogScoreAndPartialDerivation( seq, start, seq.getLength()-1, indices, dList );
		}
	}
	
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
}