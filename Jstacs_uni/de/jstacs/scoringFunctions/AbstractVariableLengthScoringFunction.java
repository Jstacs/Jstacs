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

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This abstract class implements some methods declared in {@link NormalizableScoringFunction} based on the declaration
 * of methods in {@link VariableLengthScoringFunction}.
 * 
 * It allows to score subsequences of arbitrary length. This {@link NormalizableScoringFunction} should be the
 * super class for non-motif {@link NormalizableScoringFunction}s like homogeneous Markov
 * models, cyclic Markov models, ... etc.
 * 
 * @author Jens Keilwagen
 */
public abstract class AbstractVariableLengthScoringFunction extends
		AbstractNormalizableScoringFunction implements VariableLengthScoringFunction {

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link VariableLengthScoringFunction} that models sequences of arbitrary
	 * length.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer} of this
	 *            {@link VariableLengthScoringFunction}
	 * 
	 * @see AbstractVariableLengthScoringFunction#AbstractVariableLengthScoringFunction(AlphabetContainer,
	 *      int)
	 */
	protected AbstractVariableLengthScoringFunction(AlphabetContainer alphabets) {
		this(alphabets, 0);
	}

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link VariableLengthScoringFunction} that models sequences of a given
	 * length.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer} of this
	 *            {@link VariableLengthScoringFunction}
	 * @param length
	 *            the length of the modeled sequences
	 */
	protected AbstractVariableLengthScoringFunction(AlphabetContainer alphabets,
			int length) {
		super(alphabets, length);
		if (!alphabets.isSimple() || !alphabets.isDiscrete()) {
			throw new IllegalArgumentException(
					"The AlphabetContainer has to be simple and discrete.");
		}
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link VariableLengthScoringFunction} out of its XML
	 * representation.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed.
	 */
	protected AbstractVariableLengthScoringFunction(StringBuffer source)
			throws NonParsableException {
		super(source);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * getLogNormalizationConstant()
	 */
	public double getLogNormalizationConstant() {
		return getLogNormalizationConstant(length);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * getLogPartialNormalizationConstant(int)
	 */
	public double getLogPartialNormalizationConstant(int parameterIndex)
			throws Exception {
		return getLogPartialNormalizationConstant(parameterIndex, length);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#getLogScore(de.jstacs.data
	 * .Sequence, int)
	 */
	public double getLogScoreFor(Sequence seq, int start) {
		if (length != 0) {
			return getLogScore(seq, start, length);
		} else {
			return getLogScore(seq, start, seq.getLength() - start);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#getLogScoreAndPartialDerivation
	 * (de.jstacs.data.Sequence, int, de.jstacs.utils.IntList,
	 * de.jstacs.utils.DoubleList)
	 */
	public double getLogScoreAndPartialDerivation(Sequence seq, int start,
			IntList indices, DoubleList dList) {
		if (length != 0) {
			return getLogScoreAndPartialDerivation(seq, start, length, indices,
					dList);
		} else {
			return getLogScoreAndPartialDerivation(seq, start, seq.getLength()
					- start, indices, dList);
		}
	}
}