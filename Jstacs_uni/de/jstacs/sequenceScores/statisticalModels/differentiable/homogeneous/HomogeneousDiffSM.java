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

package de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractVariableLengthDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.SamplingDifferentiableStatisticalModel;

/**
 * This is the main class for all homogeneous {@link de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore}s.
 * 
 * @author Jens Keilwagen
 */
public abstract class HomogeneousDiffSM extends
		AbstractVariableLengthDiffSM implements SamplingDifferentiableStatisticalModel {

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link HomogeneousDiffSM} that models sequences of arbitrary
	 * length.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 */
	protected HomogeneousDiffSM(AlphabetContainer alphabets) {
		super(alphabets);
	}

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link HomogeneousDiffSM} that models sequences of a given
	 * length.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param length
	 *            the length of the modeled sequences
	 */
	protected HomogeneousDiffSM(AlphabetContainer alphabets, int length) {
		super(alphabets, length);
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link HomogeneousDiffSM} out of its XML representation.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	protected HomogeneousDiffSM(StringBuffer source)
			throws NonParsableException {
		super(source);
	}

	/**
	 * Returns the maximal used markov oder.
	 * 
	 * @return the maximal used markov oder
	 */
	public abstract byte getMaximalMarkovOrder();
	
	/**
	 * This method allows to initialize the instance with an uniform distribution.
	 * 
	 * @param freeParams a switch whether to take only free parameters or to take all
	 */
	public abstract void initializeUniformly( boolean freeParams );
}
