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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class implements a {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.State} based on {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission} that allows to reuse {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission}s for different {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.State}s.
 * 
 * @author Jens Keilwagen
 */
public class SimpleDifferentiableState extends SimpleState implements DifferentiableState {
	
	/**
	 * This is the constructor of a {@link SimpleState}.
	 * 
	 * @param e the emission that is internally used for scoring subsequences
	 * @param name the name of the state
	 * @param forward a switch that decides whether to use the forward or the reverse complementary strand of a sequence
	 */
	public SimpleDifferentiableState( DifferentiableEmission e, String name, boolean forward ) {
		super( e, name, forward );
	}

	public double getLogScoreAndPartialDerivation( int startPos, int endPos, IntList indices, DoubleList partDer, Sequence seq ) throws WrongLengthException,
			OperationNotSupportedException {
		return ((DifferentiableEmission)e).getLogProbAndPartialDerivationFor( forward, startPos, endPos, indices, partDer, seq );
	}
}
