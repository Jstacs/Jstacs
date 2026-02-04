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
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableCombinedWrapperEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class implements a {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.State} based on {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission} that allows to reuse {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission}s for different {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.State}s.
 * 
 * @author Jens Keilwagen
 */
public class SimpleDifferentiableState extends SimpleState implements DifferentiableState {
	
	private int phase;
	
	/**
	 * This is the constructor of a {@link SimpleDifferentiableState}.
	 * 
	 * @param e the emission that is internally used for scoring subsequences
	 * @param name the name of the state
	 * @param forward a switch that decides whether to use the forward or the reverse complementary strand of a sequence
	 */
	public SimpleDifferentiableState( DifferentiableEmission e, String name, boolean forward ) {
		super( e, name, forward );
		determinePhase();
	}
	
	protected void determinePhase() {
		phase = e instanceof DifferentiableCombinedWrapperEmission ? -1 : -2;
		if( phase>-2 ) {
			if( name.equals("C_start") || name.equals("C_stop") ) {
				phase = 0;
			} else {
				char first = name.charAt(0);
				if(first =='A' || first =='D' || first=='C') {
					phase = Integer.parseInt(name.substring(2,3));
				}
			}
		}
	}

	public double getLogScoreFor( int startPos, int endPos, Sequence seq ) throws WrongLengthException,
		OperationNotSupportedException {
		setStartPhase();
		return super.getLogScoreFor( startPos, endPos, seq );
	}
	
	public double getLogScoreAndPartialDerivation( int startPos, int endPos, IntList indices, DoubleList partDer, Sequence seq ) throws WrongLengthException,
			OperationNotSupportedException {
		setStartPhase();
		return ((DifferentiableEmission)e).getLogProbAndPartialDerivationFor( forward, startPos, endPos, indices, partDer, seq );
	}
	
	protected void setStartPhase() {
		if( phase>=0 ) {
			((DifferentiableCombinedWrapperEmission)e).setStartPhase(phase);
		}
	}
	
	public int getPhase() {
		return phase;
	}
}