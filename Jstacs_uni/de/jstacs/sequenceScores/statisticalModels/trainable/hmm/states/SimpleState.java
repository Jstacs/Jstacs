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

import java.text.NumberFormat;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.SilentEmission;

/**
 * This class implements a {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.State} based on {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission} that allows to reuse {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission}s for different {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.State}s.
 * 
 * @author Jens Keilwagen
 */
public class SimpleState implements TrainableState {

	/**
	 * The emission that is internally used for scoring subsequences.
	 */
	protected Emission e;
	/**
	 * The name of the state.
	 */
	protected String name;
	/**
	 * A switch that decides whether to use the forward or the reverse complementary strand of a sequence.
	 */
	protected boolean forward;
	
	/**
	 * This is the constructor of a {@link SimpleState}.
	 * 
	 * @param e the emission that is internally used for scoring subsequences
	 * @param name the name of the state
	 * @param forward a switch that decides whether to use the forward or the reverse complementary strand of a sequence
	 */
	public SimpleState( Emission e, String name, boolean forward ) {
		this.e = e;
		this.name = name;
		this.forward = forward;
	}

	public String getGraphvizNodeOptions( double weight, double maxWeight, NumberFormat nf ) {
		String res = "shape="+e.getNodeShape(forward)+", label="+e.getNodeLabel( weight < 0 ? weight : ((maxWeight-weight)/maxWeight), name, nf );
		if(weight > 0){
			res += ", width=2,height=2,style=filled,fillcolor=\"0 0 "+((maxWeight-weight)/maxWeight)+"\"";
		}
		return res;
	}

	public double getLogScoreFor(int startPos, int endPos, Sequence seq)
			throws WrongLengthException, OperationNotSupportedException {
		return e.getLogProbFor( forward, startPos, endPos, seq );
	}

	public void addToStatistic( int startPos, int endPos, double weight, Sequence seq ) throws OperationNotSupportedException {
		e.addToStatistic( forward, startPos, endPos, weight, seq );
	}

	public boolean isSilent() {
		return e instanceof SilentEmission;
	}
	
	public String toString( NumberFormat nf ) {
		AlphabetContainer con = e.getAlphabetContainer();
		return "state " + name + ( con == null || !con.isDiscrete() ? "" : " " + (forward?"forward":"reverse")) + "\n" + e.toString(nf); 
	}

	
	public String getName() {
		return name;
	}
}
