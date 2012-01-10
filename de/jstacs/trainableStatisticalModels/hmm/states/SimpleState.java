package de.jstacs.trainableStatisticalModels.hmm.states;

import java.text.NumberFormat;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;
import de.jstacs.data.WrongLengthException;
import de.jstacs.trainableStatisticalModels.hmm.states.emissions.Emission;
import de.jstacs.trainableStatisticalModels.hmm.states.emissions.SilentEmission;

/**
 * This class implements a {@link de.jstacs.trainableStatisticalModels.hmm.State} based on {@link de.jstacs.trainableStatisticalModels.hmm.states.emissions.Emission} that allows to reuse {@link de.jstacs.trainableStatisticalModels.hmm.states.emissions.Emission}s for different {@link de.jstacs.trainableStatisticalModels.hmm.State}s.
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
	
	public String toString() {
		AlphabetContainer con = e.getAlphabetContainer();
		return "state " + name + ( con == null || !con.isDiscrete() ? "" : " " + (forward?"forward":"reverse")) + "\n" + e; 
	}
}
