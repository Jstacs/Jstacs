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
