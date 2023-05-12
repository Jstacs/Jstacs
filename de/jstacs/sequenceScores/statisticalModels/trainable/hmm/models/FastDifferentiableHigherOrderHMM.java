package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.SimpleDifferentiableState;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.MaxHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement;

/**
 * DifferentiableHigherOrderHMM that tries to compute the likelihood and gradient faster.
 * This especially reasonable if the same emissions are used in several states.
 * 
 * @author Jens Keilwagen
 */
public class FastDifferentiableHigherOrderHMM extends DifferentiableHigherOrderHMM {

	RefSimpleDifferentiableState[] dStates;
	DifferentiableEmission[] dEmission;
	
	public FastDifferentiableHigherOrderHMM(MaxHMMTrainingParameterSet trainingParameterSet, String[] name,
			int[] emissionIdx, DifferentiableEmission[] emission, double ess,
			TransitionElement... te) throws Exception {
		super(trainingParameterSet, name, emissionIdx, null, emission, ess, te);
	}

	public FastDifferentiableHigherOrderHMM(String type, int[][] statesGroups,
			MaxHMMTrainingParameterSet trainingParameterSet, String[] name, int[] emissionIdx,
			DifferentiableEmission[] emission, double ess, TransitionElement... te) throws Exception {
		super(type, statesGroups, trainingParameterSet, name, emissionIdx, null, emission, ess, te);
	}

	public FastDifferentiableHigherOrderHMM(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	protected void createStates() {
		this.states = dStates = new RefSimpleDifferentiableState[emissionIdx.length];
		dEmission = new DifferentiableEmission[emission.length];
		for( int i = 0; i < emission.length; i++ ) {
			dEmission[i] = (DifferentiableEmission) emission[i];
		}
		for( int i = 0; i < emissionIdx.length; i++ ) {
			dStates[i] = new RefSimpleDifferentiableState( emissionIdx[i], dEmission[emissionIdx[i]], name[i], forward[i] );
		}
	}
	
	protected void fillLogEmission( int startPos, Sequence seq ) throws OperationNotSupportedException, WrongLengthException {
		for( int e = 0; e < dEmission.length; e++ ) {
			logEmission[e] = dEmission[e].getLogProbFor(true, startPos, startPos, seq); 
		}
	}
	
	protected void fillLogEmissionAndPartialDer( int endPos, Sequence seq ) throws OperationNotSupportedException, WrongLengthException {
		for( int e = 0; e < dEmission.length; e++ ) {
			indicesState[e].clear();
			partDerState[e].clear();
			logEmission[e] = dEmission[e].getLogProbAndPartialDerivationFor(true, endPos, endPos, indicesState[e], partDerState[e], seq );
		}
	}

	protected int getIndex( int i ) {
		return dStates[i].idx;
	}
	
	/**
	 * This state is used for faster computation of likelihood and gradient.
	 * 
	 * @author Jens Keilwagen
	 *
	 * @see FastDifferentiableHigherOrderHMM#getIndex(int)
	 */
	private static class RefSimpleDifferentiableState extends SimpleDifferentiableState {
		int idx;
		public RefSimpleDifferentiableState( int idx, DifferentiableEmission e, String name, boolean forward ) {
			super( e, name, forward );
			this.idx=idx;
		}
	}
}
