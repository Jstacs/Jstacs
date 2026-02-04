package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models;

import java.util.ArrayList;
import java.util.HashMap;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.SimpleDifferentiableState;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableCombinedWrapperEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.filter.Filter;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.MaxHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement;
import de.jstacs.utils.IntList;

/**
 * DifferentiableHigherOrderHMM that tries to compute the likelihood and gradient faster.
 * This especially reasonable if the same emissions are used in several states.
 * 
 * @author Jens Keilwagen
 */
public class FastDifferentiableHigherOrderHMM extends DifferentiableHigherOrderHMM {

	RefSimpleDifferentiableState[] dStates;
	DifferentiableEmission[] dEmission;
	int[][] combis;
	
	public FastDifferentiableHigherOrderHMM(MaxHMMTrainingParameterSet trainingParameterSet, String[] name,
			int[] emissionIdx, DifferentiableEmission[] emission, double ess,
			TransitionElement... te) throws Exception {
		this(null, null, trainingParameterSet, name, null, emissionIdx, emission, ess, null, te);
	}

	public FastDifferentiableHigherOrderHMM( String type, int[][] statesGroups, 
			MaxHMMTrainingParameterSet trainingParameterSet, String[] name, Filter[] filter, int[] emissionIdx,
			DifferentiableEmission[] emission, double ess, int[] transIndex, TransitionElement... te) throws Exception {
		super(type, statesGroups, trainingParameterSet, name, filter, emissionIdx, null, emission, ess, transIndex, te);
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
		IntList phase = new IntList();
		HashMap<String,Integer> hash = new HashMap<String, Integer>();
		ArrayList<int[]> combis = new ArrayList<int[]>();
		for( int i = 0; i < emissionIdx.length; i++ ) {
			dStates[i] = new RefSimpleDifferentiableState( -100, dEmission[emissionIdx[i]], name[i], forward[i] );
			int[] combi = { emissionIdx[i], dStates[i].getPhase() };
			String key = combi[0] + "_" + combi[1];
			Integer idx = hash.get(key);
			if( idx == null ) {
				idx=combis.size();
				hash.put(key, idx);
				combis.add( combi );
			}
			dStates[i].idx = idx;
		}
		this.combis = combis.toArray( new int[0][] );
	}
	
	protected void fillLogEmission( int startPos, Sequence seq ) throws OperationNotSupportedException, WrongLengthException {
		fillLogEmissionAndPartialDer(startPos, seq, false);
	}
	
	protected void fillLogEmissionAndPartialDer( int endPos, Sequence seq, boolean grad ) throws OperationNotSupportedException, WrongLengthException {
		boolean forward = true; //XXX
		for( int c = 0; c < combis.length; c++ ) {
			if( combis[c][1]>=0 ) {
				((DifferentiableCombinedWrapperEmission) dEmission[combis[c][0]]).setStartPhase(combis[c][1]);
			}
			if( grad ) {
				indicesState[c].clear();
				partDerState[c].clear();
				logEmission[c] = dEmission[combis[c][0]].getLogProbAndPartialDerivationFor(forward, endPos, endPos, indicesState[c], partDerState[c], seq );
			} else {
				logEmission[c] = dEmission[combis[c][0]].getLogProbFor(forward, endPos, endPos, seq );
			}
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
