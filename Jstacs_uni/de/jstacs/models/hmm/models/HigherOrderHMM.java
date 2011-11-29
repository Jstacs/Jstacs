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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.models.hmm.models;

import java.util.Arrays;

import javax.naming.OperationNotSupportedException;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.WrongLengthException;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;
import de.jstacs.models.hmm.AbstractHMM;
import de.jstacs.models.hmm.HMMTrainingParameterSet;
import de.jstacs.models.hmm.State;
import de.jstacs.models.hmm.Transition;
import de.jstacs.models.hmm.states.SimpleState;
import de.jstacs.models.hmm.states.TrainableState;
import de.jstacs.models.hmm.states.emissions.Emission;
import de.jstacs.models.hmm.training.BaumWelchParameterSet;
import de.jstacs.models.hmm.training.MaxHMMTrainingParameterSet;
import de.jstacs.models.hmm.training.ViterbiParameterSet;
import de.jstacs.models.hmm.transitions.TrainableTransition;
import de.jstacs.models.hmm.transitions.BasicHigherOrderTransition.AbstractTransitionElement;
import de.jstacs.models.mixture.AbstractMixtureModel;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.RealTime;
import de.jstacs.utils.ToolBox;

/**
 * This class implements a higher order hidden Markov model.
 * Currently, the modeling of the transitions is higher order, but is easily possible to extend this to emissions.
 * 
 * <br>
 * <a name="final">
 * This implementation allows to have a set of final states {@latex.inline $S_F$}.
 * A state is denoted final states if it is allowed at the end of a path. Hence, any valid path always ends with a final state.
 * Using the method {@link #getLogProbFor(Sequence)} for sequence {@latex.inline $\\underline{x}$} returns the value
 * 
 * {@latex.ilb \\[
 * P(\\underline{x}, u_F \\in S_F | \\underline{\\lambda}) = \\sum_{i \\in S_F} P(\\underline{x}, u_F = i | \\underline{\\lambda})
 * .\\]}
 * 
 * Setting {@latex.inline $S_F$} to all states leads to the computation of the likelihood.
 * </a>
 * 
 * @author Jens Keilwagen
 */
public class HigherOrderHMM extends AbstractHMM {

	//helper variables
	/**
	 * Helper variable = only for internal use. This array is used for {@link Transition#fillTransitionInformation(int, int, int, int[])}.
	 */
	protected int[][] container;
	/**
	 * Helper variable = only for internal use. This array is used for compute the emissions at each position of a sequence only once,
	 * which might be beneficial for higher order models.
	 * 
	 * @see HigherOrderHMM#emission
	 */
	protected double[][] logEmission;
	
	/**
	 * Helper variable = only for internal use. This array is used to compute the forward matrix. It stores intermediate results.
	 * 
	 * #see {@link #numberOfSummands}
	 */
	private double[][][][] forwardIntermediate;

	/**
	 * Helper variable = only for internal use. This array is used to compute the backward matrix. It stores intermediate results.
	 * 
	 * #see {@link #numberOfSummands}
	 */
	protected double[][] backwardIntermediate;
	
	/**
	 * Helper variable = only for internal use. This array is used to compute the forward and backward matrix. It stores the number of intermediate results.
	 */
	protected int[][][] numberOfSummands;
	
	/**
	 * Helper variable = only for internal use. This field is used in the method {@link #samplePath(IntList, int, int, Sequence)}.
	 */
	protected IntList[] stateList;
	protected boolean skipInit;
	
	/**
	 * This is a convenience constructor. It assumes that state <code>i</code> used emission <code>i</code> on the forward strand.
	 * 
	 * @param trainingParameterSet the {@link de.jstacs.parameters.ParameterSet} that determines the training algorithm and contains the necessary {@link de.jstacs.parameters.Parameter}s
	 * @param name the names of the states
	 * @param emission the emissions
	 * @param te the {@link AbstractTransitionElement}s building a transition
	 * 
	 * @throws Exception if 
	 * 	<ul>
	 *  <li>some component could not be cloned</li> 
	 *  <li>some the length of <code>name, emissionIdx,</code> or <code>forward</code> is not equal to the number of states</li>
	 *  <li>not all emissions use the same {@link de.jstacs.data.AlphabetContainer}</li>
	 *  <li>the states can not be handled by the transition
	 *  </ul>
	 * 
	 * @see de.jstacs.models.hmm.models.HigherOrderHMM#HigherOrderHMM(HMMTrainingParameterSet, String[], int[], boolean[], Emission[], BasicHigherOrderTransition.AbstractTransitionElement...)
	 */
	public HigherOrderHMM( int threads, HMMTrainingParameterSet trainingParameterSet, String[] name, Emission[] emission, AbstractTransitionElement... te ) throws Exception {
		this( threads, trainingParameterSet, name, null, null, emission, te );
	}
	
	/**
	 * This is the main constructor.
	 * 
	 * @param trainingParameterSet the {@link de.jstacs.parameters.ParameterSet} that determines the training algorithm and contains the necessary {@link de.jstacs.parameters.Parameter}s
	 * @param name the names of the states
	 * @param emissionIdx the indices of the emissions that should be used for each state, if <code>null</code> state <code>i</code> will use emission <code>i</code>
	 * @param forward a boolean array that indicates whether the symbol on the forward or the reverse complementary strand should be used,
	 * 				  if <code>null</code> all states use the forward strand
	 * @param emission the emissions
	 * @param te the {@link AbstractTransitionElement}s building a transition
	 * 
	 * @throws Exception if 
	 * 	<ul>
	 *  <li>some component could not be cloned</li> 
	 *  <li>some the length of <code>name, emissionIdx,</code> or <code>forward</code> is not equal to the number of states</li>
	 *  <li>not all emissions use the same {@link de.jstacs.data.AlphabetContainer}</li>
	 *  <li>the states can not be handled by the transition
	 *  </ul>
	 */
	public HigherOrderHMM( int threads, HMMTrainingParameterSet trainingParameterSet, String[] name, int[] emissionIdx, boolean[] forward, Emission[] emission, AbstractTransitionElement... te ) throws Exception {
		super( threads, trainingParameterSet, name, emissionIdx, forward, emission );
		createStates();
		initTransition( te );
		determineFinalStates();
	}
	
	protected void createHelperVariables(int thread) {
		if(container == null){
			container = new int[threads][];
			logEmission = new double[threads][];
			forwardIntermediate = new double[threads][][][];
			backwardIntermediate = new double[threads][];
			numberOfSummands = new int[threads][][];
			stateList = new IntList[threads];
		}
		if( container[thread] == null ) {
			container[thread] = new int[3];
			
			logEmission[thread] = new double[states[thread].length];
			int m = 0, max = transition[thread].getMaximalMarkovOrder();
			for( int i = 0; i <= max; i++ ) {
				m = Math.max( m, transition[thread].getNumberOfIndexes( i ) );
			}
			forwardIntermediate[thread] = new double[2][m][transition[thread].getMaximalInDegree()+1];
			backwardIntermediate[thread] = new double[states[thread].length+1];
			numberOfSummands[thread] = new int[2][forwardIntermediate[thread][0].length];
			stateList[thread] = new IntList();
		}
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link HigherOrderHMM} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link HigherOrderHMM} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public HigherOrderHMM( StringBuffer xml ) throws NonParsableException {
		super( xml );
		for(int i=0;i<threads;i++){
			createHelperVariables(i);
		}
	}
	
	private static final String XML_TAG = "HigherOrderHMM";

	protected String getXMLTag() {
		return XML_TAG;
	}
		
	protected void appendFurtherInformation( StringBuffer xml ) {
		XMLParser.appendObjectWithTags( xml, skipInit, "skipInit" );
	}

	/**
	 * This method extracts further information from the XML representation. It allows subclasses to cast further parameters that are not defined in the superclass.
	 * 
	 * @param xml the XML representation
	 *  
	 * @throws NonParsableException if the information could not be reconstructed out of the {@link StringBuffer} <code>xml</code>
	 */
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		try{
			skipInit = XMLParser.extractObjectForTags( xml, "skipInit",boolean.class );
		}catch(NonParsableException ex){
			skipInit = false;//TODO
		}
		System.gc();
	}	
	
	public HigherOrderHMM clone() throws CloneNotSupportedException {
		HigherOrderHMM clone = (HigherOrderHMM) super.clone();
		clone.container = null;
		for(int i=0;i<threads;i++){
			clone.createHelperVariables(i);
		}
		return clone;
	}
	
	protected void createStates() {
		this.states = new State[threads][];
		for(int j=0;j<threads;j++){
			this.states[j] = new SimpleState[emissionIdx.length];
			for( int i = 0; i < emissionIdx.length; i++ ) {
				this.states[j][i] = new SimpleState( emission[j][emissionIdx[i]], name[i], forward[i] );
			}
		}
	}
	
	public double getLogPriorTerm() {
		double res = transition[0].getLogPriorTerm();
		for( int e = 0; e < emission[0].length; e++ ) {
			res += emission[0][e].getLogPriorTerm();
		}
		return res;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.models.hmm.AbstractHMM#getLogProbForPath(int[], int, de.jstacs.data.Sequence[])
	 */
	public double getLogProbForPath( IntList path, int startPos, Sequence seq ) throws Exception {
		if( !finalState[path.get( path.length()-1 )] ) {
			throw new IllegalArgumentException( "The last state of the path is no final state. Hence the path is not valid." );
		}
		double res = 0;
		int l = 0, layer = 0, state;
		container[0][1] = 0;
		while( l < path.length() ) {
			state = path.get( l );
			
			int childIdx = transition[0].getChildIdx( layer, container[0][1], state );
			if( childIdx < 0 ) {
				throw new IllegalArgumentException( "Impossible path" );
			}
			res += transition[0].getLogScoreFor( layer, container[0][1], childIdx, seq, startPos ) //transition
				+ states[0][state].getLogScoreFor( startPos, startPos, seq ); //emission
			transition[0].fillTransitionInformation( layer, container[0][1], childIdx, container[0] );
			if( container[0][2] == 1 ) {
				startPos++;
				layer++;
			}
			l++;
		}
		return res;
	}
	
	protected void fillLogStatePosteriorMatrix( double[][] statePosterior, int startPos, int endPos, Sequence seq, boolean silentZero ) throws Exception {
		int len = endPos-startPos+1;
		if( transition[0].getMaximalMarkovOrder() == 0 ) {
			for( int l = 0; l < len; l++ ) {
				for( int s = 0; s < states[0].length; s++ ) {
					transition[0].fillTransitionInformation( 0, 0, s, container[0] );
					logEmission[0][s] = statePosterior[container[0][0]][l+1] = 
						transition[0].getLogScoreFor( 0, 0, s, seq, startPos+l )
						+ states[0][container[0][0]].getLogScoreFor( startPos+l, startPos+l, seq );
				}
				double d = Normalisation.getLogSum( logEmission[0] );
				for( int s = 0; s < states[0].length; s++ ) {
					statePosterior[s][l+1] -= d;
				}
			}
		} else {
			fillFwdMatrix( 0, startPos, endPos, seq );
			fillBwdMatrix( 0, startPos, endPos, seq );
			double logProb = bwdMatrix[0][0][0];
			for( int l = 0; l <= len; l++ ) {
				for( int s = 0; s < states[0].length; s++ ) {
					statePosterior[s][l] = Double.NEGATIVE_INFINITY;
				}
				for( int c = 0; c < fwdMatrix[0][l].length; c++ ) {
					int state = transition[0].getLastContextState( l, c );
					if( state >= 0 && !(silentZero && states[0][state].isSilent()) ) {
						statePosterior[state][l] = Normalisation.getLogSum( statePosterior[state][l], fwdMatrix[0][l][c] + bwdMatrix[0][l][c] );
					}
				}
				for( int s = 0; s < states[0].length; s++ ) {
					statePosterior[s][l] -= logProb;
				}
				startPos++;
			}
		}
	}

	@Override
	protected void fillFwdMatrix( int thread, int startPos, int endPos, Sequence seq ) throws OperationNotSupportedException, WrongLengthException {
		int l = 0, stateID, context, n, h, hh;
		double logTransition;
		provideMatrix( 0, endPos-startPos+1,thread );
		
		//init
		Arrays.fill( numberOfSummands[thread][0], 0 );
		numberOfSummands[thread][0][0] = 1;
		forwardIntermediate[thread][0][0][0] = 0;
		
		//iterate
		while( startPos <= endPos ) {
			for( stateID = 0; stateID < states[thread].length; stateID++ ) {
				logEmission[thread][stateID] = states[thread][stateID].getLogScoreFor(startPos, startPos, seq); 
			}
			
			h=l%2;
			Arrays.fill( numberOfSummands[thread][1-h], 0 );
			
			for( context = 0; context < fwdMatrix[thread][l].length; context++ ) {
				n = transition[thread].getNumberOfChildren(l, context);
	
				if( numberOfSummands[thread][h][context] > 0 ) {
					fwdMatrix[thread][l][context] = Normalisation.getLogSum( 0, numberOfSummands[thread][h][context], forwardIntermediate[thread][h][context] );
				
					for( stateID = 0; stateID < n; stateID++ ) {
						transition[thread].fillTransitionInformation( l, context, stateID, container[thread] );			
						logTransition = transition[thread].getLogScoreFor( l, context, stateID, seq, startPos );

						hh = (h + container[thread][2]) % 2;						
						forwardIntermediate[thread][hh][container[thread][1]][numberOfSummands[thread][hh][container[thread][1]]] = fwdMatrix[thread][l][context] //old part
						       + logEmission[thread][container[thread][0]] //emission
						       + logTransition; //transition
						
						numberOfSummands[thread][hh][container[thread][1]]++;
					}
				} else {
					fwdMatrix[thread][l][context] = Double.NEGATIVE_INFINITY;
				}
			}
			// System.out.println( (l==0?" ":seq.toString(startPos-1, startPos)) + "\t" + l + "\t" + Arrays.toString( fwdMatrix[l] ) );
			l++;
			startPos++;
		}
		
		//final summing and silent states
		h=l%2;
		for( context = 0; context < fwdMatrix[thread][l].length; context++ ) {
			n = transition[thread].getNumberOfChildren(l, context);
			
			if( numberOfSummands[thread][h][context] > 0 ) {
				fwdMatrix[thread][l][context] = Normalisation.getLogSum( 0, numberOfSummands[thread][h][context], forwardIntermediate[thread][h][context] ); 
			
				for( stateID = 0; stateID < n; stateID++ ) {
					transition[thread].fillTransitionInformation( l, context, stateID, container[thread] );
					if( states[thread][container[thread][0]].isSilent() ) {
						logTransition = transition[thread].getLogScoreFor( l, context, stateID, seq, startPos );
	
						//hh=h
						forwardIntermediate[thread][h][container[thread][1]][numberOfSummands[thread][h][container[thread][1]]++] = fwdMatrix[thread][l][context] //old part
						       // there is no emission (silent state)
						       + logTransition; //transition
					}
				} 
			} else {
				fwdMatrix[thread][l][context] = Double.NEGATIVE_INFINITY;
			}
		}
		// System.out.println( (l==0?" ":seq.toString(startPos-1, startPos)) + "\t" + l + "\t" + Arrays.toString( fwdMatrix[l] ) );
	}

	@Override
	protected void fillBwdMatrix( int thread, int startPos, int endPos, Sequence seq )  throws Exception {
		fillBwdOrViterbiMatrix( Type.LIKELIHOOD, thread, startPos, endPos, 1, seq );
	}
	
	/**
	 * This method computes the entries of the backward or the viterbi matrix.
	 * Additionally it allows to modify the sufficient statistics as needed for Baum-Welch training.
	 * 
	 * @param t a switch to decide which computation mode 
	 * @param startPos start position of the sequence 
	 * @param endPos end position of the sequence
	 * @param weight the given external weight of the sequence (only used for Baum-Welch)
	 * @param seq the sequence
	 * 
	 * @throws Exception forwarded from {@link TrainableState#addToStatistic} and {@link de.jstacs.models.hmm.State#getLogScoreFor(int, int, Sequence)}
	 */
	protected void fillBwdOrViterbiMatrix( Type t, int thread, int startPos, int endPos, double weight, Sequence seq )  throws Exception {
		int l = endPos-startPos+1, stateID, context, n;
		boolean zero = transition[thread].getMaximalMarkovOrder() == 0;
		provideMatrix( 1, endPos-startPos+1, thread );

		double val, newWeight, res = t != Type.BAUM_WELCH ? Double.NaN : computeLogScoreFromForward( thread, l );
		
		//init
		for( context = bwdMatrix[thread][l].length-1; context >= 0; context-- ) {
			n = transition[thread].getNumberOfChildren( l, context );			
			numberOfSummands[thread][0][0] = 0;
			
			if( zero || finalState[transition[thread].getLastContextState( l, context )] ) {
				val = 0;
			} else {
				val = Double.NEGATIVE_INFINITY;
			}
			
			//for all different children states
			for( stateID = 0; stateID < n; stateID++ ) {
				transition[thread].fillTransitionInformation( l, context, stateID, container[thread] );
				if( states[thread][container[thread][0]].isSilent() ) {
					backwardIntermediate[thread][numberOfSummands[thread][0][0]] =
						bwdMatrix[thread][l][container[thread][1]] //backward score until next position
						//there is no emission (silent state)
					    + transition[thread].getLogScoreFor( l, context, stateID, seq, endPos ); //transition
					
					if( t == Type.BAUM_WELCH ) {
						newWeight = weight * Math.exp( fwdMatrix[thread][l][context] + backwardIntermediate[thread][numberOfSummands[thread][0][0]] - res );
						
						((TrainableTransition)transition[thread]).addToStatistic( l, context, stateID, newWeight, seq, endPos );
					}
					
					numberOfSummands[thread][0][0]++;
				}
			}
			
			if( numberOfSummands[thread][0][0] == 0 ) {
				bwdMatrix[thread][l][context] = val;
			} else {
				bwdMatrix[thread][l][context] = t==Type.VITERBI
						? Math.max( val, ToolBox.max( 0, numberOfSummands[thread][0][0], backwardIntermediate[thread] ) )
						: Normalisation.getLogSum( val, Normalisation.getLogSum( 0, numberOfSummands[thread][0][0], backwardIntermediate[thread] ) );
			}
		}
		//System.out.println( seq.toString(endPos, endPos+1) + "\t" + l + "\t" + Arrays.toString( bwdMatrix[l] ) );
		
		//compute scores for all positions backward
		while( --l >= 0 ) {
			for( stateID = 0; stateID < states[thread].length; stateID++ ) {
				logEmission[thread][stateID] = states[thread][stateID].getLogScoreFor(endPos, endPos, seq); 
			}
			//for all different contexts
			for( context = bwdMatrix[thread][l].length-1; context >= 0; context-- ) {
				n = transition[thread].getNumberOfChildren( l, context );			
				//for all different children states
				for( stateID = 0; stateID < n; stateID++ ) {
					transition[thread].fillTransitionInformation( l, context, stateID, container[thread] );
					
					backwardIntermediate[thread][stateID] =
						bwdMatrix[thread][l+container[thread][2]][container[thread][1]] //backward score until next position
						+ logEmission[thread][container[thread][0]] //emission
					    + transition[thread].getLogScoreFor( l, context, stateID, seq, endPos ); //transition
					
					if( t == Type.BAUM_WELCH ) {
						newWeight = weight * Math.exp( fwdMatrix[thread][l][context] + backwardIntermediate[thread][stateID] - res );
						((TrainableState)states[thread][container[thread][0]]).addToStatistic( endPos, endPos, newWeight, seq );
						((TrainableTransition)transition[thread]).addToStatistic( l, context, stateID, newWeight, seq, endPos );
					}
				}
				if( n > 0 ) {
					bwdMatrix[thread][l][context] = t == Type.VITERBI
						? ToolBox.max( 0, n, backwardIntermediate[thread] )
						: Normalisation.getLogSum( 0, n, backwardIntermediate[thread] );
				} else {
					bwdMatrix[thread][l][context] = Double.NEGATIVE_INFINITY;
				}
			}
			endPos--;

			//System.out.println( (l==0?" ":seq.toString(endPos, endPos+1)) + "\t" + l + "\t" + Arrays.toString( bwdMatrix[l] ) );
		}
	}
	
	/**
	 * This enum defined different types of computations that will be done using the backward algorithm.
	 * 
	 * @author Jens Keilwagen
	 */
	protected static enum Type {
		/**
		 * Defines to compute the likelihood: {@latex.inline $P(\\underline{x}|\\underline{\\lambda}) = \\sum_{\\underline{u}} P(\\underline{x},\\underline{u}|\\underline{\\lambda})$}.
		 */
		LIKELIHOOD,
		/**
		 * Defines to compute the viterbi score: {@latex.inline $\\max_{\\underline{u}} P(\\underline{x},\\underline{u}|\\underline{\\lambda})$}.
		 */
		VITERBI,
		/**
		 * Defines to do Baum-Welch for the given sequence. In this case the forward matrix has to be filled.
		 */
		BAUM_WELCH;
	}
		
	public Pair<IntList,Double> getViterbiPathFor(int startPos, int endPos, Sequence seq ) throws Exception {
		IntList path = new IntList(endPos-startPos+1);
		double score = viterbi( 0, path, startPos, endPos, 0, seq );
		return new Pair<IntList, Double>( path, score );
	}
	
	
	/**
	 * This method computes the viterbi score of a given sequence <code>seq</code>.
	 * Furthermore, it allows either to modify the sufficient statistics according
	 * to the viterbi training algorithm or to compute the viterbi path, which will
	 * in this case be returned in <code>path</code>.
	 * 
	 * @param path if <code>null</code> viterbi training, otherwise computation of the viterbi path
	 * @param startPos the start position
	 * @param endPos the end position
	 * @param weight the sequence weight, in most cases this is 1
	 * @param seq the sequence
	 * 
	 * @return the viterbi score of the sequence
	 * 
	 * @throws Exception an error occurs during the computation
	 */
	protected double viterbi( int thread, IntList path, int startPos, int endPos, double weight, Sequence seq ) throws Exception {
		fillBwdOrViterbiMatrix( Type.VITERBI, thread, startPos, endPos, 0, seq );
		int l = endPos-startPos+1, n, layer = 0, stateID, context = 0, add, state, newContext, childIdx;
		double current, dist, bestDist;

		if( path != null ) {
			path.clear();
		}
			
		//fill
		while( layer < l ) {
			n = transition[thread].getNumberOfChildren( layer, context );
			
			bestDist = Double.POSITIVE_INFINITY;
			childIdx = state = newContext = add = -1000;			
			for( stateID = 0; stateID < n; stateID++ ) {
				transition[thread].fillTransitionInformation( layer, context, stateID, container[thread] );
				
				current =
					bwdMatrix[thread][layer+container[thread][2]][container[thread][1]] //score until next position
					+ states[thread][container[thread][0]].getLogScoreFor( startPos, startPos, seq ) //emission
				    + transition[thread].getLogScoreFor( layer, context, stateID, seq, startPos ); //transition
				
				dist = current - bwdMatrix[thread][layer][context];
				
				dist*=dist;
				if( dist < bestDist ) {
					childIdx = stateID;
					state = container[thread][0];
					newContext = container[thread][1];
					add = container[thread][2];
					bestDist = dist;
				}
			}
			
			//System.out.println( layer + "\t" + state + "\t" + bwdMatrix[layer][context] + "\t" + n + "\t" + bestDist );
			
			if( path == null ) {
				((TrainableTransition)transition[thread]).addToStatistic( layer, context, childIdx, weight, seq, startPos );
				((TrainableState)states[thread][state]).addToStatistic( startPos, startPos, weight, seq );
			} else {
				path.add( state );
			}

			startPos += add;
			layer += add;
			context = newContext;
		}
		
		//add silent sates at the end
		do {
			n = transition[thread].getNumberOfChildren( layer, context );
			
			dist = finalState[transition[thread].getLastContextState(layer, context)] ? 0 - bwdMatrix[thread][layer][context] : Double.NEGATIVE_INFINITY;
			bestDist = dist*dist;
			childIdx = state = newContext = add = -1000;
			for( stateID = 0; stateID < n; stateID++ ) {
				transition[thread].fillTransitionInformation( layer, context, stateID, container[thread] );
				
				if( container[thread][2] == 0 ) {
					current =
						bwdMatrix[thread][layer][container[thread][1]] //score until next position
						+ states[thread][container[thread][0]].getLogScoreFor( startPos, startPos, seq ) //emission
					    + transition[thread].getLogScoreFor( layer, context, stateID, seq, startPos ); //transition
					
					dist = current - bwdMatrix[thread][layer][context];
					
					dist*=dist;
					if( dist < bestDist ) {
						childIdx = stateID;
						state = container[thread][0];
						newContext = container[thread][1];
						bestDist = dist;
					}
				}
			}
			
			if( state >= 0 ) {
				if( path == null ) {
					((TrainableTransition)transition[thread]).addToStatistic( layer, context, childIdx, weight, seq, startPos );
					((TrainableState)states[thread][state]).addToStatistic( startPos, startPos, weight, seq );
				} else {
					path.add( state );
				}
				context = newContext;
			} else {
				break;
			}
		} while ( true );
		
		return bwdMatrix[thread][0][0];
	}
	
	/**
	 * This method computes the likelihood and modifies the sufficient statistics according to the Baum-Welch algorithm.
	 * 
	 * @param startPos the start position
	 * @param endPos the end position
	 * @param weight the sequence weight, in most cases this is 1
	 * @param seq the sequence
	 * 
	 * @return the likelihood of the sequence
	 * 
	 * @throws Exception an error occurs during the computation
	 */
	protected double baumWelch( int thread, int startPos, int endPos, double weight, Sequence seq ) throws Exception {
		fillFwdMatrix( thread, startPos, endPos, seq );
		fillBwdOrViterbiMatrix( Type.BAUM_WELCH, thread, startPos, endPos, weight, seq );
		return bwdMatrix[thread][0][0];
	}
	
	public void setSkipInit(boolean skipInit){
		this.skipInit = skipInit;
	}
	
	public synchronized void train(Sample data, double[] weights) throws Exception {
		if( !(trainingParameter instanceof MaxHMMTrainingParameterSet) ) {
			throw new IllegalArgumentException( "This kind of training is currently not supported." );
		} else {
			Transition bestTransition = null;
			Emission[] bestEmissions = null;
			double best = Double.NEGATIVE_INFINITY;
			
			int N = data.getNumberOfElements();
			
			int numberOfStarts = trainingParameter.getNumberOfStarts();
			AbstractTerminationCondition tc = ((MaxHMMTrainingParameterSet) trainingParameter).getTerminantionCondition();
			WorkerThread[] workers = new WorkerThread[threads];
			int last = 0;
			for(int i=0;i<threads-1;i++){
				workers[i] = new WorkerThread( i, last, (i+1)*N/workers.length, data, weights );
				workers[i].start();
				last = workers[i].end;
			}
			workers[ workers.length-1 ] = new WorkerThread( workers.length-1, last, N, data, weights );
			workers[ workers.length-1 ].start();
			
			RealTime time = new RealTime();
			for( int it, start = 0; start < numberOfStarts; start++ ) {
				sostream.writeln( "start " + start + " ============================" );
				//init
				if(!skipInit){
					initialize( data, weights );
				}
				
				//iterate
				double old_value, new_value = Double.NEGATIVE_INFINITY;
				it = 0;
				time.reset();
				do {
					resetStatistics();
					old_value = new_value;
					new_value = getLogPriorTerm();
					
					for(int i=0;i<workers.length;i++){
						workers[i].setState( WorkerState.COMPUTE );
					}
					waitUntilWorkersFinished(workers);
					for(int i=0;i<workers.length;i++){
						new_value += workers[i].getScore();
					}
					
					sostream.writeln( it++ + "\t" + time.getElapsedTime() + "\t" + new_value + "\t" + (new_value - old_value) );
					if( tc.doNextIteration( it, old_value, new_value, null, null, Double.NaN, time) ) {
						estimateFromStatistics();
					} else {
						break;
					}					
				}while( true );
				
				//check for best
				if( new_value > best ) {
					best = new_value;
					if( numberOfStarts > 1 ) {
						bestEmissions = ArrayHandler.clone( emission[0] );
						bestTransition = transition[0].clone();
					}
				}
			}
			sostream.writeln( "best result: " + best );
			if( bestEmissions != null ) {
				emission = ArrayHandler.createArrayOf( bestEmissions,emission.length);
				createStates();
				transition = ArrayHandler.createArrayOf( bestTransition, transition.length );
			}
			for(int i=0;i<workers.length;i++){
				workers[i].setState( WorkerState.STOP );
			}
		}
	}
	
	private void waitUntilWorkersFinished(WorkerThread[] workers){
		int i;
		while( true )
		{
			i = 0;
			while( i < workers.length && workers[i].isWaiting() ){
				i++;
			}
			if( i == workers.length ){
				break;
			}else{
				try{
					wait();
				} catch( InterruptedException e ) { }
			}
		}
	}
	
	private double doOneStep(Sample data, double[] weights, int start, int end, int thread) throws Exception{
		Sequence seq;
		double weight = 1;
		double newValue = 0;
		for( int n = start; n < end; n++ ) {
			seq = data.getElementAt( n );
			if( weights != null ) {
				weight = weights[n];
			}
			
			if( trainingParameter instanceof ViterbiParameterSet ) {
				newValue += viterbi( thread, null, 0, seq.getLength()-1, weight, seq ); //viterbi
			} else if ( trainingParameter instanceof BaumWelchParameterSet ) {
				newValue += baumWelch( thread, 0, seq.getLength()-1, weight, seq ); //Baum-Welch
			} else {
				throw new IllegalArgumentException( "Training mode not available." );
			}
		}
		return newValue;
	}
	
	/**
	 * This method initializes all emissions and the transition.
	 * The initialization might use the data, but the default
	 * implementation refers to {@link #initializeRandomly()}.
	 * 
	 * @throws Exception if an error occurs during the initialization 
	 */
	protected void initialize( Sample data, double[] weight ) throws Exception {
		initializeRandomly();
	}
	
	/**
	 * This method initializes all emissions and the transition randomly.
	 */
	protected void initializeRandomly() {
		transition[0].initializeRandomly();
		for( int e = 0; e < emission[0].length; e++ ) {
			emission[0][e].initializeFunctionRandomly();
		}
		prepare();
	}
	
	protected void prepare() {
		for(int i=1;i<threads;i++){
			try{
				transition[i] = transition[0].clone();
				emission[i] = ArrayHandler.clone( emission[0] );
			}catch(CloneNotSupportedException e){
				throw new RuntimeException( e );
			}
		}
		createStates();
	}
	
	/**
	 * This method resets all sufficient statistics of all emissions and the transition.
	 */
	protected void resetStatistics() {
		for(int i=0;i<transition.length;i++){
			((TrainableTransition) transition[i]).resetStatistic();
			for( int e = 0; e < emission[i].length; e++ ) {
				emission[i][e].resetStatistic();
			}
		}
	}
	
	/**
	 * This method estimates the parameters of all emissions and the transition using their sufficient statistics.
	 */
	protected void estimateFromStatistics() {
		((TrainableTransition)transition[0]).joinStatistics( transition );
		Emission[] temp = new Emission[threads];
		for(int e = 0; e < emission[0].length; e++){
			for(int j=0;j<threads;j++){
				temp[j] = emission[j][e];
			}
			emission[0][e].joinStatistics( temp );
		}
		for(int i=0;i<threads;i++){
			((TrainableTransition) transition[i]).estimateFromStatistic();
			for( int e = 0; e < emission[i].length; e++ ) {
				emission[i][e].estimateFromStatistic();
			}
		}
	}

	public final byte getMaximalMarkovOrder() throws UnsupportedOperationException {
		// TODO Auto-generated method stub
		return Byte.MAX_VALUE;
	}
	
	public ResultSet getCharacteristics() throws Exception {
		return new ResultSet( getNumericalCharacteristics().getResults(),
				new Result[] { new StorableResult("model", "the xml representation of the model", this) } );
	}
	
	public String getInstanceName() {
		return "HMM(" + transition[0].getMaximalMarkovOrder() + ") " + trainingParameter.getClass().getSimpleName();
	}
	
	public double[] getLogProbFor(Sample data) throws Exception {
		double[] logProb = new double[data.getNumberOfElements()];
		getLogProbFor(data, logProb);
		return logProb;
	}

	public void getLogProbFor(Sample data, double[] res) throws Exception {
		if( !data.getAlphabetContainer().checkConsistency(getAlphabetContainer()) ) {
			throw new WrongAlphabetException( "The AlphabetContainer of the sample and the model do not match." );
		}
		int len = getLength(), l = data.getElementLength();
		if( len != 0 && l != len ) {
			throw new WrongLengthException( "The length of the sample and the model do not match." );
		}
		Sequence seq;
		for( int n = 0; n < data.getNumberOfElements(); n++ ) {
			seq = data.getElementAt(n);
			res[n] = logProb(0,seq.getLength()-1,seq);
		}
	}
	
	public Sample emitSample(int numberOfSequences, int... seqLength)
		throws NotTrainedException, Exception {
		// TODO Auto-generated method stub
		return null;
	}
	
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean isInitialized() {
		return true;
	}
	
	protected void finalize() throws Throwable {
		emission = null;
		emissionIdx = null;
		forward = null;
		name = null;
		container = null;
		numberOfSummands = null;
		logEmission = null;
		forwardIntermediate = null;
		backwardIntermediate = null;
		super.finalize();
	}
	
	/**
	 * This method samples a valid path for the given sequence <code>seq</code> using the internal parameters.
	 *  
	 * @param path an {@link IntList} containing the path after using this method 
	 * @param startPos the start position
	 * @param endPos the end position
	 * @param seq the sequence
	 * 
	 * @throws Exception if an error occurs during computation
	 */
	public void samplePath( int thread, IntList path, int startPos, int endPos, Sequence seq ) throws Exception {		
		fillBwdMatrix( thread, startPos, endPos, seq );
		
		int l = 0, stateID, context = 0, n;
		double logTransition;
		provideMatrix( 0, endPos-startPos+1,thread );
		path.clear();
		
		//iterate
		while( startPos <= endPos ) {
			
			//compute
			n = transition[thread].getNumberOfChildren( l, context );
			for( stateID = 0; stateID < n; stateID++ ) {
				transition[thread].fillTransitionInformation( l, context, stateID, container[thread] );			
				logTransition = transition[thread].getLogScoreFor( l, context, stateID, seq, startPos );
						
				backwardIntermediate[thread][stateID] = bwdMatrix[thread][l+container[thread][2]][container[thread][1]] //old part
				       + states[thread][container[thread][0]].getLogScoreFor( startPos, startPos, seq ) //emission
				       + logTransition; //transition
			}
			Normalisation.logSumNormalisation( backwardIntermediate[thread], 0, n );
			
			//draw
			stateID = AbstractMixtureModel.draw( backwardIntermediate[thread], 0 );
			transition[thread].fillTransitionInformation( l, context, stateID, container[thread] );
			path.add( container[thread][0] );
			context = container[thread][1];
			l += container[thread][2];
			startPos += container[thread][2];
		}
		
		//final silent states
		int add, last = path.get( path.length()-1 );
		do {
			//compute
			n = transition[thread].getNumberOfChildren( l, context );
			stateList[thread].clear();
			for( stateID = 0; stateID < n; stateID++ ) {
				transition[thread].fillTransitionInformation( l, context, stateID, container[thread] );
				if( states[thread][container[thread][0]].isSilent() ) {
					logTransition = transition[thread].getLogScoreFor( l, context, stateID, seq, startPos );
					backwardIntermediate[thread][stateList[thread].length()] = bwdMatrix[thread][l][container[thread][1]] //old part
							       // there is no emission (silent state)
							       + logTransition; //transition
					stateList[thread].add( stateID );
				}
			}
			if( finalState[last] ) {
				backwardIntermediate[thread][stateList[thread].length()]=0;
				add = 1;
			} else {
				add = 0;
			}
			Normalisation.logSumNormalisation( backwardIntermediate[thread], 0, stateList[thread].length()+add );
			
			//draw
			n = AbstractMixtureModel.draw( backwardIntermediate[thread], 0 );
			if( add == 1 && n == stateList[thread].length() ) {
				break;
			} else {
				transition[thread].fillTransitionInformation( l, context, stateList[thread].get(n), container[thread] );
				path.add( container[thread][0] );
				context = container[thread][1];
				l += container[thread][2];
				startPos += container[thread][2];
				
				last = container[thread][0];
			}
		} while( true );
	}
		
	private double computeLogScoreFromForward( int thread, int l ) {
		double res = Double.NEGATIVE_INFINITY;
		if( transition[thread].getMaximalMarkovOrder() > 0 ) {
			for( int i = 0; i < fwdMatrix[thread][l].length; i++ ) {
				if( finalState[transition[thread].getLastContextState(l, i)] ) {
					res = Normalisation.getLogSum( res, fwdMatrix[thread][l][i] );
				}
			}
		} else {
			res = Normalisation.getLogSum( fwdMatrix[thread][l] );
		}
		return res;
	}
		
	/*
	public void checkLogProbFor( Sequence sequence, int startpos, int endpos ) //currently only for HMMs with silent states
			throws Exception {
		int l = endpos - startpos+1, len = getLength();
		if( !sequence.getAlphabetContainer().checkConsistency(getAlphabetContainer()) ) {
			throw new WrongAlphabetException( "The AlphabetContainer of the sequence and the model do not match." );
		}
		if( len != 0 && l != len ) {
			throw new WrongLengthException( "The given start position ("+ startpos + ") and end position (" + endpos + ") yield an length of " + (l+1) + " which is not possible for the current model that models sequences of length " + len + "." );
		}
		l++;
		
		fillBwdMatrix( startpos, endpos, sequence );
		fillFwdMatrix( startpos, endpos, sequence );
		for( len = 0; len < l; len++ ) {
			double res = Double.NEGATIVE_INFINITY;
			for( int c = 0; c < fwdMatrix[len].length; c++ ) {
				res = Normalisation.getLogSum( res, fwdMatrix[len][c] + bwdMatrix[len][c] );
			}
			System.out.println( len + "\t" + res );
		}
	}
	
	/*protected double logProb( int startpos, int endpos, Sequence sequence ) throws Exception {
		try {
			fillBwdMatrix(startpos, endpos, sequence);
		} catch( Exception e ) {
			throw getRunTimeException( e );
		}
		
		System.out.println("-------------");
		fillFwdMatrix(startpos, endpos, sequence);
		double f = computeLogScoreFromForward(endpos-startpos+1);
		System.out.println( f + "\t" + bwdMatrix[0][0] );
		return bwdMatrix[0][0];
	}
	/**/
	
	private enum WorkerState{
		COMPUTE,
		WAIT,
		STOP
	}
	
	private class WorkerThread extends Thread{

		private WorkerState state;
		private int idx;
		private int start, end;
		private double score;
		private Sample data;
		private double[] weights;

		
		
		public WorkerThread(int idx, int start, int end, Sample data, double[] weights){
			this.idx = idx;
			this.start = start;
			this.end = end;
			this.score = 0;
			this.state = WorkerState.WAIT;
			this.data = data;
			this.weights = weights;
			this.setDaemon( true );
		}
		
		public double getScore() {
			return score;
		}

		public synchronized void setState(WorkerState state){
			this.state = state;
			notify();
		}
		
		@Override
		public synchronized void run() {
			//System.out.println(idx+" : "+state);
			while(state != WorkerState.STOP){
				//System.out.println(idx+" ; "+state);
				if(state == WorkerState.WAIT){
					try {
						wait();
					} catch ( InterruptedException e ) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}else{
					try{
						if(state == WorkerState.COMPUTE){
							score = doOneStep( data, weights, start, end, idx );
						}
						
						synchronized( HigherOrderHMM.this )
						{
							state = WorkerState.WAIT;
							HigherOrderHMM.this.notify();
						}
					}catch( Exception e ){
						RuntimeException re = new RuntimeException( e.getClass().getName() + ": " + e.getMessage() );
						re.setStackTrace( e.getStackTrace() );
						throw re;
					}
				}
			}
		}
		
		public boolean isWaiting(){
			return state == WorkerState.WAIT;
		}
		
	}
	
	
}
