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

package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models;

import java.util.Arrays;

import javax.naming.OperationNotSupportedException;

import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.SimpleState;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.TrainableState;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.BaumWelchParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.HMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.MaxHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.ViterbiParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.TrainableTransition;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.Transition;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.BasicHigherOrderTransition.AbstractTransitionElement;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.Time;
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
	protected int[] container;
	/**
	 * Helper variable = only for internal use. This array is used for compute the emissions at each position of a sequence only once,
	 * which might be beneficial for higher order models.
	 * 
	 * @see HigherOrderHMM#emission
	 */
	protected double[] logEmission;
	
	/**
	 * Helper variable = only for internal use. This array is used to compute the forward matrix. It stores intermediate results.
	 * 
	 * #see {@link #numberOfSummands}
	 */
	private double[][][] forwardIntermediate;

	/**
	 * Helper variable = only for internal use. This array is used to compute the backward matrix. It stores intermediate results.
	 * 
	 * #see {@link #numberOfSummands}
	 */
	protected double[] backwardIntermediate;
	
	/**
	 * Helper variable = only for internal use. This array is used to compute the forward and backward matrix. It stores the number of intermediate results.
	 */
	protected int[][] numberOfSummands;
	
	/**
	 * Helper variable = only for internal use. This field is used in the method {@link #samplePath(IntList, int, int, Sequence)}.
	 */
	protected IntList stateList;
	
	/**
	 * Indicates if the model should be initialized (randomly) before optimization
	 */
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
	 * @see HigherOrderHMM#HigherOrderHMM(HMMTrainingParameterSet, String[], int[], boolean[], Emission[], de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.BasicHigherOrderTransition.AbstractTransitionElement...)
	 */
	public HigherOrderHMM( HMMTrainingParameterSet trainingParameterSet, String[] name, Emission[] emission, AbstractTransitionElement... te ) throws Exception {
		this( trainingParameterSet, name, null, null, emission, te );
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
	public HigherOrderHMM( HMMTrainingParameterSet trainingParameterSet, String[] name, int[] emissionIdx, boolean[] forward, Emission[] emission, AbstractTransitionElement... te ) throws Exception {
		super( trainingParameterSet, name, emissionIdx, forward, emission );
		createStates();
		initTransition( te );
		determineFinalStates();
	}
	
	protected void createHelperVariables() {
		if( container == null ) {
			container = new int[3];
			
			logEmission = new double[states.length];
			int m = 0, max = transition.getMaximalMarkovOrder();
			for( int i = 0; i <= max; i++ ) {
				m = Math.max( m, transition.getNumberOfIndexes( i ) );
			}
			forwardIntermediate = new double[2][m][transition.getMaximalInDegree()+1];
			backwardIntermediate = new double[states.length+1];
			numberOfSummands = new int[2][forwardIntermediate[0].length];
			stateList = new IntList();
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
		createHelperVariables();
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
		skipInit = XMLParser.extractObjectForTags( xml, "skipInit",boolean.class );
	}	
	
	public HigherOrderHMM clone() throws CloneNotSupportedException {
		HigherOrderHMM clone = (HigherOrderHMM) super.clone();
		clone.container = null;
		clone.createHelperVariables();
		return clone;
	}
	
	protected void createStates() {
		this.states = new SimpleState[emissionIdx.length];
		for( int i = 0; i < emissionIdx.length; i++ ) {
			this.states[i] = new SimpleState( emission[emissionIdx[i]], name[i], forward[i] );
		}
	}
	
	public double getLogPriorTerm() {
		double res = transition.getLogPriorTerm();
		for( int e = 0; e < emission.length; e++ ) {
			res += emission[e].getLogPriorTerm();
		}
		return res;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM#getLogProbForPath(int[], int, de.jstacs.data.Sequence[])
	 */
	public double getLogProbForPath( IntList path, int startPos, Sequence seq ) throws Exception {
		if( !finalState[path.get( path.length()-1 )] ) {
			throw new IllegalArgumentException( "The last state of the path is no final state. Hence the path is not valid." );
		}
		double res = 0;
		int l = 0, layer = 0, state;
		container[1] = 0;
		while( l < path.length() ) {
			state = path.get( l );
			
			int childIdx = transition.getChildIdx( layer, container[1], state );
			if( childIdx < 0 ) {
				throw new IllegalArgumentException( "Impossible path" );
			}
			res += transition.getLogScoreFor( layer, container[1], childIdx, seq, startPos ) //transition
				+ states[state].getLogScoreFor( startPos, startPos, seq ); //emission
			transition.fillTransitionInformation( layer, container[1], childIdx, container );
			if( container[2] == 1 ) {
				startPos++;
				layer++;
			}
			l++;
		}
		return res;
	}
	
	protected void fillLogStatePosteriorMatrix( double[][] statePosterior, int startPos, int endPos, Sequence seq, boolean silentZero ) throws Exception {
		int len = endPos-startPos+1;
		if( transition.getMaximalMarkovOrder() == 0 ) {
			for( int l = 0; l < len; l++ ) {
				for( int s = 0; s < states.length; s++ ) {
					transition.fillTransitionInformation( 0, 0, s, container );
					logEmission[s] = statePosterior[container[0]][l+1] = 
						transition.getLogScoreFor( 0, 0, s, seq, startPos+l )
						+ states[container[0]].getLogScoreFor( startPos+l, startPos+l, seq );
				}
				double d = Normalisation.getLogSum( logEmission );
				for( int s = 0; s < states.length; s++ ) {
					statePosterior[s][l+1] -= d;
				}
			}
		} else {
			fillFwdMatrix( startPos, endPos, seq );
			fillBwdMatrix( startPos, endPos, seq );
			double logProb = bwdMatrix[0][0];
			for( int l = 0; l <= len; l++ ) {
				for( int s = 0; s < states.length; s++ ) {
					statePosterior[s][l] = Double.NEGATIVE_INFINITY;
				}
				for( int c = 0; c < fwdMatrix[l].length; c++ ) {
					int state = transition.getLastContextState( l, c );
					if( state >= 0 && !(silentZero && states[state].isSilent()) ) {
						statePosterior[state][l] = Normalisation.getLogSum( statePosterior[state][l], fwdMatrix[l][c] + bwdMatrix[l][c] );
					}
				}
				for( int s = 0; s < states.length; s++ ) {
					statePosterior[s][l] -= logProb;
				}
				startPos++;
			}
		}
	}

	@Override
	protected void fillFwdMatrix( int startPos, int endPos, Sequence seq ) throws OperationNotSupportedException, WrongLengthException {
		int l = 0, stateID, context, n, h, hh;
		double logTransition;
		provideMatrix( 0, endPos-startPos+1 );
		
		//init
		Arrays.fill( numberOfSummands[0], 0 );
		numberOfSummands[0][0] = 1;
		forwardIntermediate[0][0][0] = 0;
		
		//iterate
		while( startPos <= endPos ) {
			for( stateID = 0; stateID < states.length; stateID++ ) {
				logEmission[stateID] = states[stateID].getLogScoreFor(startPos, startPos, seq); 
			}
			
			h=l%2;
			Arrays.fill( numberOfSummands[1-h], 0 );
			
			for( context = 0; context < fwdMatrix[l].length; context++ ) {
				n = transition.getNumberOfChildren(l, context);
	
				if( numberOfSummands[h][context] > 0 ) {
					fwdMatrix[l][context] = Normalisation.getLogSum( 0, numberOfSummands[h][context], forwardIntermediate[h][context] );
				
					for( stateID = 0; stateID < n; stateID++ ) {
						transition.fillTransitionInformation( l, context, stateID, container );			
						logTransition = transition.getLogScoreFor( l, context, stateID, seq, startPos );

						hh = (h + container[2]) % 2;						
						forwardIntermediate[hh][container[1]][numberOfSummands[hh][container[1]]] = fwdMatrix[l][context] //old part
						       + logEmission[container[0]] //emission
						       + logTransition; //transition
						
						numberOfSummands[hh][container[1]]++;
					}
				} else {
					fwdMatrix[l][context] = Double.NEGATIVE_INFINITY;
				}
			}
			// System.out.println( (l==0?" ":seq.toString(startPos-1, startPos)) + "\t" + l + "\t" + Arrays.toString( fwdMatrix[l] ) );
			l++;
			startPos++;
		}
		
		//final summing and silent states
		h=l%2;
		for( context = 0; context < fwdMatrix[l].length; context++ ) {
			n = transition.getNumberOfChildren(l, context);
			
			if( numberOfSummands[h][context] > 0 ) {
				fwdMatrix[l][context] = Normalisation.getLogSum( 0, numberOfSummands[h][context], forwardIntermediate[h][context] ); 
			
				for( stateID = 0; stateID < n; stateID++ ) {
					transition.fillTransitionInformation( l, context, stateID, container );
					if( states[container[0]].isSilent() ) {
						logTransition = transition.getLogScoreFor( l, context, stateID, seq, startPos );
	
						//hh=h
						forwardIntermediate[h][container[1]][numberOfSummands[h][container[1]]++] = fwdMatrix[l][context] //old part
						       // there is no emission (silent state)
						       + logTransition; //transition
					}
				} 
			} else {
				fwdMatrix[l][context] = Double.NEGATIVE_INFINITY;
			}
		}
		// System.out.println( (l==0?" ":seq.toString(startPos-1, startPos)) + "\t" + l + "\t" + Arrays.toString( fwdMatrix[l] ) );
	}

	@Override
	protected void fillBwdMatrix( int startPos, int endPos, Sequence seq )  throws Exception {
		fillBwdOrViterbiMatrix( Type.LIKELIHOOD, startPos, endPos, 1, seq );
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
	 * @throws Exception forwarded from {@link TrainableState#addToStatistic} and {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.State#getLogScoreFor(int, int, Sequence)}
	 */
	protected void fillBwdOrViterbiMatrix( Type t, int startPos, int endPos, double weight, Sequence seq )  throws Exception {
		int l = endPos-startPos+1, stateID, context, n;
		boolean zero = transition.getMaximalMarkovOrder() == 0;
		provideMatrix( 1, endPos-startPos+1 );

		double val, newWeight, res = t != Type.BAUM_WELCH ? Double.NaN : computeLogScoreFromForward( l );
		
		//init
		for( context = bwdMatrix[l].length-1; context >= 0; context-- ) {
			n = transition.getNumberOfChildren( l, context );			
			numberOfSummands[0][0] = 0;
			
			if( zero || finalState[transition.getLastContextState( l, context )] ) {
				val = 0;
			} else {
				val = Double.NEGATIVE_INFINITY;
			}
			
			//for all different children states
			for( stateID = 0; stateID < n; stateID++ ) {
				transition.fillTransitionInformation( l, context, stateID, container );
				if( states[container[0]].isSilent() ) {
					backwardIntermediate[numberOfSummands[0][0]] =
						bwdMatrix[l][container[1]] //backward score until next position
						//there is no emission (silent state)
					    + transition.getLogScoreFor( l, context, stateID, seq, endPos ); //transition
					
					if( t == Type.BAUM_WELCH ) {
						newWeight = weight * Math.exp( fwdMatrix[l][context] + backwardIntermediate[numberOfSummands[0][0]] - res );
						
						((TrainableTransition)transition).addToStatistic( l, context, stateID, newWeight, seq, endPos );
					}
					
					numberOfSummands[0][0]++;
				}
			}
			
			if( numberOfSummands[0][0] == 0 ) {
				bwdMatrix[l][context] = val;
			} else {
				bwdMatrix[l][context] = t==Type.VITERBI
						? Math.max( val, ToolBox.max( 0, numberOfSummands[0][0], backwardIntermediate ) )
						: Normalisation.getLogSum( val, Normalisation.getLogSum( 0, numberOfSummands[0][0], backwardIntermediate ) );
			}
		}
		//System.out.println( seq.toString(endPos, endPos+1) + "\t" + l + "\t" + Arrays.toString( bwdMatrix[l] ) );
		
		//compute scores for all positions backward
		while( --l >= 0 ) {
			for( stateID = 0; stateID < states.length; stateID++ ) {
				logEmission[stateID] = states[stateID].getLogScoreFor(endPos, endPos, seq); 
			}
			//for all different contexts
			for( context = bwdMatrix[l].length-1; context >= 0; context-- ) {
				n = transition.getNumberOfChildren( l, context );			
				//for all different children states
				for( stateID = 0; stateID < n; stateID++ ) {
					transition.fillTransitionInformation( l, context, stateID, container );
					
					backwardIntermediate[stateID] =
						bwdMatrix[l+container[2]][container[1]] //backward score until next position
						+ logEmission[container[0]] //emission
					    + transition.getLogScoreFor( l, context, stateID, seq, endPos ); //transition
					
					if( t == Type.BAUM_WELCH ) {
						newWeight = weight * Math.exp( fwdMatrix[l][context] + backwardIntermediate[stateID] - res );
						((TrainableState)states[container[0]]).addToStatistic( endPos, endPos, newWeight, seq );
						((TrainableTransition)transition).addToStatistic( l, context, stateID, newWeight, seq, endPos );
					}
				}
				if( n > 0 ) {
					bwdMatrix[l][context] = t == Type.VITERBI
						? ToolBox.max( 0, n, backwardIntermediate )
						: Normalisation.getLogSum( 0, n, backwardIntermediate );
				} else {
					bwdMatrix[l][context] = Double.NEGATIVE_INFINITY;
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
		double score = viterbi( path, startPos, endPos, 0, seq );
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
	protected double viterbi( IntList path, int startPos, int endPos, double weight, Sequence seq ) throws Exception {
		fillBwdOrViterbiMatrix( Type.VITERBI, startPos, endPos, 0, seq );
		int l = endPos-startPos+1, n, layer = 0, stateID, context = 0, add, state, newContext, childIdx;
		double current, dist, bestDist;

		if( path != null ) {
			path.clear();
		}
			
		//fill
		while( layer < l ) {
			n = transition.getNumberOfChildren( layer, context );
			
			bestDist = Double.POSITIVE_INFINITY;
			childIdx = state = newContext = add = -1000;			
			for( stateID = 0; stateID < n; stateID++ ) {
				transition.fillTransitionInformation( layer, context, stateID, container );
				
				current =
					bwdMatrix[layer+container[2]][container[1]] //score until next position
					+ states[container[0]].getLogScoreFor( startPos, startPos, seq ) //emission
				    + transition.getLogScoreFor( layer, context, stateID, seq, startPos ); //transition
				
				dist = current - bwdMatrix[layer][context];
				
				dist*=dist;
				if( dist < bestDist ) {
					childIdx = stateID;
					state = container[0];
					newContext = container[1];
					add = container[2];
					bestDist = dist;
				}
			}
			
			//System.out.println( layer + "\t" + state + "\t" + bwdMatrix[layer][context] + "\t" + n + "\t" + bestDist );
			
			if( path == null ) {
				((TrainableTransition)transition).addToStatistic( layer, context, childIdx, weight, seq, startPos );
				((TrainableState)states[state]).addToStatistic( startPos, startPos, weight, seq );
			} else {
				path.add( state );
			}

			startPos += add;
			layer += add;
			context = newContext;
		}
		
		//add silent sates at the end
		do {
			n = transition.getNumberOfChildren( layer, context );
			
			dist = finalState[transition.getLastContextState(layer, context)] ? 0 - bwdMatrix[layer][context] : Double.NEGATIVE_INFINITY;
			bestDist = dist*dist;
			childIdx = state = newContext = add = -1000;
			for( stateID = 0; stateID < n; stateID++ ) {
				transition.fillTransitionInformation( layer, context, stateID, container );
				
				if( container[2] == 0 ) {
					current =
						bwdMatrix[layer][container[1]] //score until next position
						+ states[container[0]].getLogScoreFor( startPos, startPos, seq ) //emission
					    + transition.getLogScoreFor( layer, context, stateID, seq, startPos ); //transition
					
					dist = current - bwdMatrix[layer][context];
					
					dist*=dist;
					if( dist < bestDist ) {
						childIdx = stateID;
						state = container[0];
						newContext = container[1];
						bestDist = dist;
					}
				}
			}
			
			if( state >= 0 ) {
				if( path == null ) {
					((TrainableTransition)transition).addToStatistic( layer, context, childIdx, weight, seq, startPos );
					((TrainableState)states[state]).addToStatistic( startPos, startPos, weight, seq );
				} else {
					path.add( state );
				}
				context = newContext;
			} else {
				break;
			}
		} while ( true );
		
		return bwdMatrix[0][0];
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
	protected double baumWelch( int startPos, int endPos, double weight, Sequence seq ) throws Exception {
		fillFwdMatrix( startPos, endPos, seq );
		fillBwdOrViterbiMatrix( Type.BAUM_WELCH, startPos, endPos, weight, seq );
		return bwdMatrix[0][0];
	}
/*	
	public void setSkipInit(boolean skipInit){
		this.skipInit = skipInit;
	}
*/	
	
	public void train(DataSet data, double[] weights) throws Exception {
		if( !(trainingParameter instanceof MaxHMMTrainingParameterSet) ) {
			throw new IllegalArgumentException( "This kind of training is currently not supported." );
		} else {
			Transition bestTransition = null;
			Emission[] bestEmissions = null;
			double best = Double.NEGATIVE_INFINITY;
					
			int numberOfStarts = trainingParameter.getNumberOfStarts();
			AbstractTerminationCondition tc = ((MaxHMMTrainingParameterSet) trainingParameter).getTerminationCondition();
			
			Compute compute = new Compute( threads, this );
			compute.setDataSet( data, weights );
			
			Time time = Time.getTimeInstance( sostream );
			for( int it, start = 0; start < numberOfStarts; start++ ) {
				sostream.writeln( "start " + start + " ============================" );
				//init
				if(!skipInit){
					initialize( data, weights );
					compute.setParameters();
				}
				
				//iterate
				double old_value, new_value = Double.NEGATIVE_INFINITY;
				it = 0;
				time.reset();
				do {
					old_value = new_value;
					new_value = getLogPriorTerm() + compute.oneIteration();
										
					sostream.writeln( it++ + "\t" + time.getElapsedTime() + "\t" + new_value + "\t" + (new_value - old_value) );
					if( tc.doNextIteration( it, old_value, new_value, null, null, Double.NaN, time) ) {
						compute.estimateFromStatistics();
					} else {
						break;
					}					
				}while( true );
				
				//check for best
				if( new_value > best ) {
					best = new_value;
					if( numberOfStarts > 1 ) {
						bestEmissions = ArrayHandler.clone( emission );
						bestTransition = transition.clone();
					}
				}
			}
			sostream.writeln( "best result: " + best );
			if( bestEmissions != null ) {
				emission = bestEmissions;
				transition = bestTransition;
				createStates();
			}
			compute.stopThreads();
		}
	}
	
	private double doOneStep(DataSet data, double[] weights, int start, int end ) throws Exception{
		Sequence seq;
		double weight = 1;
		double newValue = 0;
		for( int n = start; n < end; n++ ) {
			seq = data.getElementAt( n );
			if( weights != null ) {
				weight = weights[n];
			}
			
			if( trainingParameter instanceof ViterbiParameterSet ) {
				newValue += viterbi( null, 0, seq.getLength()-1, weight, seq ); //viterbi
			} else if ( trainingParameter instanceof BaumWelchParameterSet ) {
				newValue += baumWelch( 0, seq.getLength()-1, weight, seq ); //Baum-Welch
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
	 * @param data the data set
	 * @param weight the weights for each sequence of the data set
	 * 
	 * @throws Exception if an error occurs during the initialization 
	 */
	protected void initialize( DataSet data, double[] weight ) throws Exception {
		initializeRandomly();
	}
	
	/**
	 * This method initializes all emissions and the transition randomly.
	 */
	protected void initializeRandomly() {
		transition.initializeRandomly();
		for( int e = 0; e < emission.length; e++ ) {
			emission[e].initializeFunctionRandomly();
		}
	}
	
	/**
	 * This method resets all sufficient statistics of all emissions and the transition.
	 */
	protected void resetStatistics() {
		((TrainableTransition) transition).resetStatistic();
		for( int e = 0; e < emission.length; e++ ) {
				emission[e].resetStatistic();
		}
	}
	
	/**
	 * This method estimates the parameters of all emissions and the transition using their sufficient statistics.
	 */
	protected void estimateFromStatistics() {
		((TrainableTransition) transition).estimateFromStatistic();
		for( int e = 0; e < emission.length; e++ ) {
			emission[e].estimateFromStatistic();
		}
	}

	public final byte getMaximalMarkovOrder() throws UnsupportedOperationException {
		return Byte.MAX_VALUE;
	}
	
	public ResultSet getCharacteristics() throws Exception {
		return new ResultSet( getNumericalCharacteristics().getResults(),
				new Result[] { new StorableResult("model", "the xml representation of the model", this) } );
	}
	
	public String getInstanceName() {
		return "HMM(" + transition.getMaximalMarkovOrder() + ") " + trainingParameter.getClass().getSimpleName();
	}
	
	@Override
	public double[] getLogScoreFor(DataSet data) throws Exception {
		double[] logProb = new double[data.getNumberOfElements()];
		getLogScoreFor(data, logProb);
		return logProb;
	}

	@Override
	public void getLogScoreFor(DataSet data, double[] res) throws Exception {
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
	
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
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
	public void samplePath( IntList path, int startPos, int endPos, Sequence seq ) throws Exception {
		fillBwdMatrix( startPos, endPos, seq );
		
		int l = 0, stateID, context = 0, n;
		double logTransition;
		provideMatrix( 0, endPos-startPos+1 );
		path.clear();
		
		//iterate
		while( startPos <= endPos ) {
			
			//compute
			n = transition.getNumberOfChildren( l, context );
			for( stateID = 0; stateID < n; stateID++ ) {
				transition.fillTransitionInformation( l, context, stateID, container );			
				logTransition = transition.getLogScoreFor( l, context, stateID, seq, startPos );
						
				backwardIntermediate[stateID] = bwdMatrix[l+container[2]][container[1]] //old part
				       + states[container[0]].getLogScoreFor( startPos, startPos, seq ) //emission
				       + logTransition; //transition
			}
			Normalisation.logSumNormalisation( backwardIntermediate, 0, n );
			
			//draw
			stateID = AbstractMixtureTrainSM.draw( backwardIntermediate, 0 );
			transition.fillTransitionInformation( l, context, stateID, container );
			path.add( container[0] );
			context = container[1];
			l += container[2];
			startPos += container[2];
		}
		
		//final silent states
		int add, last = path.get( path.length()-1 );
		do {
			//compute
			n = transition.getNumberOfChildren( l, context );
			stateList.clear();
			for( stateID = 0; stateID < n; stateID++ ) {
				transition.fillTransitionInformation( l, context, stateID, container );
				if( states[container[0]].isSilent() ) {
					logTransition = transition.getLogScoreFor( l, context, stateID, seq, startPos );
					backwardIntermediate[stateList.length()] = bwdMatrix[l][container[1]] //old part
							       // there is no emission (silent state)
							       + logTransition; //transition
					stateList.add( stateID );
				}
			}
			if( finalState[last] ) {
				backwardIntermediate[stateList.length()]=0;
				add = 1;
			} else {
				add = 0;
			}
			Normalisation.logSumNormalisation( backwardIntermediate, 0, stateList.length()+add );
			
			//draw
			n = AbstractMixtureTrainSM.draw( backwardIntermediate, 0 );
			if( add == 1 && n == stateList.length() ) {
				break;
			} else {
				transition.fillTransitionInformation( l, context, stateList.get(n), container );
				path.add( container[0] );
				context = container[1];
				l += container[2];
				startPos += container[2];
				
				last = container[0];
			}
		} while( true );
	}
		
	private double computeLogScoreFromForward( int l ) {
		double res = Double.NEGATIVE_INFINITY;
		if( transition.getMaximalMarkovOrder() > 0 ) {
			for( int i = 0; i < fwdMatrix[l].length; i++ ) {
				if( finalState[transition.getLastContextState(l, i)] ) {
					res = Normalisation.getLogSum( res, fwdMatrix[l][i] );
				}
			}
		} else {
			res = Normalisation.getLogSum( fwdMatrix[l] );
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
	
	private static class Compute {
		private static enum WorkerState{
			TRAIN,
			WAIT,
			STOP
		}
		
		private class WorkerThread extends Thread{

			private boolean exception;
			private WorkerState state;
			private int idx;
			private int start, end;
			private double score;
			private DataSet data;
			private double[] weights;
			private HigherOrderHMM hmm;

			public WorkerThread( int idx, HigherOrderHMM hmm ) {
				this.idx = idx;
				this.hmm = hmm;
				this.setDaemon( true );
				this.state = WorkerState.WAIT;
				start();
			}
			
			private void set( int start, int end, DataSet data, double[] weights){
				this.start = start;
				this.end = end;
				this.score = 0;
				this.data = data;
				
				/*
				int n = 0;
				for( int i = start; i < end; i++){
					n+= data.getElementAt(i).getLength();
				}
				System.out.println( idx + "\t" + start + "\t" + end + "\t" + n );
				System.out.flush();
				/**/
				this.weights = weights;
				
				this.state = WorkerState.WAIT;
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
				exception = false;
				//System.out.println(idx+" : "+state);
				while(state != WorkerState.STOP){
					//System.out.println(idx+" ; "+state);
					if(state == WorkerState.WAIT){
						try {
							wait();
						} catch ( InterruptedException e ) {}
					}else{
						try{
							score = hmm.doOneStep( data, weights, start, end );
						}catch( Exception e ){
							exception = true;
							e.printStackTrace();
						}
						synchronized( Compute.this )
						{
							state = WorkerState.WAIT;
							Compute.this.notify();
						}
					}
				}
			}
			
			public boolean isWaiting(){
				return state == WorkerState.WAIT;
			}	
		}
		
		WorkerThread[] workers;
		TrainableTransition[] transition;
		Emission[] emission;
		
		public Compute( int threads, HigherOrderHMM hmm ) throws CloneNotSupportedException {
			workers = new WorkerThread[threads];
			workers[0] = new WorkerThread(0,hmm);
			for( int i = 1 ; i < threads; i++ ) {
				workers[i] = new WorkerThread(i,hmm.clone());
			}
			transition = new TrainableTransition[threads];
			emission = new Emission[threads];
		}
		
		private synchronized void waitUntilWorkersFinished(){
			int i, t = -1;
			boolean exception = false;
			while( true )
			{
				i = 0;
				while( i < workers.length && workers[i].isWaiting() ){
					if( workers[i].exception ) {
						t = i;
						exception = true;
					}
					i++;
				}
				if( i == workers.length ){
					if( exception ) {
						for( i = 0; i < workers.length; i++ ) {
							workers[i].interrupt();
						}
						stopThreads();
						throw new RuntimeException( "Terminate program, since at least thread " + t + " throws an exception." );
					} else {
						//System.out.println( "raus" );
						break;
					}
				}else{
					try{
						wait();
					} catch( InterruptedException e ) { }
				}
			}
		}
		
		/**
		 * This method can and should be used to stop all threads if they are not needed any longer.
		 */
		private void stopThreads()
		{
			for( int i = 0; i < workers.length; i++ )
			{
				workers[i].setState( WorkerState.STOP );
				workers[i] = null;
			}
		}
		
		private void setDataSet( DataSet data, double[] weights ) {
			int last = 0, N = data.getNumberOfElements();
			for(int i=0;i<workers.length-1;i++){
				workers[i].set( last, (i+1)*N/workers.length, data, weights );
				last = workers[i].end;
			}
			workers[ workers.length-1 ].set( last, N, data, weights );
		}
		
		private double oneIteration() {
			for(int i=0;i<workers.length;i++){
				workers[i].hmm.resetStatistics();
				workers[i].setState( WorkerState.TRAIN );
			}
			waitUntilWorkersFinished();
			double res = 0;
			for(int i=0;i<workers.length;i++){
				res += workers[i].getScore();
			}
			return res;
		}
		
		private void estimateFromStatistics() {
			if( workers.length > 1 ) {
				for( int i = 0; i < workers.length; i++ ) {
					transition[i] = (TrainableTransition)workers[i].hmm.transition;
				}
				((TrainableTransition)workers[0].hmm.transition).joinStatistics( transition );
				
				for(int e = 0; e < workers[0].hmm.emission.length; e++){
					for(int j=0;j<workers.length;j++){
						emission[j] = workers[j].hmm.emission[e];
					}
					workers[0].hmm.emission[e].joinStatistics( emission );
				}
			}
			workers[0].hmm.estimateFromStatistics();
			setParameters();
		}
		
		private void setParameters() {
			for( int i = 1; i < workers.length; i++ ) {
				workers[i].hmm.transition.setParameters( workers[0].hmm.transition );
				for(int e = 0; e < workers[0].hmm.emission.length; e++){
					workers[i].hmm.emission[e].setParameters( workers[0].hmm.emission[e] );
				}
			}
		}
	}
}