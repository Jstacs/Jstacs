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
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.models.hmm.models;

import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.NonParsableException;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;
import de.jstacs.models.NormalizableScoringFunctionModel;
import de.jstacs.models.hmm.State;
import de.jstacs.models.hmm.states.DifferentiableState;
import de.jstacs.models.hmm.states.SimpleDifferentiableState;
import de.jstacs.models.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.models.hmm.training.MaxHMMTrainingParameterSet;
import de.jstacs.models.hmm.training.NumericalHMMTrainingParameterSet;
import de.jstacs.models.hmm.transitions.DifferentiableTransition;
import de.jstacs.models.hmm.transitions.elements.TransitionElement;
import de.jstacs.scoringFunctions.SamplingScoringFunction;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;


/**
 * This class combines an {@link HigherOrderHMM} and a {@link de.jstacs.scoringFunctions.NormalizableScoringFunction} by implementing some of the declared methods. 
 * 
 * @author Jens Keilwagen
 */
public class DifferentiableHigherOrderHMM extends HigherOrderHMM implements SamplingScoringFunction {	
	
	/**
	 * The number of parameters of this HMM
	 */
	protected int numberOfParameters;
	
	/**
	 * The equivalent sample size used for the prior
	 */
	protected double ess;
	
	/**
	 * The type of the score that is evaluated
	 */
	protected Type score;
	
	/**
	 * Index array used for computing the gradient
	 */
	protected int[][] index;
	
	/**
	 * Help array for the gradient
	 */
	protected double[][][] gradient;
	
	/**
	 * Help array for the indexes of the parameters of the states
	 */
	protected IntList[] indicesState; 
	/**
	 * Help array for the indexes of the parameters of the transition
	 */
	protected IntList[] indicesTransition;
	/**
	 * Help array for the derivatives of the parameters of the states
	 */
	protected DoubleList[] partDerState; 
	/**
	 * Help array for the derivatives of the parameters of the transition
	 */
	protected DoubleList[] partDerTransition;
	
	/**
	 * This is the main constructor.
	 * 
	 * @param trainingParameterSet the {@link de.jstacs.parameters.ParameterSet} that determines the training algorithm and contains the necessary {@link de.jstacs.parameters.Parameter}s
	 * @param name the names of the states
	 * @param emissionIdx the indices of the emissions that should be used for each state, if <code>null</code> state <code>i</code> will use emission <code>i</code>
	 * @param forward a boolean array that indicates whether the symbol on the forward or the reverse complementary strand should be used,
	 * 				  if <code>null</code> all states use the forward strand
	 * @param emission the emissions
	 * @param likelihood if <code>true</code> the likelihood is return  by {@link #getLogScoreFor(Sequence)} otherwise the viterbi score
	 * @param ess the ess of the model
	 * @param te the {@link TransitionElement}s used for creating a {@link de.jstacs.models.hmm.Transition}
	 * 
	 * @throws Exception if 
	 * 	<ul>
	 *  <li>some component could not be cloned</li> 
	 *  <li>some the length of <code>name, emissionIdx,</code> or <code>forward</code> is not equal to the number of states</li>
	 *  <li>not all emissions use the same {@link de.jstacs.data.AlphabetContainer}</li>
	 *  <li>the states can not be handled by the transition
	 *  </ul>
	 */
	public DifferentiableHigherOrderHMM( MaxHMMTrainingParameterSet trainingParameterSet, String[] name, int[] emissionIdx, boolean[] forward,
			DifferentiableEmission[] emission, boolean likelihood, double ess, TransitionElement... te ) throws Exception {
		super( trainingParameterSet, name, emissionIdx, forward, emission, te );
		getOffsets();
		this.score = likelihood ? Type.LIKELIHOOD : Type.VITERBI;
		if( ess < 0 ) {
			throw new IllegalArgumentException();
		}
		this.ess = ess;
	}	
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link DifferentiableHigherOrderHMM} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DifferentiableHigherOrderHMM} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public DifferentiableHigherOrderHMM( StringBuffer xml ) throws NonParsableException {
		super( xml );
		getOffsets();
	}
	
	@Override
	protected void appendFurtherInformation( StringBuffer xml ) {
		super.appendFurtherInformation( xml );
		XMLParser.appendObjectWithTags( xml, ess, "ess" );
		XMLParser.appendObjectWithTags( xml, score, "score" );
	}

	@Override
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		super.extractFurtherInformation( xml );
		ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
		score = XMLParser.extractObjectForTags( xml, "score", Type.class );

	}

	protected void createHelperVariables() {
		if( container == null ) {
			int maxOrder = transition.getMaximalMarkovOrder(), anz = 0, i;
			for( i = 0; i <= maxOrder; i++ ) {
				anz = Math.max( anz, transition.getNumberOfIndexes( i ) );
			}
			if( gradient == null || gradient[0].length != anz || gradient[0][0].length != numberOfParameters ) {
				gradient = new double[2][anz][numberOfParameters];
				index = new int[3][anz];
			}
			if( indicesState == null ) {
				anz = transition.getMaximalNumberOfChildren();
				try {
					indicesState = ArrayHandler.createArrayOf( new IntList(), states.length );
					partDerState = ArrayHandler.createArrayOf( new DoubleList(), states.length );
					
					indicesTransition = ArrayHandler.createArrayOf( new IntList(), anz );
					partDerTransition = ArrayHandler.createArrayOf( new DoubleList(), anz );
				} catch( CloneNotSupportedException cnse ) {
					throw getRunTimeException( cnse );
				}
			}
		}
		super.createHelperVariables();
	}
	
	protected void createStates() {
		this.states = new SimpleDifferentiableState[emissionIdx.length];
		for( int i = 0; i < emissionIdx.length; i++ ) {
			this.states[i] = new SimpleDifferentiableState( (DifferentiableEmission) emission[emissionIdx[i]], name[i], forward[i] );
		}
	}

	public DifferentiableHigherOrderHMM clone() throws CloneNotSupportedException {
		//prepare for clone
		double[][][] grad = gradient;
		gradient = null;
		IntList[] ind = indicesState;
		indicesState = null;
		//clone
		DifferentiableHigherOrderHMM clone = (DifferentiableHigherOrderHMM) super.clone();
		//reverse
		gradient = grad;
		indicesState = ind;
		return clone;
	}
	
	public double getESS() {
		return ess;
	}

	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception {
		for( int e = 0; e < emission.length; e++ ) {
			((DifferentiableEmission)emission[e]).addGradientOfLogPriorTerm( grad, start );
		}
		((DifferentiableTransition) transition).addGradientForLogPriorTerm( grad, start );
	}
	
	protected void getOffsets() {
		numberOfParameters = 0;
		for( int e = 0; e < emission.length; e++ ) {
			numberOfParameters = ((DifferentiableEmission)emission[e]).setParameterOffset( numberOfParameters );
			if( numberOfParameters == UNKNOWN ) {
				return;
			}
		}
		numberOfParameters = ((DifferentiableTransition)transition).setParameterOffset( numberOfParameters );
		if( numberOfParameters == UNKNOWN ) {
			return;
		}
		createHelperVariables();
	}
	
	public int getNumberOfParameters() {
		return numberOfParameters;
	}

	public int getNumberOfRecommendedStarts() {
		return trainingParameter.getNumberOfStarts();
	}

	public double[] getCurrentParameterValues() throws Exception {
		int i = 0, n = getNumberOfParameters();
		
		if( n != UNKNOWN ) {
			double[] params = new double[n];
			for( int e = 0; e < emission.length; e++, i++ ) {
				((DifferentiableEmission)emission[e]).fillCurrentParameter( params );
			}
			((DifferentiableTransition)transition).fillParameters( params );
			return params;
		} else {
			throw new IllegalArgumentException();
		}
	}
	
	public boolean isInitialized() {
		return true;
	}
/*	
	public void setTrainingParameters(MaxHMMTrainingParameterSet params) throws CloneNotSupportedException{
		this.trainingParameter = (HMMTrainingParameterSet)params.clone();
	}
*/
	public void setParameters( double[] params, int start ) {
		for( int e = 0; e < emission.length; e++ ) {
			((DifferentiableEmission)emission[e]).setParameter( params, start );
		}
		((DifferentiableTransition)transition).setParameters( params, start );
	}
	
	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {
		if(skipInit){
			return;
		}
		initializeRandomly();
		getOffsets();
	}
	
	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		if(skipInit){
			return;
		}
		if( trainingParameter instanceof NumericalHMMTrainingParameterSet ) {
			initializeFunctionRandomly( freeParams );
		} else {
			train( data[index], weights==null? null : weights[index] );
			getOffsets();
		}
	}
	
	public void train( DataSet data, double[] weights ) throws Exception {
		if( trainingParameter instanceof NumericalHMMTrainingParameterSet ) {
			NumericalHMMTrainingParameterSet params = (NumericalHMMTrainingParameterSet) trainingParameter;
			NormalizableScoringFunctionModel model = new NormalizableScoringFunctionModel( this, params.getNumberOfThreads(), params.getAlgorithm(), params.getTerminationCondition(), params.getLineEps(), params.getStartDistance() );
			model.setOutputStream( sostream );
			model.train( data, weights );
			
			DifferentiableHigherOrderHMM hmm = (DifferentiableHigherOrderHMM) model.getFunction();
			this.emission = hmm.emission;
			createStates();
			this.transition = hmm.transition;	
		} else {
			super.train( data, weights );
		}
	}

//XXX is normalized? start
	public boolean isNormalized() {
		return true;
	}
	
	public double getLogNormalizationConstant() {
		return 0;
	}

	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception {
		return Double.NEGATIVE_INFINITY;
	}
	
	public double getInitialClassParam(double classProb) {
		return Math.log( classProb );
	}
//end

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#getLogScoreFor(de.jstacs.data.Sequence)
	 */
	@Override
	public double getLogScoreFor( Sequence seq ) {
		return getLogScoreFor( seq, 0 );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#getLogScoreFor(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor( Sequence seq, int start ) {
		return getLogScoreFor( seq, start, seq.getLength()-1 );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#getLogScoreFor(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor( Sequence seq, int start, int end ) {
		return logProb( start, seq.getLength()-1, seq );
	}
	
	protected double logProb( int startpos, int endpos, Sequence sequence ) {
		try {
			fillBwdOrViterbiMatrix( score, startpos, endpos, 0, sequence );
		} catch( Exception e ) {
			throw getRunTimeException( e );
		}

		return bwdMatrix[0][0];
	}
	
	public double getLogScoreAndPartialDerivation( Sequence seq, IntList indices, DoubleList partialDer ) {
		return getLogScoreAndPartialDerivation( seq, 0, indices, partialDer );
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int startPos, IntList indices, DoubleList partialDer ) {
		return getLogScoreAndPartialDerivation( seq, startPos, seq.getLength()-1, indices, partialDer);
	}
		
	public double getLogScoreAndPartialDerivation( Sequence seq, int startPos, int endPos, IntList indices, DoubleList partialDer ) {
		try {
			boolean zero = transition.getMaximalMarkovOrder() == 0;
			int l = endPos-startPos+1, stateID, context, n, children;
			
			provideMatrix( 1, endPos-startPos+1 );
			
			for( int idx2 = 0; idx2 < gradient[1].length; idx2++ ) {
				Arrays.fill( gradient[0][idx2], 0 );
				Arrays.fill( gradient[1][idx2], 0 );
			}
			DifferentiableTransition diffTransition = (DifferentiableTransition) transition;
			
			//init
			double val;
			for( stateID = 0; stateID < states.length; stateID++ ) {
				indicesState[stateID].clear();
				partDerState[stateID].clear();
			}
			for( context = bwdMatrix[l].length-1; context >= 0; context-- ) {
				n = transition.getNumberOfChildren( l, context );

				children = 0;
				if( zero || finalState[transition.getLastContextState( l, context )] ) {
					val = 0;
				} else {
					val = Double.NEGATIVE_INFINITY;
				}
				//for all different children states
				for( stateID = 0; stateID < n; stateID++ ) {
					transition.fillTransitionInformation( l, context, stateID, container );
					if( states[container[0]].isSilent() ) {
						indicesTransition[children].clear();
						partDerTransition[children].clear();
						
						backwardIntermediate[children] =
							bwdMatrix[l][container[1]] //backward score until next position
							//there is no emission (silent state)
						    + diffTransition.getLogScoreAndPartialDerivation( l, context, stateID, indicesTransition[children], partDerTransition[children], seq, endPos ); //transition
						
						if( backwardIntermediate[children] != Double.NEGATIVE_INFINITY ) {	
							index[0][children] = container[0];
							index[1][children] = container[1];
							index[2][children] = container[2];
							children++;
						}
					}
				}
				if( children == 0 ) {
					bwdMatrix[l][context] = val;
					resetGradient( l, context, 0 );
				} else {
					merge( children, l, context, val );
				}
			}
			//System.out.println( seq.toString(endPos, endPos+1) + "\t" + l + "\t" + Arrays.toString( bwdMatrix[l] ) );
			
			//compute scores for all positions backward
			while( --l >= 0 ) {
				for( stateID = 0; stateID < states.length; stateID++ ) {
					indicesState[stateID].clear();
					partDerState[stateID].clear();
					logEmission[stateID] = ((DifferentiableState) states[stateID]).getLogScoreAndPartialDerivation( endPos, endPos, indicesState[stateID], partDerState[stateID], seq );
				}
				//for all different contexts
				for( context = bwdMatrix[l].length-1; context >= 0; context-- ) {
					n = transition.getNumberOfChildren( l, context );
					//for all different children states
					children = 0;
					for( stateID = 0; stateID < n; stateID++ ) {
						indicesTransition[children].clear();
						partDerTransition[children].clear();
						
						transition.fillTransitionInformation( l, context, stateID, container );
						
						backwardIntermediate[children] =
							bwdMatrix[l+container[2]][container[1]] //backward score until next position
							+ logEmission[container[0]] //emission
							+ diffTransition.getLogScoreAndPartialDerivation( l, context, stateID, indicesTransition[children], partDerTransition[children], seq, endPos ); //transition
				
						if( backwardIntermediate[children] != Double.NEGATIVE_INFINITY ) {
							index[0][children] = container[0];
							index[1][children] = container[1];
							index[2][children] = container[2];
							children++;
						}
					}
					if( children == 0 ) {
						bwdMatrix[l][context] = Double.NEGATIVE_INFINITY;
						resetGradient( l, context, 0 );
					} else {
						merge( children, l, context, Double.NEGATIVE_INFINITY );
					}
				}
				endPos--;

				//System.out.println( (l==0?" ":seq.toString(endPos, endPos+1)) + "\t" + l + "\t" + Arrays.toString( bwdMatrix[l] ) );
			}
			
			for( int p = 0; p < numberOfParameters; p++ ) {
				if( gradient[0][0][p] != 0 ) {
					indices.add( p );
					partialDer.add( gradient[0][0][p] );
				}
			}
			return bwdMatrix[0][0];
		} catch( Exception e ) {
			throw getRunTimeException( e );
		}
	}
	
	private void merge( int anz, int layer, int context, double extra ) {
		int h = layer % 2;
		if( score == Type.VITERBI ) {
			int idx = ToolBox.getMaxIndex( 0, anz, backwardIntermediate );
			if( backwardIntermediate[idx] > extra ) {
				System.arraycopy( gradient[(layer+index[2][idx]) % 2][index[1][idx]], 0, gradient[h][context], 0, numberOfParameters );
				miniMerge( idx, 1, h, context );
				
				bwdMatrix[layer][context] = backwardIntermediate[idx];
			} else {
				bwdMatrix[layer][context] = extra;
			}
		} else { //LIKELIHOOD
			if( extra != Double.NEGATIVE_INFINITY ) {
				bwdMatrix[layer][context] = Normalisation.logSumNormalisation( backwardIntermediate, 0, anz, new double[]{extra}, backwardIntermediate, 0 );
			} else {
				bwdMatrix[layer][context] = Normalisation.logSumNormalisation( backwardIntermediate, 0, anz, backwardIntermediate, 0 );
			}
			
			// S = \sum_i u_i v_i
			// \frac{\partial \log(S)}{\partial \lambda}
			// = \sum_i 
			//		(u_i v_i / S) \frac{\partial \log u_i}{\partial \lambda}   (AAA)
			//	 	+ (u_i v_i / S) \frac{\partial \log v_i}{\partial \lambda} (BBB)
			
			//help[0][0][i] = u_i v_i / S
			
			Arrays.fill( gradient[h][context], 0 );
			
			//slow version
			/*
			// old = (AAA)
			for( int p = 0; p < numberOfParameters; p++ ) {
				for( int i = 0; i < anz; i++ ) {
					int x = (layer+index[2][i]) % 2;
					gradient[h][context][p] += help[0][0][i] * gradient[x][index[1][i]][p];
				}
			}	
			// transition & emission = (BBB)
			for( int i = 0; i < anz; i++ ) {
				miniMerge( i, help[0][0][i], h, context );
			}
			*/
			
			//fast version
			for( int i = 0; i < anz; i++ ) {
				// old = (AAA)
				int x = (layer+index[2][i]) % 2;
				for( int p = 0; p < numberOfParameters; p++ ) {
					gradient[h][context][p] += backwardIntermediate[i] * gradient[x][index[1][i]][p];
				}

				// transition & emission = (BBB)
				miniMerge( i, backwardIntermediate[i], h, context );
			}
		}
	}
	
	//add the partial derivation with given weight
	private void miniMerge( int i, double weight, int h, int context ) {
		for( int p = 0; p < indicesTransition[i].length(); p++ ) {
			gradient[h][context][indicesTransition[i].get(p)] += weight * partDerTransition[i].get(p); 
		}
		for( int p = 0; p < indicesState[index[0][i]].length(); p++ ) {
			gradient[h][context][indicesState[index[0][i]].get(p)] += weight * partDerState[index[0][i]].get(p); 
		}
	}
	
	private void resetGradient( int layer, int context, double val ) {
		Arrays.fill( gradient[layer % 2][context], val );
	}
	
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		int off = 0;
		for(int i=0;i<emission.length;i++){
			int num = ((DifferentiableEmission)emission[i]).getNumberOfParameters();
			if( num > 0){
				if(index >= off && index < off + num){
					return ((DifferentiableEmission)emission[i]).getSizeOfEventSpace();
				}
			}
			off += num;
		}
		return ((DifferentiableTransition)transition).getSizeOfEventSpace(index);
	}
	
	@Override
	public int[][] getSamplingGroups( int parameterOffset ) {
		LinkedList<int[]> list = new LinkedList<int[]>();
		for(int i=0;i<emission.length;i++){
			((DifferentiableEmission)emission[i]).fillSamplingGroups(parameterOffset, list);
		}
		((DifferentiableTransition)transition).fillSamplingGroups(parameterOffset, list);
		return list.toArray( new int[0][0] );
	}
	
	public String getInstanceName() {
		return "differentiable HMM(" + transition.getMaximalMarkovOrder() + ", " + score + ")";
	}
}