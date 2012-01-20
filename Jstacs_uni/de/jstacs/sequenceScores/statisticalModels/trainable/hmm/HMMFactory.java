package de.jstacs.sequenceScores.statisticalModels.trainable.hmm;

import java.lang.reflect.Constructor;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.io.ArrayHandler;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.SequenceIterator;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.DifferentiableHigherOrderHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.HigherOrderHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.SamplingHigherOrderHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.SamplingPhyloHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.SamplingEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.SilentEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.UniformEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete.AbstractConditionalDiscreteEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete.DiscreteEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete.PhyloDiscreteEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete.ReferenceSequenceDiscreteEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.HMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.MaxHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.SamplingHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement;
import de.jstacs.sequenceScores.statisticalModels.trainable.phylo.PhyloTree;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;


/**
 * This class allows to create some frequently used HMMs.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class HMMFactory {

	/**
	 * This enum defines some standard architecture of profile HMMs.
	 * 
	 * @author Jan Grau, Jens Keilwagen
	 * 
	 * @see HMMFactory#createProfileHMM(MaxHMMTrainingParameterSet, HMMType, boolean, int, int, AlphabetContainer, double, boolean, boolean, double[][])
	 */
	public static enum HMMType {
		/**
		 * The PLAN7 architecture without connections between delete and insert states
		 */
		PLAN7,
		/**
		 * The PLAN9 architecture, where each state has transitions to the insert state of the same layer and
		 * the delete and match states of the next layer
		 */
		PLAN9,
		/**
		 * Architecture similar to PLAN7 but with connections from delete to insert states of the same layer
		 */
		PLAN8I,
		/**
		 * Architecture similar to PLAN7 but with connections from insert to delete states of the next layer
		 */
		PLAN8D;
	}
	
	private static AbstractHMM getHMM( HMMTrainingParameterSet pars, String[] name, int[] emissionIdx, boolean[] forward, Emission[] emission, TransitionElement[] te, double ess, boolean likelihood ) throws Exception {
		if( pars instanceof SamplingHMMTrainingParameterSet ) {
			if(emission instanceof PhyloDiscreteEmission[]) {
			    return new SamplingPhyloHMM( (SamplingHMMTrainingParameterSet) pars, name, emissionIdx, forward, ArrayHandler.cast( PhyloDiscreteEmission.class, emission ), te );
			} else {
			    return new SamplingHigherOrderHMM( (SamplingHMMTrainingParameterSet) pars, name, emissionIdx, forward, ArrayHandler.cast( SamplingEmission.class, emission ), te );
			}
		} else {
			boolean diff = true;//transition element okay?
			for( int e = 0; diff && e < emission.length; e++ ) {
				diff = emission[e] instanceof DifferentiableEmission;					
			}
			if( diff ) {
				return new DifferentiableHigherOrderHMM( (MaxHMMTrainingParameterSet)pars, name, emissionIdx, forward, ArrayHandler.cast( DifferentiableEmission.class, emission ), likelihood, ess, te );
			} else {
				return new HigherOrderHMM( pars, name, emission, te );
			}
		}
	}	
	
	private static void addTransitions( int[] states, double ess, double selfTranistionPart, int o, AbstractList<TransitionElement> list ) {
		int[] context = new int[o];
		SequenceIterator it = new SequenceIterator(o);
		Arrays.fill( context, states.length );
		it.setBounds( context );
		
		double[] hyperParams = new double[states.length];
		double e = ess*(1-selfTranistionPart)/(states.length-1);
		Arrays.fill( hyperParams, e );
		do {
			for( int i = 0; i < context.length; i++ ) {
				context[i] = it.discreteValAt(i);
			}
			hyperParams[context[o-1]] = ess*selfTranistionPart;
			list.add( new TransitionElement( context, states, hyperParams ) );
			hyperParams[context[o-1]] = e;
		} while( it.next() );
	}
	
	/**
	 * This method creates an ergodic, i.e. a completely connected, HMM using the given emissions.
	 * 
	 * @param pars the parameter determining the training procedure of the HMM
	 * @param order the Markov order of the HMM 
	 * @param ess the ess to be used for the HMM
	 * @param selfTranistionPart the a-priori probability of any self transition
	 * @param expectedSequenceLength the expected length of sequences to be modeled
	 * @param emission the emissions of the states
	 * 
	 * @return an ergodic HMM
	 * 
	 * @throws Exception if the HMM can not be created properly
	 */
	public static AbstractHMM createErgodicHMM( HMMTrainingParameterSet pars, int order, double ess, double selfTranistionPart, double expectedSequenceLength, Emission... emission ) throws Exception {
		int n = emission.length;
		
		String[] name = new String[n];
		int[] states = new int[n];
		for( int i = 0; i < n; i++ ) {
			name[i] = "" + i;
			states[i] = i;
			if( emission[i] instanceof SilentEmission ) {
				throw new IllegalArgumentException( "An ergodic HMM can not contain silent states." );
			}
		}

		double[] hyperParams = new double[n];
		LinkedList<TransitionElement> list = new LinkedList<TransitionElement>();
		double e = ess/n;
		Arrays.fill( hyperParams, e );
		list.add( new TransitionElement( null, states, hyperParams ) );
		for( int o = 1 ; o < order; o++ ) {
			addTransitions( states, e, selfTranistionPart, o, list );
			e /= n;
		}
		addTransitions( states, ess*(expectedSequenceLength-order), selfTranistionPart, order, list );
		
		return getHMM( pars, name, null, null, emission, list.toArray( new TransitionElement[0] ), ess, true );
	}

	/**
	 * Creates an HMM with <code>numStates+1</code> states, where <code>numStates</code> emitting build a clique and each of those states is connected to the absorbing silent final state.
	 * Such an HMM models the length of the input sequences, i.e. summing the likelihood over all input sequences (of different lengths) will give 1,
	 * whereas an ergodic HMM does not model the input length and summing the likelihood of each possible input length will give 1.
	 * 
	 * @param pars the parameters of the algorithm for learning the model parameters
	 * @param ess the equivalent sample size, is propagated between states to obtain consistent hyper-parameters for all parameters
	 * @param selfTranistionPart the a-priori probability of a self transition for each emitting state
	 * @param finalTranistionPart the a-priori probability of the transition to the final state from each emitting state 
	 * @param con the {@link AlphabetContainer} of the HMM
	 * @param numStates the number of emitting states
	 * @param insertUniform if <code>true</code> the emitting states will use {@link UniformEmission}s
	 * 
	 * @return an HMM with <code>numStates+1</code> states, where <code>numStates</code> emitting build a clique and each of those states is connected to the absorbing silent final state 
	 * 
	 * @throws Exception if the HMM could not be created properly
	 * 
	 * @see #propagateESS(double, ArrayList)
	 */
	public static AbstractHMM createPseudoErgodicHMM( HMMTrainingParameterSet pars, double ess, double selfTranistionPart, double finalTranistionPart, AlphabetContainer con, int numStates, boolean insertUniform ) throws Exception {
		Emission[] emission = new Emission[numStates+1];
		String[] name = new String[numStates+1];
		int[] states;
		
		ArrayList<PseudoTransitionElement> list = new ArrayList<PseudoTransitionElement>();
		double[] hyperParams;
		
		//start
		states = new int[numStates];
		for( int i = 0; i < states.length; i++ ){
			states[i] = i;
		}
		hyperParams = new double[states.length];
		double e = ess/hyperParams.length;
		Arrays.fill( hyperParams, e );
		list.add( new PseudoTransitionElement( null, states, hyperParams ) );
		
		//main
		states = new int[numStates+1];
		for( int i = 0; i < states.length; i++ ){
			states[i] = i;
		}
		hyperParams = new double[states.length];
		e = (1-selfTranistionPart-finalTranistionPart)*ess/(hyperParams.length-2);
		Arrays.fill( hyperParams, e );
		hyperParams[numStates] = finalTranistionPart*ess;
		for( int i = 0; i < numStates; i++ ) {
			hyperParams[i] = selfTranistionPart*ess;
			list.add( new PseudoTransitionElement( new int[]{i}, states, hyperParams ) );
			hyperParams[i] = e;			
		}
		
		//final
		name[numStates] = "F";
		states[numStates] = numStates-1;
		emission[numStates] = new SilentEmission();
		
		Pair<double[][], double[]> p = propagateESS(ess, list);
		double[] stateEss = p.getSecondElement(); 
		for( int i = 0; i < numStates; i++ ) {
			name[i] = "" + i;
			if( insertUniform ) {
				emission[i] = new UniformEmission( con );
			} else {
				emission[i] = new DiscreteEmission( con, stateEss[i] );
			}
		}
		return getHMM( pars, name, null, null, emission, createTransition( p.getFirstElement(), list ), ess, false );
	}
	
	/**
	 * This method creates a first order sunflower HMM.
	 * <b>The current implementation does not set any hyper parameters for the prior.</b> 
	 * 
	 * @param pars the parameter determining the training procedure of the HMM
	 * @param con the {@link AlphabetContainer} of the HMM
	 * @param ess the equivalent sample size (ess) of this model
	 * @param expectedSequenceLength the expected sequence length to be modeled; this parameter is used to determine the prior
	 * @param startCentral a switch for deciding between starting in the central state or in all states
	 * @param motifLength the length of the motifs building the petals of the sunflower
	 * 
	 * @return a first order sunflower HMM
	 * 
	 * @throws Exception if the HMM can not be created properly
	 * 
	 * @see #createSunflowerHMM(HMMTrainingParameterSet, AlphabetContainer, double, int, boolean, PhyloTree[], double[], int[])
	 */
	public static AbstractHMM createSunflowerHMM( HMMTrainingParameterSet pars, AlphabetContainer con, double ess, int expectedSequenceLength, boolean startCentral, int... motifLength ) throws Exception {
		return createSunflowerHMM( pars, con, ess, expectedSequenceLength, startCentral, null, null, motifLength );
	}

   /**
	 * This method creates a first order sunflower HMM allowing phylogenetic emissions.
	 * <b>The current implementation does not set any hyper parameters for the prior.</b>
	 *
	 * @param pars the parameter determining the training procedure of the HMM
	 * @param con the {@link AlphabetContainer} of the HMM
	 * @param ess the equivalent sample size (ess) of this model
	 * @param expectedSequenceLength the expected sequence length to be modeled; this parameter is used to determine the prior
	 * @param startCentral a switch for deciding between starting in the central state or in all states
	 * @param t an array of length two that contains a {@link PhyloTree} for the background and the motif, can be <code>null</code> than a normal sunflower HMM is returned
	 * @param motifProb the a-priori probabilities for each motif, i.e., the a-priori probabilities for the edges from the central node to the first motif states
	 * @param motifLength the length of the motifs building the petals of the sunflower
	 *
	 * @return a first order sunflower HMM
	 *
	 * @throws Exception if the HMM can not be created properly
	 */
	public static AbstractHMM createSunflowerHMM( HMMTrainingParameterSet pars, AlphabetContainer con, double ess, int expectedSequenceLength, boolean startCentral, PhyloTree[] t, double[] motifProb, int[] motifLength ) throws Exception {
		if( motifProb == null ) {
			motifProb = new double[motifLength.length];
			Arrays.fill( motifProb, 0.1/motifLength.length );
		}
		
		int anz = 1;
		int[] states = new int[1+motifLength.length];
		for( int i = 0; i < motifLength.length; i++ ) {
			states[1+i] = anz;
			anz += motifLength[i];
		}
		Emission[] e = t == null ? new Emission[anz] : new PhyloDiscreteEmission[anz];
		String[] name = new String[anz];
		LinkedList<TransitionElement> list = new LinkedList<TransitionElement>();
		
		//prepare propagation of ESS
		double[][] hyperNext = new double[motifLength.length+1][], stateESS;
		double self = 0;
		for( int i = 0; i < motifLength.length; i++ ) {
			hyperNext[i] = new double[motifLength[i]];
			self += motifProb[i];
		}
		self = 1-self;
		hyperNext[motifLength.length] = new double[1];
		stateESS = ArrayHandler.clone( hyperNext );
		
		double[] hyperParams;
		if( startCentral ) {
			hyperNext[motifLength.length][0] = ess;
			hyperParams = new double[]{ess};
		} else {
			hyperNext[motifLength.length][0] = ess*self;
			hyperParams = new double[anz];
			hyperParams[0] = hyperNext[motifLength.length][0];
			for( int i = 0; i < motifLength.length; i++ ) {
				Arrays.fill( hyperNext[i], motifProb[i]*ess/hyperNext[i].length );
				Arrays.fill( hyperParams, states[1+i], states[1+i]+motifLength[i], hyperNext[i][0] );
			}			
		}
		
		//start transition
		int[] children;
		if( startCentral ) {
			children = new int[]{0};
		} else {
			children = new int[anz];
			for( int i = 0; i < anz; i++) {
				children[i] = i;
			}
		}
		list.add( new TransitionElement( null, children, hyperParams ) );
		/*
		System.out.println( "S\t" + Arrays.toString(hyperParams) );
		System.out.println();
		*/
		//propagate ESS
		hyperParams = new double[motifLength.length+1];
		for( int j, i, l = 0; l < expectedSequenceLength; l++ ) {
			double d = 0;
			for( i = 0; i < motifLength.length; i++ ) {
				d += hyperNext[i][motifLength[i]-1];
				for( j = motifLength[i]-1; j > 0; j-- ) {
					stateESS[i][j] += hyperNext[i][j];
					hyperNext[i][j] = hyperNext[i][j-1];
				}
				//j==0
				stateESS[i][j] += hyperNext[i][j];
				hyperNext[i][j] = motifProb[i] * hyperNext[motifLength.length][j];
				hyperParams[1+i] += hyperNext[i][j];
			}
			stateESS[i][0] += hyperNext[i][0];
			hyperParams[0] += self*hyperNext[i][0];
			hyperNext[i][0] = d+self*hyperNext[i][0];
		}
		
		/*
		for( int i = 0; i < motifLength.length; i++ ) {
			System.out.println(i+"\t"+Arrays.toString(stateESS[i]));
		}
		System.out.println(motifLength.length+"\t"+stateESS[motifLength.length][0]);
		System.out.println();
		System.out.println( "X\t" + Arrays.toString(hyperParams) );
		*/
		
		//create states and transition
		e[0] = getEmission( con, stateESS[motifLength.length][0], t == null ? null : t[0] );
		name[0] = "bg";

		//transition to the petals
		list.add( new TransitionElement( new int[]{0}, states, hyperParams ) );

		//states and transitions for each motif 
		int idx = 1;
		hyperParams = new double[1];
		for( int m = 0 ; m < motifLength.length; m++ ) {
			hyperParams[0] = ess; //does not matter 
			for( int p = 0; p < motifLength[m]; p++, idx++ ) {
				e[idx] = getEmission( con, stateESS[m][p],  t == null ? null : t[1] );
				name[idx] = "motif " + m + " position " + p;
				if( p == motifLength[m]-1 ) {
					list.add( new TransitionElement( new int[]{idx}, new int[]{0}, hyperParams ) );
				} else {
					list.add( new TransitionElement( new int[]{idx}, new int[]{idx+1}, hyperParams ) );
				}
			}
		}

		return getHMM( pars, name, null, null, e, list.toArray( new TransitionElement[0] ), ess, true );
	}
	
	private static Emission getEmission( AlphabetContainer con, double ess, PhyloTree t ) {
		if( t != null ) {
			return new PhyloDiscreteEmission( con, ess, t );
		} else {
			return new DiscreteEmission( con, ess );
		}
	}
	
	/**
	 * Creates a new profile HMM for a given architecture and number of layers.
	 * @param trainingParameterSet the parameters of the algorithm for learning the model parameters
	 * @param type the type of the HMM, i.e., its architecture
	 * @param likelihood if <code>true</code>, the likelihood is considered as score of a sequence, and the probability of the Viterbi path otherwise
	 * @param order the order of the HMM, i.e., the number of previous states that are considered for a transition probability
	 * @param numLayers the number of layers of the profile HMM
	 * @param con the alphabet of the profile HMM
	 * @param ess the equivalent sample size, is propagated between states to obtain consistent hyper-parameters for all parameters
	 * @param conditionalMain if <code>true</code>, the match states have {@link ReferenceSequenceDiscreteEmission}s, and {@link DiscreteEmission}s otherwise
	 * @param closeCircle if <code>true</code> the circle from end to initial state is closed, i.e., the HMM can be traversed several times
	 * @param conditionInitProbs the hyper-parameters for initializing the match states if <code>conditionalMain</code> is <code>true. May be <code>null</code> for using the hyper-parameters of the prior
	 * @return the profile HMM
	 * @throws Exception if the profile HMM could not be created
	 */
	public static AbstractHMM createProfileHMM(MaxHMMTrainingParameterSet trainingParameterSet, HMMType type, boolean likelihood, int order, int numLayers, AlphabetContainer con, double ess, boolean conditionalMain, boolean closeCircle, double[][] conditionInitProbs) throws Exception{
		double[][] initFromTo = getInitFromTo(type, ess);
		return createProfileHMM( trainingParameterSet, initFromTo, likelihood, order, numLayers, con, ess, conditionalMain, closeCircle, conditionInitProbs, false );
	}

	/**
	 * Creates a new profile HMM for a given architecture and number of layers.
	 * 
	 * @param trainingParameterSet the parameters of the algorithm for learning the model parameters
	 * @param initFromTo hyper-parameters of the transition from each state (first dimension) of the current layer to each other state in the same layer (first three entries of the second dimension)
	 * 						and the next layer (next three entries in the second dimension). If a hyper-parameter is set to {@link Double#NaN}, the corresponding transition is not allowed
	 * @param likelihood if <code>true</code>, the likelihood is considered as score of a sequence, and the probability of the Viterbi path otherwise
	 * @param order the order of the HMM, i.e., the number of previous states that are considered for a transition probability
	 * @param numLayers the number of layers of the profile HMM
	 * @param con the alphabet of the profile HMM
	 * @param ess the equivalent sample size, is propagated between states to obtain consistent hyper-parameters for all parameters
	 * @param conditionalMain if <code>true</code>, the match states have {@link ReferenceSequenceDiscreteEmission}s, and {@link DiscreteEmission}s otherwise
	 * @param closeCircle if <code>true</code> the circle from end to initial state is closed, i.e., the HMM can be traversed several times
	 * @param conditionInitProbs the hyper-parameters for initializing the match states if <code>conditionalMain</code> is <code>true. May be <code>null</code> for using the hyper-parameters of the prior
	 * @param insertUniform if <code>true</code> the insert states will use {@link UniformEmission}s
	 * 
	 * @return the profile HMM
	 * @throws Exception if the profile HMM could not be created
	 */
	public static AbstractHMM createProfileHMM(MaxHMMTrainingParameterSet trainingParameterSet, double[][] initFromTo, boolean likelihood, int order, int numLayers, AlphabetContainer con, double ess, boolean conditionalMain, boolean closeCircle, double[][] conditionInitProbs, boolean insertUniform ) throws Exception{
		return createProfileHMM(trainingParameterSet, initFromTo, likelihood, order, numLayers, con, ess, conditionalMain, closeCircle?1:0, conditionInitProbs, insertUniform );
	}
	
	/**
	 * Creates a new profile HMM for a given architecture and number of layers.
	 * 
	 * @param trainingParameterSet the parameters of the algorithm for learning the model parameters
	 * @param initFromTo hyper-parameters of the transition from each state (first dimension) of the current layer to each other state in the same layer (first three entries of the second dimension)
	 * 						and the next layer (next three entries in the second dimension). If a hyper-parameter is set to {@link Double#NaN}, the corresponding transition is not allowed
	 * @param likelihood if <code>true</code>, the likelihood is considered as score of a sequence, and the probability of the Viterbi path otherwise
	 * @param order the order of the HMM, i.e., the number of previous states that are considered for a transition probability
	 * @param numLayers the number of layers of the profile HMM
	 * @param con the alphabet of the profile HMM
	 * @param ess the equivalent sample size, is propagated between states to obtain consistent hyper-parameters for all parameters
	 * @param conditionalMain if <code>true</code>, the match states have {@link ReferenceSequenceDiscreteEmission}s, and {@link DiscreteEmission}s otherwise
	 * @param joiningStates the number of states used in the joining arc, if not positive the profile HMM does not contain any joining states (i.e. the circle is not closed)
	 * @param conditionInitProbs the hyper-parameters for initializing the match states if <code>conditionalMain</code> is <code>true. May be <code>null</code> for using the hyper-parameters of the prior
	 * @param insertUniform if <code>true</code> the insert states will use {@link UniformEmission}s
	 * 
	 * @return the profile HMM
	 * @throws Exception if the profile HMM could not be created
	 */
	public static AbstractHMM createProfileHMM(MaxHMMTrainingParameterSet trainingParameterSet, double[][] initFromTo, boolean likelihood, int order, int numLayers, AlphabetContainer con, double ess, boolean conditionalMain, int joiningStates, double[][] conditionInitProbs, boolean insertUniform ) throws Exception{
		PseudoTransitionElement[] coreTransitionTemplate = null;
		AbstractList<Class<? extends DifferentiableEmission>> emList = new LinkedList<Class<? extends DifferentiableEmission>>();
		ArrayList<PseudoTransitionElement> list = new ArrayList<PseudoTransitionElement>();
		AbstractList<String> nameList = new LinkedList<String>();
		
		Class<? extends DifferentiableEmission> insertClass = insertUniform ? UniformEmission.class : DiscreteEmission.class;
		
		//add states of silent chain (end)
		int endContext[] = new int[order];
		for( int i = 0; i < order; i++ ) {
			nameList.add( "E" + i );
			emList.add( SilentEmission.class );
			endContext[i] = i;
		}
		
		//add silent chain (start)
		int[] startContext = new int[order], context;
		for( int i = 0; i < order; i++ ) {
			nameList.add( "S" + i );
			emList.add( SilentEmission.class );
			context = new int[i];
			for( int k = 0; k < i; k++ ) {
				context[k]=order+k;
			}
			list.add( new PseudoTransitionElement( context, new int[]{order+i}, new double[]{ess} ) );
			startContext[i] = order+i;
		}		

		ArrayList<int[]> lastLayer = new ArrayList<int[]>();
		int lastStartNodeIndex = emList.size()-1;
		int offset = emList.size();
		//TODO multiple cores?
		
		//createProfileHMMCore( numLayers, initFromTo, order, emList, nameList, list, offset, conditionalMain, coreTransitionTemplate, "" );
		createProfileHMMCore( numLayers, lastLayer, initFromTo, order, emList, nameList, list, offset, lastStartNodeIndex, conditionalMain, coreTransitionTemplate, "", insertClass );
		
		//connect with end chain
		int[] states = new int[6];
		Arrays.fill(states, -1);
		states[3] = emList.size()-1;
		
		ArrayList<int[]> newLayer = new ArrayList<int[]>();
		for( int i = 0; i < order; i++ ) {
			shiftContext( states, 3 );
			states[3] = i;
			addProfileTransitions(initFromTo, order, states, list, lastLayer, newLayer );
		}		
		
		//end state
		nameList.add( "F" );
		emList.add( SilentEmission.class );
				
		if(joiningStates>0){
			int e = emList.size();
			list.add( new PseudoTransitionElement( endContext, new int[]{e-1, e}, new double[]{ess/2.0,ess/2.0} ) );
			
			for( int i = 0; i < joiningStates; i++ ) {
				nameList.add( "J" + i );
				emList.add( insertClass );//TODO?
			}
			int[] myContext = endContext.clone();
			System.arraycopy( myContext, 1, myContext, 0, order-1 );
			myContext[order-1] = e;
			
			ContextContainer[] possible = new ContextContainer[joiningStates];
			for( int i = 0; i < joiningStates; i++ ) {
				possible[i] = new ContextContainer();
			}
			possible[0].addConditional( myContext );
			ContextContainer extra = new ContextContainer();
			
			//cycle in state
			double[] hyperParams = new double[]{ess/3.0,ess/3.0,ess/3.0};
			int[] children;
			for( int f = e, j = 0; j < joiningStates; j++, f++ ) {
				if( j+1 == joiningStates ) {
					children = new int[]{f,order};
					hyperParams = new double[]{ess/2.0,ess/2.0};
				} else {
					children = new int[]{f,f+1,order};
				}
				for( int k = 0; k < possible[j].size(); k++ ) {
					myContext = possible[j].get(k);
					list.add( new PseudoTransitionElement( myContext, children, hyperParams ) );
					for( int l = 0; l < children.length; l++ ) {
						int[] next = myContext.clone();
						System.arraycopy( next, 1, next, 0, order-1 );
						next[order-1] = children[l];
						if( children[l] == order ) {
							extra.addConditional( next );
						} else {
							possible[children[l]-e].addConditional( next );
						}
					}
				}
			}
			//go to start
			for( int k = 0; k < extra.size(); k++ ) {
				myContext = extra.get(k).clone();
				for( int i = 0; i < order-1; i++ ) {
					list.add( new PseudoTransitionElement( myContext, new int[]{order+i+1}, new double[]{ess} ) );
					System.arraycopy( myContext, 1, myContext, 0, order-1 );
					myContext[order-1] = order+i+1;
				}	
			}
		} else {
			list.add( new PseudoTransitionElement( endContext, new int[]{emList.size()-1}, new double[]{ess} ) );
		}
		
		Pair<double[][], double[]> p = propagateESS( ess, list );
		
		return getHMM( trainingParameterSet, nameList.toArray( new String[0] ), null, null, 
				getEmissions( nameList, emList, p.getSecondElement(), con, conditionInitProbs ),
				createTransition( p.getFirstElement(), list ), ess, likelihood );
	}
	
	private static class ContextContainer extends ArrayList<int[]> {
		public boolean addConditional( int[] newContext ) {
			int i = 0, j; 
			while( i < size() ) {
				int[] test = get(i);
				j = 0;
				while( j < newContext.length && newContext[j] == test[j] ) {
					j++;
				}
				if( j == newContext.length ) {
					return false;
				}
				i++;
			}
			add( newContext );
			return true;
		}
	}

	private static double[][] getInitFromTo( HMMType type, double ess ) {
		double[][] initFromTo = new double[3][];
		if(type == HMMType.PLAN9){
			initFromTo[0] = new double[]{Double.NaN,ess/3.0,Double.NaN,ess/3.0,Double.NaN,ess/3.0};
			initFromTo[1] = initFromTo[0].clone();
			initFromTo[2] = initFromTo[1].clone();
		}else if(type == HMMType.PLAN7){
			initFromTo[0] = new double[]{Double.NaN,Double.NaN,Double.NaN,ess/2.0,Double.NaN,ess/2.0};
			initFromTo[1] = new double[]{Double.NaN,ess/2.0,Double.NaN,Double.NaN,Double.NaN,ess/2.0};
			initFromTo[2] = new double[]{Double.NaN,ess/3.0,Double.NaN,ess/3.0,Double.NaN,ess/3.0};
		}else if(type == HMMType.PLAN8I){
			initFromTo[0] = new double[]{Double.NaN,ess/3.0,Double.NaN,ess/3.0,Double.NaN,ess/3.0};
			initFromTo[1] = new double[]{Double.NaN,ess/2.0,Double.NaN,Double.NaN,Double.NaN,ess/2.0};
			initFromTo[2] = new double[]{Double.NaN,ess/3.0,Double.NaN,ess/3.0,Double.NaN,ess/3.0};
		}else if(type == HMMType.PLAN8D){
			initFromTo[0] = new double[]{Double.NaN,Double.NaN,Double.NaN,ess/2.0,Double.NaN,ess/2.0};
			initFromTo[1] = new double[]{Double.NaN,ess/3.0,Double.NaN,ess/3.0,Double.NaN,ess/3.0};
			initFromTo[2] = new double[]{Double.NaN,ess/3.0,Double.NaN,ess/3.0,Double.NaN,ess/3.0};
		}
		return initFromTo;
	}
	
	private static DifferentiableEmission[] getEmissions(AbstractList<String> nameList, AbstractList<Class<? extends DifferentiableEmission>> emList, double[] ess, AlphabetContainer con, double[][] conditionInitAPrioriProbs) throws Exception{
		DifferentiableEmission[] em = new DifferentiableEmission[emList.size()];
		Iterator<Class<? extends DifferentiableEmission>> it = emList.iterator();
		Iterator<String> n = nameList.iterator();
		int i=0;
		int refIdx = 0;
		while( it.hasNext() ) {
			Class<? extends DifferentiableEmission> cl = it.next();
			String name = n.next(), shape;
			if(cl == SilentEmission.class){
				em[i] = new SilentEmission();
			}else if( AbstractConditionalDiscreteEmission.class.isAssignableFrom( cl ) ) {
				if(cl == ReferenceSequenceDiscreteEmission.class){
					if(conditionInitAPrioriProbs != null){
						double[][] conditionInitHyperpars = ArrayHandler.clone( conditionInitAPrioriProbs );
						for(int j=0;j<conditionInitHyperpars.length;j++){
							for(int k=0;k<conditionInitHyperpars[j].length;k++){
								conditionInitHyperpars[j][k] *= ess[i];
							}
						}
						em[i] = new ReferenceSequenceDiscreteEmission( con, con, refIdx, ess[i], conditionInitHyperpars );
					}else{
						em[i] = new ReferenceSequenceDiscreteEmission( con, con, refIdx, ess[i] );
					}
					refIdx++;
				} else {
					em[i] = cl.getConstructor( AlphabetContainer.class, double.class ).newInstance( con, ess[i] );
				}
				if( name.charAt(0) == 'M' ) {
					shape = "rect";
				} else {
					shape = "diamond";
				}
				((AbstractConditionalDiscreteEmission)em[i]).setShape( shape );
			}else{ 
				Constructor<? extends DifferentiableEmission>[] cons = (Constructor<? extends DifferentiableEmission>[])cl.getConstructors();
				for(int j=0;em[i] == null && j<cons.length;j++){
					Class[] pars = cons[j].getParameterTypes();
					switch( pars.length ) {
						case 1:
							if( pars[0] == AlphabetContainer.class) {
								em[i] = cons[j].newInstance( con );
							}
							break;
						case 2:
							if( pars[0] == AlphabetContainer.class && pars[1] == double.class){
								em[i] = cons[j].newInstance( con, ess[i] );
							}
					}
				}
				if(em[i] == null){
					throw new Exception("Unsupported emission class.");
				}
			}
			i++;
		}
		return em;
	}
	
	private static void shiftContext( int[] context ) {
		shiftContext( context, 1 );
	}
	
	private static void shiftContext( int[] context, int s ) {
		System.arraycopy( context, s, context, 0, context.length-s );
	}

	private static void createProfileHMMCore(int numLayers, ArrayList<int[]> lastLayer, double[][] initFromTo, int order, AbstractList<Class<? extends DifferentiableEmission>> emList, AbstractList<String> nameList, AbstractList<PseudoTransitionElement> list, int totalOffset, int lastStartNodeIndex, boolean conditionalMain, PseudoTransitionElement[] coreTransitionTemplate, String suffix, Class<? extends DifferentiableEmission> insertClass ){
		if( order < 1 ) {
			throw new IllegalArgumentException("The order of a profile HMM has to be at least 1.");
		}
		ArrayList<int[]> newLayer = new ArrayList<int[]>();
			
		//start layer
		emList.add( SilentEmission.class );
		nameList.add( "D"+0 + suffix );

		emList.add( insertClass );
		nameList.add( "I"+0 + suffix );
		
		//connect start layer with initial states
		int[] context = new int[order];
		for( int i = 0; i < order; i++ ) {
			context[i] = lastStartNodeIndex-order+i+1;
		}
		//create PTE
		if( Double.isNaN( initFromTo[0][1] ) ) {
			list.add( new PseudoTransitionElement( context, new int[]{emList.size()-2, emList.size()-1}, new double[]{0.5, 0.5} )  );
		} else {
			list.add( new PseudoTransitionElement( context, new int[]{emList.size()-2}, new double[]{1} )  );
		}
		
		
		//fill lastLayer
		lastLayer.clear();
		shiftContext( context );
		context[order-1] = emList.size()-2;
		lastLayer.add( context.clone() );
		if( Double.isNaN( initFromTo[0][1] ) ) {
			context[order-1] = emList.size()-1;
			lastLayer.add( context );
		}
		int[] states = {-1, -1, -1, totalOffset, totalOffset+1, -1};

		//main layers
		boolean[] createStates = new boolean[3];
		Arrays.fill( createStates, true );
		for( int i=0;i<numLayers;i++ ){
			if( i == numLayers-1 ) {
				for( int j = 0; j < createStates.length; j++ ) {
					createStates[j] = !Double.isNaN( initFromTo[j][3] ); 
				}
			}
			int last = emList.size();
			createEmissions( createStates, emList, nameList, i+1, conditionalMain, suffix, insertClass );
			shiftContext( states, 3 );
			for( int j = 3; j < states.length; j++ ) {
				states[j] = createStates[j-3] ? last++ : -1;
			}
			addProfileTransitions( initFromTo, order, states, list, lastLayer, newLayer );
		}
		
		//last layer
		emList.add( SilentEmission.class );
		nameList.add( "D"+(numLayers+1)+suffix);
		
		shiftContext( states, 3 );
		states[3] = emList.size()-1;
		states[4] = -1;
		states[5] = -1;
		addProfileTransitions( initFromTo, order, states, list, lastLayer, newLayer );	
	}
	
	private static void addProfileTransitions( double[][] initFromTo, int order, int[] states, List<PseudoTransitionElement> allTE, List<int[]> lastLayer, List<int[]> newLayer ) {		
		DoubleList hyper = new DoubleList();
		DoubleList weights = new DoubleList();
		IntList children = new IntList();
		
		int[] current, next;
		for( int o = order-1, a, k, s = 0; s < lastLayer.size(); s++ ) {
			current = lastLayer.get(s);
			
			k = 0; 
			while( current[o] != states[k] ) {
				k++;
			}
			
			hyper.clear();
			weights.clear();
			children.clear();
			
			a=0;
			for( int j = 0; j < initFromTo[k].length; j++ ) {
				if( !Double.isNaN( initFromTo[k][j] ) && states[j] >= 0 ) {
					hyper.add( initFromTo[k][j] );
					children.add( states[j] );
					if( j < 3 ) {
						weights.add( 1000 );
						a = children.length();
					} else {
						weights.add( 1 );
					}
				}
			}
			
			if( children.length() > 0 ) {
				allTE.add( new PseudoTransitionElement( current, children.toArray(), hyper.toArray(), weights.toArray() ) );
			}
			
			next = current.clone();
			shiftContext( next );
			for( int j = 0; j < children.length(); j++ ) {
				next[o] = children.get(j);
				addConditional( next, j < a ? lastLayer : newLayer );
			}
		}
		
		lastLayer.clear();
		lastLayer.addAll( newLayer );
		newLayer.clear();
	}

	
	private static boolean addConditional( int[] context, List<int[]> list ) {
		for( int j, i = 0; i < list.size(); i++ ) {
			int[] current = list.get(i);
			j = 0;
			while( j < context.length && current[j] == context[j] ) {
				j++;
			}
			if( j == context.length ) {
				return false;
			}
		}
		list.add( context.clone() );
		return true;
	}
	
	/*
	//can be used later
	private static void joinTransitions(AbstractList<PseudoTransitionElement> elements){
		LinkedList<PseudoTransitionElement> tes = new LinkedList<PseudoTransitionElement>();
		PseudoTransitionElement[] els = elements.toArray( new PseudoTransitionElement[0] );
		for(int i=0;i<els.length;i++){
			if(els[i] != null){
				int[] cont1 = els[i].context;
				double[] temp = els[i].prob;
				IntList children = new IntList();
				DoubleList hyperParams = new DoubleList();
				for(int j=0;j<els[i].states.length;j++){
					children.add( els[i].states[ j ] );
					hyperParams.add( temp[j] );
				}
				for(int j=i+1;j<els.length;j++){
					if(els[j] != null && (els[i] == els[j] || els[i].hasSameContext( els[j] ))){
						temp = els[j].prob;
						for(int k=0;k<els[j].states.length;k++){
							int ac = -1;
							if( (ac = children.addConditional( els[j].states[ k ] ) ) < 0){
								hyperParams.add( temp[k] );
							}
						}
						els[j] = null;
					}
				}
				tes.add( new PseudoTransitionElement( cont1, children.toArray(), hyperParams.toArray(), true ) );
			}
		}
		elements.clear();
		elements.addAll( tes );
	}/**/

	private static int createEmissions( boolean[] createState, AbstractList<Class<? extends DifferentiableEmission>> emList, AbstractList<String> nameList, int offset, boolean conditionalMain, String suffix, Class<? extends DifferentiableEmission> insertClass){
		int i=0, a = 0;
		if( createState[i] ) {
			emList.add( SilentEmission.class );
			nameList.add( "D"+(offset) + suffix );
			a++;
		}
		i++;

		if( createState[i] ) {
			emList.add( insertClass );
			nameList.add( "I"+(offset) + suffix );
			a++;
		}
		i++;

		if( createState[i] ) {
			Class<? extends AbstractConditionalDiscreteEmission> de;
			if(!conditionalMain){
				de = DiscreteEmission.class;
			}else{
				de = ReferenceSequenceDiscreteEmission.class;
			}
			emList.add( de );
			nameList.add( "M"+(offset) + suffix );
			a++;
		}

		return a;
	}
	
	/**
	 * Creates the real {@link TransitionElement}s that can be used to create the HMM.
	 * 
	 * @param hyperParams the hyper-parameters for {@link PseudoTransitionElement} from the <code>list</code>
	 * @param list a list of {@link PseudoTransitionElement} that is used to create the HMM
	 * 
	 * @return an array of {@link TransitionElement}s that correspond to all entries of the <code>list</code>
	 * 
	 * @see #propagateESS(double, ArrayList)
	 */
	public static TransitionElement[] createTransition( double[][] hyperParams, ArrayList<PseudoTransitionElement> list ) {
		TransitionElement[] t = new TransitionElement[list.size()];
		for( int i = 0; i < list.size(); i++ ) {
			PseudoTransitionElement current = list.get(i);
			t[i] = new TransitionElement( current.context, current.states, current.prob, current.weights );
		}
		return t;
	}
		
	/**
	 * Propagates the <code>ess</code> for an HMM with absorbing states.
	 * 
	 * @param ess the equivalent sample size of the HMM
	 * @param list a list of {@link PseudoTransitionElement} that is used to create the HMM
	 * 
	 * @return a {@link Pair}, which contains as first element the hyper-parameters for {@link PseudoTransitionElement} from the <code>list</code>, and as second element the ess of each state 
	 */
	public static Pair<double[][], double[]> propagateESS( double ess, ArrayList<PseudoTransitionElement> list ) {
		
		//prepare
		int maxOrder = 0, max = 0, startIdx = -1;
		IntList end = new IntList();
		for( int i = 0; i < list.size(); i++ ) {
			int[] current = list.get(i).states;
			for( int j = 0; j < current.length; j++ ) {
				max = Math.max( max, current[j] );
			}
			if( list.get(i).context.length == 0 ) {
				if( startIdx >= 0 ) {
					throw new IllegalArgumentException( "Multiple start transitions!" );
				} else {
					startIdx = i;
				}
			} else {
				maxOrder = Math.max( maxOrder, list.get(i).context.length );
			}
			
		}
		
		//set children
		for( int t = 0; t < list.size(); t++ ){
			PseudoTransitionElement te = list.get(t);
			if( te.states.length == 0 ) {
				end.add( t );
			}
			for( int child = 0; child < te.states.length; child++){
				int[] nextContext = new int[ Math.min( maxOrder, te.context.length+1 )];
				nextContext[nextContext.length-1] = te.states[child];
				for( int k = 1, j = nextContext.length-2; j >= 0; j--, k++ ) {
					nextContext[j] = te.context[te.context.length - k];
				}
				for( int s = 0; s < list.size(); s++ ){
					PseudoTransitionElement te2 = list.get(s);
					if(nextContext.length == te2.context.length){
						int c=0;
						while( c < nextContext.length && nextContext[c] == te2.context[c] ){
							c++;
						}
						if(c == nextContext.length){
							te.child[child] = s;
							s = list.size();
						}
					}
				}
				
				if( te.child[child] == -1 ) {//nextContext has no descendant
					if( maxOrder == 0 ) {
						te.child[child] = 0;
					} else {
						te.child[child] = list.size();
						list.add( new PseudoTransitionElement( nextContext, null, null ) );
					}
				}
			}
		}
		

		
		if( end.length() == 0 ) {
			throw new IllegalArgumentException( "No absorbing state!" );
		}
		
		//iterate
		double eps = 1E-12, sum = 0;
		
		double[] cumulatedHyper = new double[list.size()];
		double[] currentHyper = new double[list.size()];
		double[] nextHyper = new double[list.size()];
		
		currentHyper[startIdx] = ess;
		
		//int a = 0;
		do {
			//propagate one round
			for( int t = 0; t < list.size(); t++ ){
				PseudoTransitionElement te = list.get(t);
				for( int child = 0; child < te.states.length; child++){
					nextHyper[te.child[child]] += te.prob[child]*currentHyper[t];
				}
			}
			
			//System.out.println("-------------------------------");
			//System.out.println( Arrays.toString(cumulatedHyper) );
			//System.out.println( Arrays.toString(currentHyper) );
			//System.out.println( Arrays.toString(nextHyper) );
			
			//add
			for( int i = 0; i < cumulatedHyper.length; i++ ) {
				cumulatedHyper[i] += currentHyper[i];
				currentHyper[i] = nextHyper[i];
				nextHyper[i] = 0;
			}
			
			//compute remaining hyper-parameters for all states
			sum = 0;
			for( int i = 0; i < currentHyper.length; i++ ) {
				sum += currentHyper[i];
			}
			
			//System.out.println( a++ + "\t" + sum + "\t" + ess ) ;
			
		} while( sum > eps );
		
		//build results
		double[][] t = new double[list.size()][];
		double[] stateESS = new double[max+1];
		for( int i = 0; i < list.size(); i++ ) {
			PseudoTransitionElement current = list.get(i);
			t[i] = new double[current.prob.length];
			
			for( int j = 0; j < t[i].length; j++ ) {
				t[i][j] =  cumulatedHyper[i] * current.prob[j];
			}
			
			if( current.context.length > 0 ) {
				stateESS[current.context[current.context.length-1]] += cumulatedHyper[i];
			}
		}
		
		//System.out.println( Arrays.toString(stateESS) );
				
		return new Pair<double[][], double[]>( t, stateESS );
	}
	
	
	/**
	 * This class is used as place holder for a later {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.BasicHigherOrderTransition.AbstractTransitionElement}. It is used in the factory
	 * and can be used externally for using the method {@link HMMFactory#propagateESS(double, ArrayList)}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class PseudoTransitionElement {
		private int[] context;
		private int[] states;
		private int[] child;
		private double[] prob, weights;

		/**
		 * This constructor creates an new {@link PseudoTransitionElement} without edge weights.
		 * 
		 * @param context the context of this instance; can be <code>null</code>
		 * @param states the states to be visited from the instance; can be <code>null</code>
		 * @param posScore a positive score reflecting the a-priori assumption about the next visited state which is used in the method {@link HMMFactory#propagateESS(double, ArrayList)}; can be <code>null</code>
		 */
		public PseudoTransitionElement( int[] context, int[] states, double[] posScore ) {
			this( context, states, posScore, null );
		}
		
		/**
		 * This constructor creates an new {@link PseudoTransitionElement} with specific edge weights.
		 * 
		 * @param context the context of this instance; can be <code>null</code>
		 * @param states the states to be visited from the instance; can be <code>null</code>
		 * @param posScore a positive score reflecting the a-priori assumption about the next visited state which is used in the method {@link HMMFactory#propagateESS(double, ArrayList)}; can be <code>null</code>
		 * @param weights the edge weight that is used for drawing the final HMM using one oth the methods {@link AbstractHMM#getGraphvizRepresentation(java.text.NumberFormat)}, ...; can be <code>null</code>
		 */
		public PseudoTransitionElement( int[] context, int[] states, double[] posScore, double[] weights ) {
			this.context = context == null ? new int[0] : context.clone();
			this.states = states == null ? new int[0] : states.clone();
			this.prob = new double[this.states.length];
			if( posScore == null ) {
				Arrays.fill( this.prob, 1d / this.states.length );
			} else {
				Normalisation.sumNormalisation( posScore, this.prob, 0 );
			}
			if( weights == null ) {
				this.weights = new double[this.states.length];
				Arrays.fill( this.weights, 1 );
			} else {
				this.weights = weights;
			}
			
			child = new int[this.states.length];
			Arrays.fill( child, -1 );
		}
		
		public String toString(){
			String str = "";
			String cont = "";
			for(int i=0;i<context.length-1;i++){
				cont += context[i]+", ";
			}
			if(context.length > 0){
				cont += context[context.length-1];
			}
			for(int i=0;i<states.length;i++){
				str += "P("+states[i]+"|"+cont+")\t= "+prob[i]+(i==states.length-1 ? "\n" : "\t");
			}
			return str;
		}
	}

	/**
	 * This method returns a {@link HashMap} that can be used in
	 * {@link AbstractHMM#getGraphvizRepresentation(java.text.NumberFormat, de.jstacs.data.DataSet, double[], HashMap)}
	 * to create a Graphviz representation of the {@link AbstractHMM}
	 * 
	 * @return a {@link HashMap} that can be used to create a Graphviz representation
	 * 
	 * @see AbstractHMM#getGraphvizRepresentation(java.text.NumberFormat, de.jstacs.data.DataSet, double[], HashMap)
	 */
	public static HashMap<String,String> getHashMap() {
		HashMap<String, String> hash = new HashMap<String, String>();
		hash.put("J.*", "min" );
		hash.put("[EFIS].*", "same");
		hash.put("M.*", "same");
		hash.put("D.*", "max" );
		return hash;
	}
}
