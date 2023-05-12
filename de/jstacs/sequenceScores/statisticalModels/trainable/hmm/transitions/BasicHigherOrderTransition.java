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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.Storable;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * This class implements the basic transition that allows to be trained using the viterbi or the Baum-Welch algorithm.
 * 
 * Layer 0, context 0 encodes the start state.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class BasicHigherOrderTransition implements TrainableTransition {

	/**
	 * The internally used {@link AbstractTransitionElement}s.
	 */
	protected AbstractTransitionElement[] transitions;
	/**
	 * A vector indicating for each state whether it is silent or not.
	 */
	protected boolean[] isSilent;
	/**
	 * The maximal Markov order of the transition. 
	 */
	protected int maximalMarkovOrder;
	/**
	 * The maximal in-degree of any state.
	 */
	protected int maxInDegree;
	/**
	 * The lookup table for spare context en- and decoding.
	 */
	protected int[][][] lookup;
	
	/**
	 * The main constructor.
	 * 
	 * @param transitions the {@link AbstractTransitionElement}s for the internal use
	 * @param isSilent an array indicating for each state whether it is silent or not
	 * 
	 * @throws Exception if an error occurs during checking the {@link AbstractTransitionElement}s and creating internal fields 
	 */
	public BasicHigherOrderTransition( boolean[] isSilent, AbstractTransitionElement... transitions ) throws Exception {
		this.isSilent = isSilent.clone();
		this.transitions = ArrayHandler.clone( transitions );
		init();
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link BasicHigherOrderTransition} out of an XML representation.
	 * 
	 * @param xml the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link BasicHigherOrderTransition} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 *             
	 */
	public BasicHigherOrderTransition( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getXMLTag() );
		this.transitions = (AbstractTransitionElement[]) XMLParser.extractObjectForTags( xml, "transitions" );
		this.isSilent = (boolean[]) XMLParser.extractObjectForTags( xml, "isSilent" );
		extractFurtherInformation( xml );
		try {
			init();
		} catch (Exception e) {
			NonParsableException npe = new NonParsableException( e.getMessage() );
			throw npe;
		}
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, transitions, "transitions" );
		XMLParser.appendObjectWithTags( xml, isSilent, "isSilent" );
		appendFurtherInformation( xml );
		XMLParser.addTags( xml, getXMLTag() );
		return xml;
	}
	
	private static final String XML_TAG = "BasicHigherOrderTransition";
	
	/**
	 * The method returns the XML tag used during saving and loading the transition.
	 * 
	 * @return he XML tag used during saving and loading the transition.
	 */
	protected String getXMLTag() {
		return XML_TAG;
	}
	
	/**
	 * This method appends further information to the XML representation. It allows subclasses to save further parameters that are not defined in the superclass.
	 * 
	 * @param xml the XML representation
	 */
	protected void appendFurtherInformation( StringBuffer xml ) {
	}

	/**
	 * This method extracts further information from the XML representation. It allows subclasses to cast further parameters that are not defined in the superclass.
	 * 
	 * @param xml the XML representation
	 *  
	 * @throws NonParsableException if the information could not be reconstructed out of the {@link StringBuffer} <code>xml</code>
	 */
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
	}
	
	@Override
	public BasicHigherOrderTransition clone() throws CloneNotSupportedException {
		BasicHigherOrderTransition clone = (BasicHigherOrderTransition) super.clone();
		clone.isSilent = isSilent == null ? null : isSilent.clone();
		clone.lookup = ArrayHandler.clone( lookup );
		clone.transitions = ArrayHandler.clone( transitions );
		return clone;
	}
/*	
	private void showLookup() {
		for( int index = 0; index < lookup[0].length; index++ ) {
		System.out.println(index+" ----------------------------");
			for( int idx, i = 0; i < lookup[0][index].length; i++ ) {
				idx = lookup[0][index][i];
				String s = Arrays.toString( transitions[idx].states );
				s = s.substring( 1, s.length()-1 );
				System.out.println( i + " =^= " + idx + "\t" + Arrays.toString(transitions[idx].context) + " -> {" + s + "}" );
			}
		}
	}
/**/	
	private void addAndTopsortTransitionsWithSilentStates( IntList currentLayer, IntList nextLayer, IntList lookup, boolean [] canBeUsed, boolean[] used, int[] inDeg, boolean[] nextUsed ) {
		int i = 0, j, n, t, idx, d;
		Arrays.fill( canBeUsed, false );
		Arrays.fill( used, false );
		Arrays.fill( inDeg, 0 );
		LinkedList<Integer> roots = new LinkedList<Integer>();
		
		//find all transitions with last silent state
		int thresh = currentLayer.length();
		while( i < currentLayer.length() ) {
			idx = currentLayer.get(i);
			
			n = transitions[idx].getNumberOfChildren();
			for( int s = 0; s < n; s++ ) {
				d = transitions[idx].getDescendant( s );
				if( isSilent[ transitions[idx].getChild( s ) ] ) {
					if( !canBeUsed[d]  ) {						
						currentLayer.add( d );
						canBeUsed[d] = true;
					}
					if( i >= thresh ) {
						inDeg[d]++;
					}
				} else {
					if( !nextUsed[d] ) {
						nextLayer.add( d );
					}
				}
			}
			
			i++;
		}
		
		//topsort
		for( t = 0; t < inDeg.length; t++ ){
			if( canBeUsed[t] ){
				if( inDeg[t] == 0 ){
					roots.add( t );
				}
			}else{
				used[t] = true;
			}
		}
		
		while(roots.size() > 0){
			idx = roots.removeFirst();
			used[idx] = true;
			lookup.add( idx );
			int num = transitions[idx].getNumberOfChildren();
			for(j=0;j<num;j++){
				int child = transitions[idx].getDescendant(j);
				if( canBeUsed[child] ) {
					inDeg[ child ]--;
					if( inDeg[ child ] == 0 ){
						roots.addLast( child );
					}
				}
			}
		}
		
		for( t = 0; t < inDeg.length; t++ ){
			if( canBeUsed[t] && !used[t] ){
				String s = Arrays.toString( transitions[t].states );
				throw new IllegalArgumentException( "Check transition element "+ t + ": " + Arrays.toString( transitions[t].context) + " -> {" + s.substring(1,s.length()-1) + "}" );
			}
		}	
	}
	
	private void createLookup( int layer, IntList current ) {
		lookup[0][layer] = current.toArray();
		lookup[1][layer] = new int[transitions.length];
		Arrays.fill( lookup[1][layer], -1 );
		for( int i = 0; i < lookup[0][layer].length; i++ ) {
			lookup[1][layer][lookup[0][layer][i]] = i;
		}
	}
	
	private void init() throws Exception {		
		computeMaximalMarkovOrder();
		
		//link TransitionElements
		ArrayList<AbstractTransitionElement> elements = new ArrayList<AbstractTransitionElement>( transitions.length*2 );
		for( int t = 0; t < transitions.length; t++ ){
			elements.add( transitions[t] );
		}
		
		AbstractTransitionElement te, te2;
		for( int t = 0; t < elements.size(); t++ ){
			te = elements.get(t);
			for( int child = 0; child <te.getNumberOfChildren(); child++){
				int[] nextContext = te.getNextContext(child, maximalMarkovOrder);
				for( int s = 0; s < elements.size(); s++ ){
					te2 = elements.get(s);
					if(nextContext.length == te2.context.length){
						int c=0;
						while( c < nextContext.length && nextContext[c] == te2.context[c] ){
							c++;
						}
						if(c == nextContext.length){
							te.setIndexOfDescendantTransitionElement(child, s);
						}
					}
				}
				
				if( te.getDescendant(child) == -1 ) {//nextContext has no descendant
					if( maximalMarkovOrder == 0 ) {
						te.setIndexOfDescendantTransitionElement( child, 0 );
					} else {
						te.setIndexOfDescendantTransitionElement( child, elements.size() );
						elements.add( new TransitionElement( nextContext, null, null ) );
					}
				}
			}
		}
		if( elements.size() > transitions.length ) {
			transitions = elements.toArray( new AbstractTransitionElement[0] );
		}
		
		//maxInDegree & check connected & multiple times the same context 
		int[] inDeg = new int[transitions.length];
		for( int t = 0; t < transitions.length; t++ ){
			for( int d = 0; d < transitions[t].descendants.length; d++){
				inDeg[transitions[t].descendants[d]]++;
			}
			for( int s = t+1; s < transitions.length; s++ ){
				if( transitions[t].hasSameContext( transitions[s] ) ) {
					//multiple times the same context 
					throw new IllegalArgumentException( "The context " + Arrays.toString( transitions[t].context ) + " is used by more than one TransitionElements." ); 
				}
			}
		}

		maxInDegree = inDeg[0];
		int minInDegree = Integer.MAX_VALUE, j=0;
		for( int t = 1; t < transitions.length; t++ ){
			maxInDegree = Math.max( maxInDegree, inDeg[t] );
			minInDegree = Math.min( minInDegree, inDeg[t] );
			if( minInDegree > 0 ) {
				j = t;
			}
		}
		if( minInDegree <= 0 ) {
			//not connected
			throw new IllegalArgumentException( "The in-degree of some TransitionElements is zero. " + transitions[j+1] );
		}

		
		lookup = new int[2][maximalMarkovOrder+1][];
		IntList currentLayer = new IntList(), nextLayer = new IntList(), help;
		boolean[] canBeUsed = new boolean[transitions.length];
		boolean[] used = new boolean[transitions.length];
		boolean[] nextUsed = new boolean[transitions.length];
		Arrays.fill( nextUsed, false );
		IntList current = new IntList();
		
		if( maximalMarkovOrder == 0 ) {
			for( int i = 0; i < isSilent.length; i++ ) {
				if( isSilent[i] ) {
					throw new IllegalArgumentException( "A hidden Markov model of order zero is not alllowed to have any silent states." );
				}
			}
			current.add( 0 );
			createLookup( maximalMarkovOrder, current );
		} else {
			//find start TE		
			for( int t=0; t<transitions.length; t++ ){
				//transitions[t].appendGraphvizDescription(sb, null, null);//out
				if( transitions[t].context.length == 0 ) {
					currentLayer.add( t );
					current.add( t );
				}
			}
			//sb.append("}");	System.out.println(sb);//out
			
			//add and topsort contexts with silent states
			addAndTopsortTransitionsWithSilentStates( currentLayer, nextLayer, current, canBeUsed, used, inDeg, nextUsed );	
			createLookup( 0, current );
			help = nextLayer;
			nextLayer = currentLayer;
			currentLayer = help;
			
			//fill context
			for( int idx, n, p = 1; p < maximalMarkovOrder; p++ ) {
				current.clear();
				nextLayer.clear();
				//add all transitions with non-silent state as last state
				Arrays.fill( nextUsed, false);
				
				for( int i = 0; i < currentLayer.length(); i++ )  {			
					idx = currentLayer.get(i);
					current.add( idx );
					
					n = transitions[idx].getNumberOfChildren();
					
					for( int d, s = 0; s < n; s++ ) {
						d = transitions[idx].descendants[s];
						if( !isSilent[ transitions[idx].getChild( s ) ] && !nextUsed[d]) {
							nextLayer.add( d );
							nextUsed[d] = true;
							//System.out.println( "add " + d + ": " + Arrays.toString( transitions[idx].context ) + " -> " + transitions[idx].states[s] );
						}
					}
				}
				
				//add all reachable transitions with a silent state as last state
				addAndTopsortTransitionsWithSilentStates( currentLayer, nextLayer, current, canBeUsed, used, inDeg, nextUsed );
				createLookup( p, current );
				
				//prepare for next iteration
				help = nextLayer;
				nextLayer = currentLayer;
				currentLayer = help;
			}
			
			//full context
			currentLayer.clear();
			nextLayer.clear();
			current.clear();
			//add all transitions with full context and a non-silent state as last state
			Arrays.fill(used, false);
			for(int t=0;t<transitions.length;t++){
				if(transitions[t].context.length == maximalMarkovOrder && !isSilent[ transitions[t].getLastContextState()]){
					currentLayer.add( t );
					current.add( t );
				}
			}
			
			//add all reachable transitions with a silent state as last state
			addAndTopsortTransitionsWithSilentStates( currentLayer, nextLayer, current, canBeUsed, used, inDeg, nextUsed );
			createLookup( maximalMarkovOrder, current );
		}
	}
	
	@Override
	public void resetStatistic() {
		for( int t = 0; t<transitions.length; t++ ){
			transitions[t].resetStatistic();
		}
	}
	
	public void joinStatistics(Transition... transitions){
		AbstractTransitionElement[] ats = new AbstractTransitionElement[transitions.length];
		for(int j=0;j<this.transitions.length;j++){
			for(int i=0;i<transitions.length;i++){
				ats[i] = ((BasicHigherOrderTransition)transitions[i]).transitions[j];
			}
			this.transitions[j].joinStatistics( ats );
		}
	}
	
	@Override
	public void addToStatistic( int layer, int index, int childIdx, double weight, Sequence sequence, int sequencePosition ) {
		transitions[getTransitionElementIndex( layer, index )].addToStatistic( childIdx, weight, sequence, sequencePosition );
	}

	@Override
	public void estimateFromStatistic() {
		for( int t = 0; t<transitions.length; t++ ){
			transitions[t].estimateFromStatistic();
		}
	}

	/**
	 * This method allows to draw parameters from the sufficient statistic, i.e., to draw from the posterior.
	 * If no data has been added using {@link #addToStatistic(int, int, int, double, Sequence, int)}, the method draws from the prior.
	 * 
	 * @throws Exception if for instance the prior has some illegal hyper-parameters (e.g. 0)
	 * 
	 * @see AbstractTransitionElement#drawParametersFromStatistic()
	 */
	public void drawParametersFromStatistic() throws Exception {
		for( int t = 0; t<transitions.length; t++ ){
			transitions[t].drawParametersFromStatistic();
		}
	}
	
	private void computeMaximalMarkovOrder(){
		maximalMarkovOrder = 0;
		int temp = 0;
		for( int t = 0; t<transitions.length; t++ ){
			temp = transitions[t].context.length;
			if(temp > maximalMarkovOrder){maximalMarkovOrder = temp;}
		}
	}
	
	@Override
	public int getMaximalMarkovOrder() {
		return maximalMarkovOrder;
	}

	@Override
	public double getLogPriorTerm() {
		double pt = 0;
		for( int t = 0; t<transitions.length; t++ ){
			pt += transitions[t].getLogPriorTerm();
		}
		return pt;
	}

	@Override
	public int getNumberOfStates() {
		return isSilent.length;
	}

	@Override
	public void initializeRandomly() {
		for( int t = 0; t<transitions.length; t++ ){
			transitions[t].initializeRandomly();
		}
	}

	@Override
	public String getGraphizNetworkRepresentation( NumberFormat nf, String arrowOption, boolean graphical ) {
		StringBuffer res = new StringBuffer();
		for( int t = 0; t < transitions.length; t++ ) {
			transitions[t].appendGraphvizDescription( res, nf, arrowOption, graphical );
		}
		return res.toString();
	}

	@Override
	public void fillTransitionInformation( int layer, int index, int childIdx, int[] container ) {
		int idx = getTransitionElementIndex( layer, index );
		container[0] = transitions[ idx ].states[childIdx];
		container[1] = transitions[ idx ].descendants[childIdx];
		container[2] = isSilent[ container[0] ]?0:1;

		//map back
		container[1] = lookup[1][getLookupIndex(layer+container[2])][container[1]];
	}
	
	//will be inlined
	private int getLookupIndex( int layer ){
		return layer < maximalMarkovOrder ? layer : maximalMarkovOrder;
	}

	//will be inlined
	/**
	 * This method return the index of the {@link AbstractTransitionElement} using the {@link #lookup} table.
	 * 
	 * @param layer the layer of the matrix
	 * @param index the index encoding the context
	 * 
	 * @return the index of the {@link AbstractTransitionElement} with respect to {@link #transitions}
	 */
	protected final int getTransitionElementIndex( int layer, int index){
		return lookup[0][ getLookupIndex( layer ) ][index];
	}

	@Override
	public int getNumberOfChildren( int layer, int index ) {
		return transitions[ getTransitionElementIndex( layer, index ) ].getNumberOfChildren();
	}

	@Override
	public double getLogScoreFor( int layer, int index, int childIdx, Sequence sequence, int sequencePosition ) {
		return transitions[ getTransitionElementIndex( layer, index ) ].getLogScoreFor( childIdx, sequence, sequencePosition );
	}

	@Override
	public int getNumberOfIndexes( int layer ) {
		return lookup[0][ getLookupIndex( layer ) ].length;
	}

	@Override
	public boolean hasAnySelfTransitions() {
		int t = 0, lastState;
		while( t < transitions.length ) {
			if( transitions[t].context.length > 0 ) {
				lastState = transitions[t].getLastContextState();
				for( int c = 0; c < transitions[t].states.length; c++ ) {
					if( transitions[t].states[c] == lastState ) {
						return true;
					}
				}
			}
		}
		return false;
	}

	@Override
	public int getMaximalInDegree() {
		return maxInDegree;
	}
	
	@Override
	public int getMaximalNumberOfChildren() {
		int out = 0; 
		for( int t = 0; t < transitions.length; t++ ) {
			if( out < transitions[t].getNumberOfChildren() ) {
				out = transitions[t].getNumberOfChildren();
			}
		}
		return out;
	}
	
	@Override
	public int getLastContextState( int layer, int index ) {
		int idx = getTransitionElementIndex( layer, index );
		if( transitions[idx].context.length == 0 ) {
			return -1;
		} else {
			return transitions[idx].getLastContextState();
		}
	}
	

	@Override
	public int getChildIdx( int layer, int index, int state ) {
		int idx = getTransitionElementIndex( layer, index ), i = 0;
		while( i < transitions[idx].states.length && state != transitions[idx].states[i] ) {
			i++;
		}
		return i == transitions[idx].states.length ? -1 : i;
	}
	
    public double getLogGammaScoreFromStatistic() {
        double res = 0;
        for( int t = 0; t < transitions.length; t++ ) {
        	res += transitions[t].getLogGammaScoreFromStatistic();
        }
        return res;
    }
    
    public String toString() {
    	return toString( null, null );
    }
    
    public String toString( String[] stateNames, NumberFormat nf ) {
    	StringBuffer sb = new StringBuffer();
    	for( int t = 0; t < transitions.length; t++ ) {
        	sb.append( transitions[t].toString( stateNames, nf ) );
        }
    	return sb.toString();
    }
	
    public boolean[] isAbsorbing() {
    	boolean[] absorbing = new boolean[isSilent.length];
    	Arrays.fill( absorbing, true );
    	for( int t = 0; t < transitions.length; t++ ) {
        	if( transitions[t].context.length > 0 && transitions[t].states.length > 0 ) {
        		absorbing[ transitions[t].getLastContextState() ] = false;
        	}
        }
    	return absorbing;
    }
    
	@Override
	public void setParameters(Transition t) throws IllegalArgumentException {
		if( !t.getClass().equals( getClass() ) ) {
			throw new IllegalArgumentException( "The transitions are not comparable." );
		}
		BasicHigherOrderTransition tt = (BasicHigherOrderTransition) t;
		for( int i = 0; i < transitions.length; i++ ) {
			transitions[i].setParameters( tt.transitions[i] );
		}		
	}
    
	/**
	 * This class declares the probability distribution for a given context, i.e. it contains all possible
	 * transition and the corresponding probabilities for a given set offset previously visited states. 
	 * 
	 * @author Jan Grau, Jens Keilwagen
	 */
	public static abstract class AbstractTransitionElement implements Cloneable, Storable {
		
		/**
		 * The context, i.e. the visited states, of the {@link AbstractTransitionElement}
		 */
		protected int[] context;
		/**
		 * The states that can be visited
		 */
		protected int[] states;
		/**
		 * The hyperparameters of the prior over the parameters.
		 * The order of the elements is the same as for states.
		 * 
		 * @see #getLogPriorTerm()
		 * @see #parameters
		 */
		protected double[] hyperParameters;
		/**
		 * The parameters defining the distribution over all states that can be visited.
		 * The order of the elements is the same as for states.
		 * 
		 * @see #states
		 */
		protected double[] parameters;
		/**
		 * The sufficient statistic for determining the parameters during sampling, viterbi or Baum-Welch training.
		 * The order of the elements is the same as for states.
		 * 
		 * @see #addToStatistic(int, double, Sequence, int)
		 * @see #resetStatistic()
		 * @see #parameters
		 */
		protected double[] statistic;
		/**
		 * The log normalization constant based on the parameters.
		 * 
		 * @see #parameters
		 */
		protected double logNorm;
		/**
		 * The indices for the descendant transition elements that can be visited following the states.
		 * The order of the elements is the same as for states.
		 * 
		 * @see #states
		 */
		protected int[] descendants;
		/**
		 * The weights for plotting the edges with Graphviz.
		 * 
		 * @see #getArrowOption(NumberFormat, double, double, String, boolean)
		 */
		private double[] weight;
		
		/**
		 * This is the main constructor creating a new instance with given context, descendant states, and hyper parameters.
		 * 
		 * @param context the context (=previously visited state indices); last entry corresponds to the last state visited
		 * @param states the transitions to all possible states; if <code>null</code> then no transition allowed
		 * @param hyperParameters the hyper parameters for the transitions; if <code>null</code> then no prior is used
		 * 
		 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.BasicHigherOrderTransition.AbstractTransitionElement#BasicHigherOrderTransition.AbstractTransitionElement(int[], int[], double[], double[])
		 */
		public AbstractTransitionElement( int[] context, int[] states, double[] hyperParameters ){
			this( context, states, hyperParameters, null );
		}

		/**
		 * This is the main constructor creating a new instance with given context, descendant states, and hyper parameters.
		 * 
		 * @param context the context (=previously visited state indices); last entry corresponds to the last state visited
		 * @param states the transitions to all possible states; if <code>null</code> then no transition allowed
		 * @param hyperParameters the hyper parameters for the transitions; if <code>null</code> then no prior is used
		 * @param weight the weight for plotting the edges in Graphviz, enables to modify the edge length, larger weights imply shorter edges (default: 1)
		 */
		public AbstractTransitionElement( int[] context, int[] states, double[] hyperParameters, double[] weight ){
			this.context = context==null ? new int[0] : context.clone();
			int s = states==null ? 0 : states.length, h = hyperParameters == null ? s : hyperParameters.length;
			if( s != h ) {
				throw new IllegalArgumentException( "You have to provide the same number of states and hyperparameters." );
			}
			if( states == null ) {
				this.states = new int[0];
			} else {
				this.states = states.clone();
				Arrays.sort( this.states );
				int old = -1;
				for( int i = 0; i < this.states.length; i++ ) {
					if( this.states[i] == old ) {
						throw new IllegalArgumentException( "It is not allowed to have several edges to the same child. Please check: " + Arrays.toString(context) + " -> " + this.states[i] );
					}
					old = this.states[i];
				}
				System.arraycopy( states, 0, this.states, 0, s );
			}
			this.hyperParameters = new double[s];
			if( hyperParameters != null ) {
				for( int i = 0; i < hyperParameters.length; i++ ) {
					if( hyperParameters[i] < 0 ) {
						throw new IllegalArgumentException( "Please check the hyper-parameter " + i + "." );
					}
					this.hyperParameters[i] = hyperParameters[i];
				}
			}
			this.parameters = new double[s];
			init();
			this.weight = weight==null ? null : weight.clone();
		}
		
		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Constructs a {@link AbstractTransitionElement} out of an XML representation.
		 * 
		 * @param xml the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link AbstractTransitionElement} could not be reconstructed out of
		 *             the {@link StringBuffer} <code>xml</code>
		 *             
		 */
		public AbstractTransitionElement( StringBuffer xml ) throws NonParsableException {
			xml = XMLParser.extractForTag( xml, getXMLTag() );
			this.context = (int[]) XMLParser.extractObjectForTags( xml, "context" );
			this.states = (int[]) XMLParser.extractObjectForTags( xml, "states" );
			this.hyperParameters = (double[]) XMLParser.extractObjectForTags( xml, "hyperparameters" );
			this.parameters = (double[]) XMLParser.extractObjectForTags( xml, "parameters" );
			if( XMLParser.hasTag(xml, "weight", null, null) ) {
				this.weight = (double[]) XMLParser.extractObjectForTags( xml, "weight" );
			} else {
				this.weight = null;
			}
			extractFurtherInformation( xml );
			init();
		}
		
		/**
		 * This method appends further information to the XML representation.
		 * It allows subclasses to save further parameters that are not defined in the superclass.
		 * 
		 * @param xml the XML representation
		 */
		protected abstract void appendFurtherInformation( StringBuffer xml );

		/**
		 * This method extracts further information from the XML representation.
		 * It allows subclasses to cast further parameters that are not defined in the superclass.
		 * 
		 * @param xml the XML representation
		 *  
		 * @throws NonParsableException if the information could not be reconstructed out of the {@link StringBuffer} <code>xml</code>
		 */
		protected abstract void extractFurtherInformation( StringBuffer xml ) throws NonParsableException;
		
		/**
		 * This method returns the xml tag used in {@link #toXML()}.
		 * 
		 * @return the xml tag used in {@link #toXML()}
		 */
		protected String getXMLTag() {
			return XML_TAG;
		}
		
		private static final String XML_TAG = "TRANSITION_ELEMENT";
		
		@Override
		public StringBuffer toXML() {
			StringBuffer xml = new StringBuffer();
			XMLParser.appendObjectWithTags( xml, context, "context" );
			XMLParser.appendObjectWithTags( xml, states, "states" );
			XMLParser.appendObjectWithTags( xml, hyperParameters, "hyperparameters" );
			XMLParser.appendObjectWithTags( xml, parameters, "parameters" );
			XMLParser.appendObjectWithTags( xml, weight, "weight" );
			appendFurtherInformation( xml );
			XMLParser.addTags( xml, getXMLTag() );
			return xml;
		}
		
		/**
		 * This method initializes internal fields.
		 */
		protected void init() {
			this.statistic = new double[states.length];
			this.descendants = new int[states.length];
			Arrays.fill( descendants, -1 );
			precompute();
		}

		public AbstractTransitionElement clone() throws CloneNotSupportedException {
			AbstractTransitionElement clone = (AbstractTransitionElement) super.clone();
			clone.context = context.clone();
			clone.states = states.clone();
			clone.hyperParameters = hyperParameters.clone();
			clone.parameters = parameters.clone();
			clone.statistic = statistic.clone();
			clone.descendants = descendants.clone();
			clone.weight = weight==null ? null : weight.clone();
			return clone;
		}
		
		/**
		 * Returns the last state of the context
		 * @return the last state
		 */
		public int getLastContextState() {
			return context[context.length-1];
		}

		/**
		 * This method precomputes internal fields as for instance the normalization constant.
		 * 
		 * @see #logNorm
		 */
		protected void precompute() {
			logNorm = Normalisation.getLogSum( parameters );
		}
		
		/**
		 * This method appends the current transition element to a <i>Graphviz</i> representation of the structure that
		 * can be used to create an image.
		 * 
		 * @param representation the current <i>Graphviz</i> representation
		 * @param nf the {@link NumberFormat} used for the probabilities, if <code>null</code> no probabilities will we written
		 * @param arrowOption this parameter gives the possibility to set some arrow option 
		 * @param graphical represent transition probabilities as thickness of edges instead of textual output
		 * 
		 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.Transition#getGraphizNetworkRepresentation(NumberFormat, String, boolean)
		 * @see #appendTransitions(StringBuffer, String, NumberFormat, String, boolean)
		 */
		public void appendGraphvizDescription( StringBuffer representation, NumberFormat nf, String arrowOption, boolean graphical ) {
			if( states.length > 0 ) {
				String contextNode;
				switch( context.length ) {
					case 0: contextNode = AbstractHMM.START_NODE; break;
					case 1: contextNode = "" + context[0];break;
					default:
						contextNode = "" + context[0];
						for( int c = 1; c < context.length; c++ ) {
							contextNode += "_" + context[c];
						}
						representation.append( "\tp"+contextNode+"[label=\"\",fixedsize=true,width=0,height=0]\n" );
						
						representation.append( "\t{" +contextNode.replace('_', ' ') + "} -> p" + contextNode + "[arrowhead=none,style=dashed]\n" );
						/*//TODO improve layout for hyper-edges					
						for( int c = 0; c < context.length; c++ ) {
							if( c == context.length-1 ) {
								representation.append( "\t" +context[c] + " -> p" + contextNode + "[arrowhead=none,style=dashed]\n" );
							} else {
								representation.append( "\t" +context[c] + " -> p" + contextNode + "[arrowhead=none,style=dashed,constraint=false]\n" );
							}
						}*/
						contextNode = "p" + contextNode;
				}
				appendTransitions( representation, contextNode, nf, arrowOption, graphical );
			}
		}
		
		/**
		 * This method appends all edges of the transition element to a given <i>Graphviz</i> representation.
		 * 
		 * @param representation the current <i>Graphviz</i> representation
		 * @param contextNodeRepresentation the String representation of the context node(s)
		 * @param nf the {@link NumberFormat} used for the probabilities, if <code>null</code> no probabilities will we written
		 * @param arrowOption this parameter gives the possibility to set some arrow option 
		 * @param graphical represent transition probabilities as thickness of edges instead of textual output
		 * 
		 * @see #appendGraphvizDescription(StringBuffer, NumberFormat, String, boolean)
		 * @see #getArrowOption(NumberFormat, double, double, String, boolean)
		 */
		protected void appendTransitions( StringBuffer representation, String contextNodeRepresentation, NumberFormat nf, String arrowOption, boolean graphical  ) {
			for( int s = 0; s < states.length; s++ ) {
				representation.append( "\t" + contextNodeRepresentation + "->" + states[s] + getArrowOption( nf, Math.exp(parameters[s]-logNorm), getGraphvizEdgeWeight(s), arrowOption, graphical ) + "\n" );	
			}
		}
		
		/**
		 * This method returns the edge weight for plotting the edge with Graphviz.
		 * Larger weights imply shorter edges.
		 * 
		 * @param s the index of the states (child)
		 * 
		 * @return the weight for the edge
		 */
		protected final double getGraphvizEdgeWeight( int s ) {
			return weight == null ? 1 : weight[s];
		}
		
		/**
		 * This method returns the option for an edge in Graphviz. Given a probability and a format
		 * it returns a {@link String} that is used in the Graphviz notation to quantify a specific transition.
		 * 
		 * @param nf the {@link NumberFormat} for formating the probability
		 * @param prob the probability of a certain transition
		 * @param weight the weight for plotting the edges in Graphviz, enables to modify the edge length, larger weights imply shorter edges (default: 1)
		 * @param arrowOption further arrow options (e.g. color, ...)
		 * @param graphical represent transition probabilities as thickness of edges instead of textual output
		 * 
		 * @return a {@link String} that is used in the Graphviz notation to quantify a specific transition
		 * 
		 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.Transition#getGraphizNetworkRepresentation(NumberFormat, String, boolean)
		 * @see #appendTransitions(StringBuffer, String, NumberFormat, String, boolean)
		 */
		protected static String getArrowOption( NumberFormat nf, double prob, double weight, String arrowOption, boolean graphical ) {
			if( graphical ) {
				double penwidth = 1.0 + prob*20.0;
				return "[penwidth=\""+penwidth+"\", weight=\""+weight+"\"]";
			} else if( nf == null && arrowOption == null ) {
				return "[weight=\""+weight+"\"]";
			} else if( nf == null ){
				return "[" + arrowOption + ", weight=\""+weight+"\"]";
			} else {
				return "[label=" + nf.format( prob ) + (arrowOption==null?"":", " + arrowOption) + ", weight=\""+weight+"\"]";
			}
			/*
			String res = "[weight=\""+weight+"\"";
			if( nf != null ) {
				res += ", label=" + nf.format( prob );
			}
			if( graphical ) {
				double penwidth = 1.0 + prob*20.0;
				res += ", penwidth=\""+penwidth+"\""; 
			}
			if( arrowOption != null ) {
				res += ", " + arrowOption;
			}
			return res + "]";
			*/
		}
		
		/**
		 * This method returns the state index encoded by the child index.
		 * 
		 * @param index the index of the child
		 * 
		 * @return the state index
		 * 
		 * @see #states
		 */
		public int getChild( int index ) {
			return states[index];
		}
		
		/**
		 * This method returns the index of the descendant transition element when following the child with index <code>index</code>
		 * 
		 * @param index the child index
		 * 
		 * @return the index of the descendant transition element
		 * 
		 * @see #descendants
		 * @see #setIndexOfDescendantTransitionElement(int, int)
		 */
		public int getDescendant( int index ) {
			return descendants[index];
		}

		/**
		 * This method sets the index of the descendant transition element for the child with index <code>index</code>. 
		 * 
		 * @param index the index of the child
		 * @param descendant the index of the descendant transition element
		 * 
		 * @see #descendants
		 * @see #getDescendant(int)
		 */
		public void setIndexOfDescendantTransitionElement( int index, int descendant ) {
			descendants[index] = descendant;
		}

		/**
		 * This method returns the next context that will be visited when visiting the child with index <code>index</code>.
		 * 
		 * @param index the index of the child
		 * @param maximalMarkovOrder the maximal Markov order to be used
		 * 
		 * @return the next context that will be visited
		 */
		public int[] getNextContext( int index, int maximalMarkovOrder ) {
			int myOrd = context.length < maximalMarkovOrder ? context.length + 1 : maximalMarkovOrder;
			int[] newContext = new int[myOrd];
			if( myOrd > 0 ) {
				System.arraycopy( context, context.length < maximalMarkovOrder ? 0 : 1, newContext, 0, myOrd-1 );
				newContext[myOrd-1] = states[index];
			}
			return newContext;
		}

		/**
		 * This method returns the score for the transition from the current context to the state with index <code>index</code>.
		 * 
		 * @param index the index of the child
		 * @param sequence the sequence, which might be used to retrieve {@link de.jstacs.data.sequences.annotation.SequenceAnnotation}
		 * @param sequencePosition the position within the sequence
		 * 
		 * @return the score for the transition from the current context to the state with index <code>index</code>
		 */
		public double getLogScoreFor( int index, Sequence sequence, int sequencePosition ) {
			return parameters[index] - logNorm;
		}
		
		/**
		 * Returns a value that is proportional to the log of the prior. For maximum
		 * likelihood (ML) 0 should be returned.
		 * 
		 * @return a value that is proportional to the log of the prior
		 * 
		 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.Transition#getLogPriorTerm()
		 */
		public double getLogPriorTerm() {
			if( hyperParameters.length > 0 ) {
				double res = 0, sumHyper = 0;
				for( int d = 0; d < hyperParameters.length; d++ ) {
					sumHyper += hyperParameters[d];
					res += hyperParameters[d] * parameters[d];
				}
				if( sumHyper == 0 ) {
					return 0;
				} else {
					return res - sumHyper * logNorm;
				}
			} else {
				return 0;
			}
		}

		/**
		 * This method returns the number of states that can be visited.
		 *  
		 * @return the number of states that can be visited
		 * 
		 * @see #states
		 */
		public final int getNumberOfChildren() {
			return states.length;
		}

		/**
		 * This method adds a given weight to the sufficient statistic for the parameters.
		 * 
		 * @param childIdx the index of the descendant state
		 * @param weight the weight to be added
		 * @param sequence the sequence, which might be used to retrieve {@link de.jstacs.data.sequences.annotation.SequenceAnnotation}
		 * @param sequencePosition the position within the sequence
		 * 
		 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.TrainableTransition#addToStatistic(int, int, int, double, Sequence, int)
		 * @see #statistic
		 * @see #resetStatistic()
		 * @see #estimateFromStatistic()
		 * @see #drawParametersFromStatistic()
		 */
		public void addToStatistic( int childIdx, double weight, Sequence sequence, int sequencePosition ) {
			statistic[ childIdx ] += weight;
		}
		
		/**
		 * This method joins the statistics of different instances and sets this joined statistic as statistic of each instance.
		 * 
		 * This method might be used for instance in a multi-threaded optimization to join partial statistics.
		 * 
		 * @param te the transition elements to be joined
		 */
		public void joinStatistics(BasicHigherOrderTransition.AbstractTransitionElement... te){
			for(int j=0;j<te.length;j++){
				if(te[j] != this){
					for(int i=0;i<statistic.length;i++){
						statistic[i] += te[j].statistic[i];
					}
				}
			}
			for(int j=0;j<te.length;j++){
				if(te[j] != this){
					System.arraycopy( this.statistic, 0, te[j].statistic, 0, this.statistic.length );
				}
			}
		}

		/**
		 * This method returns the number of parameters in this transition element.
		 * 
		 * @return the number of parameters in this transition element
		 * 
		 * @see #parameters
		 */
		public int getNumberOfParameters() {
			return parameters.length;
		}

		/**
		 * This method draws new parameters from the sufficient statistics.
		 * 
		 * @see #statistic
		 * @see #resetStatistic()
		 * @see #addToStatistic(int, double, Sequence, int)
		 */
		public void drawParametersFromStatistic() {
			if( states.length > 1 ) {
				double sum = 0;
				for( int i = 0; i < statistic.length; i++ ) {
					parameters[i] = statistic[i] + hyperParameters[i];
					sum += parameters[i];
				}
				DirichletMRGParams par;
				if( sum == 0 ) {
					par = new DirichletMRGParams( 1, states.length );
				} else {
					par = new DirichletMRGParams( parameters );
				}
				DirichletMRG.DEFAULT_INSTANCE.generateLog(parameters, 0, parameters.length, par);
			} else {
				Arrays.fill( parameters, 0 );
			}
			precompute();
		}

		/**
		 * This method estimates the parameters from the sufficient statistic.
		 * 
		 * @see #statistic
		 * @see #resetStatistic()
		 * @see #addToStatistic(int, double, Sequence, int)
		 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.TrainableTransition#estimateFromStatistic()
		 */
		public void estimateFromStatistic() {
			if( states.length > 0 ) {
				logNorm = 0;
				for(int i=0;i<statistic.length;i++){
					statistic[i] += hyperParameters[i];
					logNorm += statistic[i];
					parameters[i] = Math.log( statistic[i] );
				}
				if( logNorm == 0 ) {
					Arrays.fill( parameters, 0 );
					logNorm = Math.log( parameters.length );
				} else {
					logNorm = Math.log( logNorm );
				}
			}
		}

		/**
		 * This method resets the sufficient statistic for the parameters.
		 * 
		 * @see #statistic
		 * @see #estimateFromStatistic()
		 * @see #drawParametersFromStatistic()
		 * @see #addToStatistic(int, double, Sequence, int)
		 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.TransitionWithSufficientStatistic#resetStatistic()
		 */
		public void resetStatistic() {
			Arrays.fill( statistic, 0 );
		}

		/**
		 * This method draws new parameters from the prior.
		 */
		public void initializeRandomly() {
			resetStatistic();
			drawParametersFromStatistic();
			resetStatistic();
		}
		
		/**
		* This method calculates a score for the current statistics, which is independent from the current parameters
		*
		* In general the gamma-score is a product of gamma-functions parameterized with the current statistics
		*
		* @return the logarithm of the gamma-score for the current statistics
		*/
		public double getLogGammaScoreFromStatistic() {
			double sum = 0, all = 0, res = 0;
			for( int state = 0; state < hyperParameters.length; state++ ) {
				sum += hyperParameters[state];
				all += statistic[state];
				res += Gamma.logOfGamma( statistic[state] ) - Gamma.logOfGamma( hyperParameters[state] );
			}
			res += Gamma.logOfGamma( sum ) - Gamma.logOfGamma( all );
	        return res;
		}
		
		public final String toString() {
			return toString( null, null );
		}
		
		/**
		 * This method returns a {@link String} representation of the transition element using the given names of the states.
		 * 
		 * @param stateNames the names of the states, can be <code>null</code>
		 * @param nf the {@link NumberFormat} for the {@link String} representation of probabilities
		 * 
		 * @return a {@link String} representation of the transition element using the given names of the states
		 */
		public String toString( String[] stateNames, NumberFormat nf ) {
			//return Arrays.toString( context ) + " -> " + Arrays.toString( states ) + " (desc: " + Arrays.toString( descendants ) + ")";
			if( parameters.length > 0 ) {
				StringBuffer sb = new StringBuffer();
				String context = getContext( stateNames );
				for( int i = 0; i < parameters.length; i++ ) {
					double v = Math.exp( parameters[i] - logNorm );
					sb.append("P(" + getLabel( stateNames, states[i] ) + context + ") \t= " + (nf==null?v:nf.format( v )) );
					sb.append("\t");
				}
				sb.append( "\n" );
				return sb.toString();
			} else {
				return "";
			}
		}
		
		/**
		 * This method returns a {@link String} representation of the context.
		 * 
		 * @param stateNames the names of the states, can be <code>null</code>
		 * 
		 * @return a {@link String} representation of the context
		 * 
		 * @see #toString()
		 * @see #toString(String[],NumberFormat)
		 */
		protected String getContext( String[] stateNames ) {
			if( context.length==0 ) {
				return "";
			} else {
				String c = "|" + getLabel( stateNames, context[0] );
				for( int i = 1; i < context.length; i++ ) {
					c += ", " + getLabel( stateNames, context[i] );
				}
				return c;
			}
		}
		
		/**
		 * This method returns a label for the state.
		 * 
		 * @param stateNames the names of the states, can be <code>null</code>
		 * @param stateIdx the index of the states
		 * @return the label for the state
		 * 
		 * @see #toString(String[], NumberFormat)
		 * @see #getContext(String[])
		 */
		protected final String getLabel( String[] stateNames, int stateIdx ) {
			return stateNames==null?""+stateIdx:stateNames[stateIdx];
		}
		
		private boolean hasSameContext( AbstractTransitionElement te ) {
			if( te.context.length == context.length ) {
				int i = 0;
				while( i < context.length && context[i] == te.context[i] ) {
					i++;
				}
				return i == context.length;
			} else {
				return false;
			}
		}
		
		/**
		 * Set values of parameters of the instance to the value of the parameters of the given instance.
		 * It can be assumed that the given instance and the current instance are from the same class.
		 * 
		 * This method might be used for instance in a multi-threaded optimization to broadcast the parameters. 
		 * 
		 * @param t the transition element with the parameters to be set 
		 * 
		 * @throws IllegalArgumentException if the assumption about the same class for given and current instance is wrong
		 */
		public void setParameters( AbstractTransitionElement t ) throws IllegalArgumentException {
			if( !t.getClass().equals( getClass() ) || t.parameters.length != parameters.length ) {
				throw new IllegalArgumentException( "The transition elements are not comparable." );
			}
			System.arraycopy( t.parameters, 0, parameters, 0, t.parameters.length );
			precompute();
		}
	}
}