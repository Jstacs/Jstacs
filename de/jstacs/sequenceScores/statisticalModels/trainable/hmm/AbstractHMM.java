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

package de.jstacs.sequenceScores.statisticalModels.trainable.hmm;

import java.io.OutputStream;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import de.jstacs.Storable;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.AbstractTrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.State;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.HMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.MultiThreadedTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.BasicHigherOrderTransition;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.HigherOrderTransition;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.Transition;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.BasicHigherOrderTransition.AbstractTransitionElement;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.ToolBox;

/**
 * This class is the super class of all implementations hidden Markov models (HMMs) in Jstacs.
 * The training algorithm of the the HMM is determined by a specialized {@link de.jstacs.parameters.ParameterSet}
 * denoted as {@link HMMTrainingParameterSet}.
 * For creating frequently used HMMs please check {@link HMMFactory}.
 * 
 * @author Jan Grau, Jens Keilwagen, Michael Scharfe
 * 
 * @see State
 * @see Transition
 * @see HMMFactory
 */
public abstract class AbstractHMM extends AbstractTrainableStatisticalModel implements Cloneable, Storable {

	/**
	 * The (hidden) states of the HMM.
	 */
	protected State[] states;

	//for the states
	/**
	 * The names of the states.
	 */
	protected String[] name;
	/**
	 * The index of the used emission of each state.
	 */
	protected int[] emissionIdx;
	/**
	 * An array of switches that contains for each state whether the emission is forward or the reverse strand.
	 * Only of interest for discrete, reverse complementable sequences.
	 * 
	 * @see de.jstacs.data.alphabets.ComplementableDiscreteAlphabet
	 */
	protected boolean[] forward;
	/**
	 * The emissions used in the states.
	 */
	protected Emission[] emission;
	
	/**
	 * The transitions between all (hidden) states of the HMM.
	 */
	protected Transition transition;
	
	/**
	 * matrix for all forward-computed variables;
	 * fwdMatrix[l][c] = log P(x_1,...,x_l,(s_{l-order+1},...,s_l)=c | parameter)
	 */
	protected double[][] fwdMatrix;
	
	/**
	 * matrix for all backward-computed variables;
	 * bwdMatrix[l][c] = log P(x_{l+1},...,x_L | (s_{l-order+1},...,s_l)=c , parameter)
	 */
	protected double[][] bwdMatrix;
	
	/**
	 * The {@link de.jstacs.parameters.ParameterSet} containing all {@link de.jstacs.parameters.Parameter}s for the training of the HMM.
	 */
	protected HMMTrainingParameterSet trainingParameter;
	
	/**
	 * This is the stream for writing information while training.
	 */
	protected SafeOutputStream sostream;
	
	/**
	 * An array of switches that contains for each state whether is is a final state or not (cf. <a href="#final">final states</a>).
	 */
	protected boolean[] finalState;

	/**
	 * The number of threads that is internally used.
	 */
	protected int threads;
	
	/**
	 * This is the main constructor for an HMM.
	 * 
	 * @param trainingParameterSet a {@link de.jstacs.parameters.ParameterSet} containing all {@link de.jstacs.parameters.Parameter}s for the training of the HMM
	 * @param name the names of the states
	 * @param emissionIdx the indices of the emissions that should be used for each state, if <code>null</code> state <code>i</code> will use emission <code>i</code>
	 * @param forward a boolean array that indicates whether the symbol on the forward or the reverse complementary strand should be used,
	 * 				  if <code>null</code> all states use the forward strand
	 * @param emission the emissions
	 * 
	 * @throws CloneNotSupportedException if <code>trainingParameterSet</code> can not be cloned
	 * @throws WrongAlphabetException if not all (non-silent) emissions have use the same {@link AlphabetContainer}
	 */
	protected AbstractHMM( HMMTrainingParameterSet trainingParameterSet, String[] name, int[] emissionIdx, boolean[] forward, Emission[] emission ) throws CloneNotSupportedException, WrongAlphabetException {
		super( getAlphabetContainer( emission ), 0 );
		if( !trainingParameterSet.hasDefaultOrIsSet() ) {
			throw new IllegalArgumentException( "Please check the training parameters." );
		}
		this.trainingParameter = (HMMTrainingParameterSet) trainingParameterSet.clone();
		setThreads();
		setOutputStream( SafeOutputStream.DEFAULT_STREAM );
		
		int n = name.length;
		this.name = new String[n];
		HashSet<String> hash = new HashSet<String>();
		for( int i = 0; i < n; i++ ) {
			if( hash.contains( name[i] ) ) {
				throw new IllegalArgumentException( "The state names should be unique. Please check: " + name[i] );
			}
			this.name[i] = name[i];
			hash.add( name[i] );
		}
		hash.clear();
		hash = null;
		
		if( emissionIdx == null ) {
			this.emissionIdx = new int[n];
			for( int e = 0; e < n; e++ ) {
				this.emissionIdx[e] = e;
			}
		} else {
			if( n != emissionIdx.length ) {
				throw new IllegalArgumentException();
			}
			this.emissionIdx = emissionIdx.clone();
		}
		if( forward == null ) {
			this.forward = new boolean[n];
			Arrays.fill( this.forward, true );
		} else {
			if( n != forward.length ) {
				throw new IllegalArgumentException();
			}
			this.forward = forward.clone();
		}
		
		if( emission.length > n ) {
			throw new IllegalArgumentException();
		}
		
		this.emission = ArrayHandler.clone(emission);
	}
	
	private void setThreads() {
		threads = trainingParameter instanceof MultiThreadedTrainingParameterSet ? ((MultiThreadedTrainingParameterSet) trainingParameter).getNumberOfThreads() : 1;
	}
	
	private static AlphabetContainer getAlphabetContainer( Emission... e ) throws WrongAlphabetException {
		AlphabetContainer con = null;
		int i = 0;
		while( con == null ) {
			con = e[i++].getAlphabetContainer();
		}
		while( i < e.length ) {
			AlphabetContainer current = e[i++].getAlphabetContainer();
			if( current != null && !current.checkConsistency( con ) ) {
				throw new WrongAlphabetException( "All emission should use the same AlphabetContainer." );
			}
		}
		if( !con.isSimple() ) {
			throw new IllegalArgumentException( "The AlphabetContainer has to be simple." );
		}
		return con;
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link AbstractHMM} out of an XML representation.
	 * 
	 * @param xml the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AbstractHMM} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	protected AbstractHMM( StringBuffer xml ) throws NonParsableException {
		super( xml );
		setOutputStream( SafeOutputStream.DEFAULT_STREAM );
	}
	
	/**
	 * This method creates the internal transition. 
	 * 
	 * @param te the individual transition elements
	 * 
	 * @throws Exception if the transition can not handle the current states
	 */
	protected void initTransition( AbstractTransitionElement... te ) throws Exception {
		boolean[] isSilent = new boolean[states.length];
		for( int i = 0; i < states.length; i++ ) {
			isSilent[i] = states[i].isSilent();
		}
		if( te instanceof TransitionElement[] ) {
			transition = new HigherOrderTransition( isSilent, (TransitionElement[]) te );
		} else {
			int t = 0;
			TransitionElement[] help = new TransitionElement[te.length];
			for( int i = 0; i < help.length; i++ ) {
				if( te[t] instanceof TransitionElement ) {
					help[t] = (TransitionElement) te[t];
				} else {
					break;
				}
			}
			if( t == te.length ) {
				transition = new HigherOrderTransition( isSilent, help );
			} else {
				transition = new BasicHigherOrderTransition( isSilent, te );
			}
		}
	}
	
	/**
	 * Returns the tag for the XML representation.
	 * 
	 * @return the tag for the XML representation
	 * 
	 * @see AbstractHMM#fromXML(StringBuffer)
	 * @see AbstractHMM#toXML()
	 */
	protected abstract String getXMLTag();
	
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, trainingParameter, "trainingParameter" );
		XMLParser.appendObjectWithTags( xml, transition, "transition" );
		
		XMLParser.appendObjectWithTags( xml, name, "name" );
		XMLParser.appendObjectWithTags( xml, emissionIdx, "emissionIdx" );
		XMLParser.appendObjectWithTags( xml, forward, "strand" );
		XMLParser.appendObjectWithTags( xml, emission, "emission" );
		
		appendFurtherInformation( xml );
		XMLParser.addTags( xml, getXMLTag() );
		return xml;
	}
	
	/**
	 * This method is used by the {@link AbstractHMM#AbstractHMM(StringBuffer)} constructor for creating an instance from an XML representation.
	 * This method should never be made <code>public</code>.
	 * 
	 * @param xml the XML representation
	 * @throws NonParsableException if the XML representation can not be parsed properly
	 */
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		length =0;
		xml = XMLParser.extractForTag( xml, getXMLTag() );
		trainingParameter = (HMMTrainingParameterSet) XMLParser.extractObjectForTags( xml, "trainingParameter" );
		setThreads();
		
		transition = ((Transition) XMLParser.extractObjectForTags( xml, "transition" ));
		
		name = XMLParser.extractObjectForTags( xml, "name", String[].class );
		emissionIdx = XMLParser.extractObjectForTags( xml, "emissionIdx", int[].class );
		forward = XMLParser.extractObjectForTags( xml, "strand", boolean[].class );
		emission = (XMLParser.extractObjectForTags( xml, "emission", Emission[].class ));
		
		extractFurtherInformation( xml );
		
		try {
			alphabets = getAlphabetContainer( emission );
		} catch (WrongAlphabetException e) {
			NonParsableException npe = new NonParsableException( e.getMessage() );
			throw npe;
		}
		createStates();
		determineFinalStates();
	}
	
	/**
	 * This method appends further information to the XML representation. It allows subclasses to save further parameters that are not defined in the superclass.
	 * 
	 * @param xml the XML representation
	 */
	protected abstract void appendFurtherInformation( StringBuffer xml );
	
	/**
	 * This method extracts further information from the XML representation. It allows subclasses to cast further parameters that are not defined in the superclass.
	 * 
	 * @param xml the XML representation
	 *  
	 * @throws NonParsableException if the information could not be reconstructed out of the {@link StringBuffer} <code>xml</code>
	 */
	protected abstract void extractFurtherInformation( StringBuffer xml ) throws NonParsableException;
	
	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	//XXX states still have to be *cloned/recreated* => createStates() after cloning emissions;
	public AbstractHMM clone() throws CloneNotSupportedException {
		AbstractHMM clone = (AbstractHMM) super.clone();
		clone.name = name.clone();
		clone.emissionIdx = emissionIdx.clone();
		clone.forward = forward.clone();
		clone.emission = ArrayHandler.clone( emission );
		clone.transition = transition.clone();
		clone.fwdMatrix = ArrayHandler.clone( fwdMatrix );
		clone.bwdMatrix = ArrayHandler.clone( bwdMatrix );
		clone.trainingParameter = (HMMTrainingParameterSet) trainingParameter.clone();
		clone.finalState = finalState.clone();
		
		clone.createStates();		
		clone.setOutputStream( sostream.doesNothing() ? null : SafeOutputStream.DEFAULT_STREAM );
		return clone;
	}
	
	/**
	 * This method creates states for the internal usage.
	 */
	protected abstract void createStates();
	
	/**
	 * This method fills the forward-matrix for a given sequence.
	 * 
	 * @param startPos the start position (inclusive) in the sequence
	 * @param endPos the end position (inclusive) in the sequence
	 * @param seq the sequence
	 * 
	 * @throws Exception if some error occurs during the computation 
	 */
	protected abstract void fillFwdMatrix( int startPos, int endPos, Sequence seq ) throws Exception;
	
	/**
	 * This method fills the backward-matrix for a given sequence.
	 * 
	 * @param startPos the start position (inclusive) in the sequence
	 * @param endPos the end position (inclusive) in the sequence
	 * @param seq the sequence
	 * 
	 * @throws Exception if some error occurs during the computation
	 */
	protected abstract void fillBwdMatrix( int startPos, int endPos, Sequence seq ) throws Exception;
	
	/**
	 * The {@link String} for the start node used in Graphviz annotation.
	 * 
	 * @see #getGraphvizRepresentation(NumberFormat)
	 */
	public static final String START_NODE = "START";
	
	/**
	 * This method returns the number of threads that is internally used.
	 * 
	 * @return the number of threads that is internally used
	 */
	public int getNumberOfThreads(){
		return threads;
	}
	
	/**
	 * This method returns a {@link String} representation of the structure that
	 * can be used in <i>Graphviz</i> to create an image.
	 * 
	 * @param nf an instance of {@link NumberFormat} for formating the probabilities of the transition
	 * 
	 * @return a {@link String} representation of the structure
	 * 
	 * @see #getGraphvizRepresentation(NumberFormat, DataSet, double[], boolean)
	 */
	public String getGraphvizRepresentation( NumberFormat nf ) {
		return getGraphvizRepresentation( nf, null, null, false );
	}
	
	/**
	 * This method returns a {@link String} representation of the structure that
	 * can be used in <i>Graphviz</i> to create an image.
	 * 
	 * @param nf an instance of {@link NumberFormat} for formating the probabilities of the transition
	 * @param sameTypeSameRank if <code>true</code>, states of the same type, i.e., having the same type of emission, are displayed on the same rank
	 * 
	 * @return a {@link String} representation of the structure
	 * 
	 * @see #getGraphvizRepresentation(NumberFormat, DataSet, double[], boolean)
	 */
	public String getGraphvizRepresentation( NumberFormat nf, boolean sameTypeSameRank ) {
		return getGraphvizRepresentation( nf, null, null, sameTypeSameRank );
	}
	
	/**
	 * This method returns a {@link String} representation of the structure that
	 * can be used in <i>Graphviz</i> to create an image.
	 * 
	 * @param nf an instance of {@link NumberFormat} for formating the probabilities of the transition
	 * @param data the data to determine the state posterior; can be <code>null</code> 
	 * @param weight the weights to weight the determined state posterior; can be <code>null</code>
	 * @param sameTypeSameRank if <code>true</code>, states of the same type, i.e., having the same type of emission, are displayed on the same rank
	 * 
	 * @return a {@link String} representation of the structure
	 */
	public String getGraphvizRepresentation( NumberFormat nf, DataSet data, double[] weight, boolean sameTypeSameRank ) {
		HashMap<String, String> map = new HashMap<String, String>();
		String regex = ".*";
		for( int i = 0; i < name.length; i++ ) {
			map.put( name[i].charAt(0) + regex, "same" );
		}
		return getGraphvizRepresentation(nf, data, weight, map );
	}
	
	/**
	 * This method returns a {@link String} representation of the structure that
	 * can be used in <i>Graphviz</i> to create an image.
	 * 
	 * @param nf an instance of {@link NumberFormat} for formating the probabilities of the transition
	 * @param data the data to determine the state posterior; can be <code>null</code> 
	 * @param weight the weights to weight the determined state posterior; can be <code>null</code>
	 * @param rankPatterns a {@link HashMap} contain regular expressions and their corresponding value for the option <code>rank</code> in Graphviz
	 * 
	 * @return a {@link String} representation of the structure
	 * 
	 * @see HMMFactory#getHashMap()
	 */
	public String getGraphvizRepresentation( NumberFormat nf, DataSet data, double[] weight, HashMap<String, String> rankPatterns ) {
		double[] freq = null;
		double maxFreq = 0;
		if(data != null){
			try {
				freq = getStateFreq( data, weight );
				maxFreq = ToolBox.max( freq );
			} catch (Exception e) {
				e.printStackTrace();
				freq = new double[states.length];
				maxFreq = 0;
			}
		}else{
			freq = new double[states.length];
			Arrays.fill( freq, -1 );
			maxFreq = -1;
		}
		StringBuffer sb = new StringBuffer();
		sb.append( "digraph G {\n\trankdir="+(rankPatterns != null ? "TB" : "LR")+"\n\n" );
		sb.append( "\t" + START_NODE + "[shape=point];\n\n" );
		for( int s = 0; s < states.length; s++ ) {
			sb.append( "\t"+s+"[" + states[s].getGraphvizNodeOptions( freq[s], maxFreq, nf ) + ",color=" + (finalState[s]?"red":"black")+ "];\n" );
		}
		
		if(rankPatterns != null){
			StringBuffer ranks = new StringBuffer();
			HashMap<String, IntList> map = new HashMap<String, IntList>();
			for(int s=0;s<name.length;s++){
				Iterator<String> it = rankPatterns.keySet().iterator();
				while( it.hasNext() ) {
					String key = it.next();
					if( name[s].matches( key ) ) {
						IntList list = map.get( key );
						if(list == null){
							map.put( key, new IntList() );
						}
						map.get( key ).add( s );
						break;
					}
				}
			}
			
			boolean startRanked = false;
			for(String key : map.keySet()){
				
				ranks.append( "{rank=" + rankPatterns.get(key) + "; " );
				if( !startRanked && START_NODE.matches( key ) ) {
					ranks.append( START_NODE + " " );
					startRanked = true;
				}
				IntList list = map.get( key );
				for(int i=0;i<list.length();i++){
					ranks.append( list.get( i )+" " );
				}
				
				ranks.append( ";}\n" );
			}
			sb.append( ranks );
		}
		
		sb.append( "\n" );
		sb.append( transition.getGraphizNetworkRepresentation( nf, null, data!=null ) );
		sb.append( "}" );
		return sb.toString();
	}
	
	private double[] getStateFreq( DataSet data, double[] weight ) throws Exception {
		double[] res = new double[states.length];
		if( data != null ) {			
			double w = 1, sum = 0;
			double[][] current = createMatrixForStatePosterior( 0, data.getMaximalElementLength()-1 );
			for( int i = 0; i < data.getNumberOfElements(); i++ ) {
				Sequence seq = data.getElementAt(i);
				fillLogStatePosteriorMatrix( current, 0, seq.getLength()-1, data.getElementAt(i), false );
				if( weight != null ) {
					w = weight[i];
				}
				sum += w;
				for( int s = 0; s < states.length; s++ ) {
					res[s] += w*Math.exp(Normalisation.getLogSum( current[s] ));
				}				
			}
			for( int s = 0; s < states.length; s++ ) {
				res[s] /= sum;
			}
		}
		return res;
	}
	
	/**
	 * This method creates an empty matrix for the log state posterior.
	 * 
	 * @param startPos the start position
	 * @param endPos the end position
	 * 
	 * @return an empty matrix for the log state posterior
	 * 
	 * @see #getLogStatePosteriorMatrixFor(int, int, Sequence)
	 * @see #fillLogStatePosteriorMatrix(double[][], int, int, Sequence, boolean)
	 */
	protected double[][] createMatrixForStatePosterior( int startPos, int endPos ) {
		return new double[states.length][endPos-startPos+1+1];
	}
	
	/**
	 * This method fills the log state posterior of Sequence <code>seq</code> in a given matrix.
	 * 
	 * @param statePosterior the matrix for the log state posterior
	 * @param startPos the start position
	 * @param endPos the end position
	 * @param seq the sequence
	 * @param silentZero <code>true</code> if the state posterior for silent states is defined to be zero, otherwise <code>false</code>
	 * 
	 * @throws Exception if an error occurs during the computation
	 * 
	 * @see #getLogStatePosteriorMatrixFor(int, int, Sequence)
	 * @see #createMatrixForStatePosterior(int, int)
	 */
	protected abstract void fillLogStatePosteriorMatrix( double[][] statePosterior, int startPos, int endPos, Sequence seq, boolean silentZero ) throws Exception;
		
	/**
	 * This method returns the log state posterior of all states for a sequence.
	 * 
	 * @param startPos the start position within the sequence
	 * @param endPos the end position within the sequence
	 * @param seq the sequence
	 * 
	 * @return the score for each state an each sequence position
	 * 
	 * @throws Exception if the state posterior could not be computed, for instance if the model is not trained, ...
	 */
	public double[][] getLogStatePosteriorMatrixFor( int startPos, int endPos, Sequence seq ) throws Exception {
		double[][] m = createMatrixForStatePosterior( startPos, endPos );
		fillLogStatePosteriorMatrix( m, startPos, endPos, seq, true );
		return getFinalStatePosterioriMatrix( m );
	}
	
	/**
	 * This method is used if {@link #fillLogStatePosteriorMatrix(double[][], int, int, Sequence, boolean)} is used with code>silentZero==true</code>
	 * to eliminate the first row.
	 * 
	 * @param intermediate the intermediate (log) state posterior matrix containing one additional row for silent states before the first emission
	 * 
	 * @return the final (log) state posterior matrix
	 */
	protected double[][] getFinalStatePosterioriMatrix( double[][] intermediate ) {
		double[][] res = new double[intermediate.length][];
		for( int i = 0; i < res.length; i++ ) {
			res[i] = new double[intermediate[i].length-1];
			System.arraycopy( intermediate[i], 1, res[i], 0, res[i].length );
		}
		return res;
	}
	
	/**
	 * This method returns the log state posterior of all states for a sequence.
	 * 
	 * @param seq the sequence
	 * 
	 * @return the score for each state an each sequence position
	 * 
	 * @throws Exception if the state posterior could not be computed, for instance if the model is not trained, ...
	 * 
	 * @see #getLogStatePosteriorMatrixFor(int, int, Sequence)
	 */
	public double[][] getStatePosteriorMatrixFor( Sequence seq ) throws Exception {
		double[][] matrix = getLogStatePosteriorMatrixFor( 0, seq.getLength()-1, seq );
		for( int i = 0; i < matrix.length; i++ ) {
			for( int j = 0; j < matrix[i].length; j++ ) {
				matrix[i][j] = Math.exp( matrix[i][j] );
			}
		}
		return matrix;
	}
	
	/**
	 * This method returns the log state posteriors for all sequences of the sample <code>data</code>.
	 * 
	 * @param data the sequences
	 * 
	 * @return the log state posterior matrices for all sequences
	 * 
	 * @throws Exception if the state posterior could not be computed, for instance if the model is not trained, ...
	 * 
	 * @see #getLogStatePosteriorMatrixFor(int, int, Sequence)
	 */
	public double[][][] getLogStatePosteriorMatrixFor( DataSet data ) throws Exception {
		double[][][] matrix = new double[data.getNumberOfElements()][][];
		for( int i = 0; i < matrix.length; i++ ) {
			Sequence s = data.getElementAt(i);
			matrix[i] = getLogStatePosteriorMatrixFor( 0, s.getLength()-1, s );
		}
		return matrix;
	}
	
	/**
	 * This method returns the state posteriors for all sequences of the sample <code>data</code>.
	 * 
	 * @param data the sequences
	 * 
	 * @return the state posterior matrices for all sequences
	 * 
	 * @throws Exception if the state posterior could not be computed, for instance if the model is not trained, ...
	 * 
	 * @see #getStatePosteriorMatrixFor(Sequence)
	 */
	public double[][][] getStatePosteriorMatrixFor( DataSet data ) throws Exception {
		double[][][] matrix = new double[data.getNumberOfElements()][][];
		for( int i = 0; i < matrix.length; i++ ) {
			matrix[i] = getStatePosteriorMatrixFor( data.getElementAt(i) );
		}
		return matrix;
	}
	
	/**
	 * @param startPos the start position within the sequence
	 * @param endPos the end position within the sequence
	 * @param seq the sequence
	 * 
	 * @return a {@link Pair} containing the viterbi state path and the corresponding score
	 * 
	 * @throws Exception if the viterbi path could not be computed, for instance if the model is not trained, ...
	 */
	public abstract Pair<IntList,Double> getViterbiPathFor(int startPos, int endPos, Sequence seq ) throws Exception;

	/**
	 * @param seq the sequence
	 * 
	 * @return a {@link Pair} containing the viterbi state path and the corresponding score
	 * 
	 * @throws Exception if the viterbi path could not be computed, for instance if the model is not trained, ...
	 * 
	 * @see #getViterbiPathFor(int, int, Sequence)
	 */
	public Pair<IntList,Double> getViterbiPathFor( Sequence seq ) throws Exception {
		return getViterbiPathFor( 0, seq.getLength()-1, seq );
	}
	
	/**
	 * This method returns the viterbi paths and scores for all sequences of the sample <code>data</code>.
	 * 
	 * @param data the sequences
	 * 
	 * @return the viterbi paths and scores for all sequences
	 * 
	 * @throws Exception if the viterbi paths and scores could not be computed, for instance if the model is not trained, ...
	 * 
	 * @see #getViterbiPathFor(Sequence)
	 */
	public Pair<IntList,Double>[] getViterbiPathsFor( DataSet data ) throws Exception {
		Pair<IntList,Double>[] matrix = new Pair[data.getNumberOfElements()];
		for( int i = 0; i < matrix.length; i++ ) {
			matrix[i] = getViterbiPathFor( data.getElementAt(i) );
		}
		return matrix;
	}
	
	/**
	 * This method decodes any path of the HMM, i.e. it converts the integer representation of the path in a String representation.
	 * @param path the path in integer representation
	 * @return the path in String representation, the i-th entry is the i-th visited state
	 * 
	 * @see #getViterbiPathFor(Sequence)
	 * @see #getViterbiPathFor(int, int, Sequence)
	 */
	public final String[] decodePath( IntList path ) {
		String[] decoded = new String[path.length()];
		for( int i = 0; i < decoded.length; i++ ) {
			decoded[i] = name[path.get(i)];
		}
		return decoded;
	}
	
	/**
	 * @param path the given state path
	 * @param startPos the start position within the sequence(s) (inclusive)
	 * @param seq the sequence(s)
	 * 
	 * @return the logarithm of the probability for the given path and the given sequence(s)
	 * 
	 * @throws Exception if the probability for the sequence given path could not be computed, for instance if the model is not trained, ... 
	 */
	public abstract double getLogProbForPath(IntList path, int startPos, Sequence seq ) throws Exception;
	
	/**
	 * This method instantiates all helper variables that are need inside the model for instance for filling forward and backward matrix, ...
	 */
	protected abstract void createHelperVariables();
	
	/**
	* This method invokes the method {@link #createHelperVariables()} and provides the matrix with given type. Type 0 stands for {@link AbstractHMM#fwdMatrix}, and type 1 stands for {@link AbstractHMM#bwdMatrix}.
	*
	* @param type the type of the matrix
	* @param length the maximal sequence length
	*/
	protected void provideMatrix( int type, int length ) {
		createHelperVariables();
		length++;//because of silent states
		double[][] matrix;
		switch( type ) {
			case 0: matrix = fwdMatrix; break;
			case 1: matrix = bwdMatrix; break;
			default:
				 throw new IllegalArgumentException( "unknown matrix type" );
		}
		if( matrix == null || matrix.length < length ) {
			matrix = new double[length][];
			int maxOrder = transition.getMaximalMarkovOrder(), dim = -1, l = 0;
			for( l = 0; l <= maxOrder && l < length; l++ ) {
				dim = transition.getNumberOfIndexes( l );
				matrix[l] = new double[dim];
			}
			while( l < length ) {
				matrix[l++] = new double[dim];
			}
		}
		for( int l = 0; l < length; l++ ) {
			Arrays.fill( matrix[l], Double.NEGATIVE_INFINITY );
		}
		switch( type ) {
			case 0: fwdMatrix = matrix; break;
			case 1: bwdMatrix = matrix; break;
		}
	}
	
	/**
	 * This method returns the number of the (hidden) states
	 * 
	 * @return the number of the states
	 */
	public int getNumberOfStates() {
		return states.length;
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getLogProbFor(de.jstacs.data.Sequence, int,
	 * int)
	 */
	public double getLogProbFor(Sequence sequence, int startpos, int endpos)
			throws Exception {
		int l = endpos - startpos+1, len = getLength();
		if( !sequence.getAlphabetContainer().checkConsistency(getAlphabetContainer()) ) {
			throw new WrongAlphabetException( "The AlphabetContainer of the sequence and the model do not match." );
		}
		if( len != 0 && l != len ) {
			throw new WrongLengthException( "The given start position ("+ startpos + ") and end position (" + endpos + ") yield an length of " + l + " which is not possible for the current model that models sequences of length " + len + "." );
		}
		return logProb( startpos, endpos, sequence );
	}
	
	/**
	 * This method creates an {@link RuntimeException} from any other {@link Exception}
	 * 
	 * @param e the {@link Exception}
	 * 
	 * @return a {@link RuntimeException}
	 */
	protected static RuntimeException getRunTimeException( Exception e ) {
		RuntimeException re;
		if( e instanceof RuntimeException ) {
			re = (RuntimeException) e;
		} else {
			re = new RuntimeException( e.getMessage() );
			re.setStackTrace( e.getStackTrace() );
		}
		return re;
	}
	
	/**
	 * This method computes the logarithm of the probability of the corresponding subsequences.
	 * The method does not check the {@link AlphabetContainer} and possible further features
	 * before starting the computation.
	 * 
	 * @param startpos the start position (inclusive)
	 * @param endpos the end position (inclusive)
	 * @param sequence the {@link Sequence}(s)
	 * 
	 * @return the logarithm of the probability
	 * 
	 * @throws Exception if the model has no parameters (for instance if it is not trained)
	 */
	protected double logProb( int startpos, int endpos, Sequence sequence ) throws Exception {
		try {
			fillBwdMatrix(startpos, endpos, sequence);
		} catch( Exception e ) {
			throw getRunTimeException( e );
		}
		return bwdMatrix[0][0];
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#train(de.jstacs.data.DataSet)
	 */
	public void train(DataSet data) throws Exception {
		train( data, null );
	}
	
	/**
	 * Sets the {@link OutputStream} that is used e.g. for writing information
	 * while training. It is possible to set <code>o=null</code>, than nothing
	 * will be written.
	 * 
	 * @param o the {@link OutputStream}
	 */
	public final void setOutputStream( OutputStream o ) {
		sostream = SafeOutputStream.getSafeOutputStream( o );
	}
	
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#finalize()
	 */
	protected void finalize() throws Throwable {
		transition = null;
		states = null;
		trainingParameter = null;
		bwdMatrix = null;
		fwdMatrix = null;
		super.finalize();
	}
	
	/**
	 * This method determines the final states of the HMM.
	 * 
	 * @see #finalState
	 */
	protected void determineFinalStates() {
		finalState = transition.isAbsoring();
		int i = 0;
		while( i < finalState.length && !finalState[i] ){
			i++;
		}
		
		if( i == finalState.length ) {
			for( i = 0; i < finalState.length; i++ ) {
				finalState[i] = !states[i].isSilent(); 
			}
		}
	}
	
	/**
	 * The method returns the decoded state posterior, i.e. a sequence of states.
	 * Be careful: For HMMs that do not allow all transitions between states this sequence does not have to be a valid path of the HMM.
	 * 
	 * @param statePosterior the (log) state posterior(s)
	 * 
	 * @return the decoded state posterior(s)
	 * 
	 * @see #getLogStatePosteriorMatrixFor(int, int, Sequence)
	 */
	public static int[][] decodeStatePosterior( double[][]... statePosterior ) {
		int[][] res = new int[statePosterior.length][];
		for( int s = 0; s < res.length; s++ ) {
			res[s] = new int[statePosterior[s][0].length];
			for( int l = 0; l < res[s].length; l++ ) {
				res[s][l] = 0;
				for( int j = 1; j < statePosterior[s].length; j++ ) {
					if( statePosterior[s][res[s][l]][l] < statePosterior[s][j][l] ) {
						res[s][l] = j;
					}
				}
			}
		}
		return res;
	}
	
	public String toString() {
		String res = "Transition:\n-----------\n" + transition.toString( name ); 
		res += "\nStates:\n-------\n";
		for( int e = 0; e < states.length; e++ ) {
			res += states[e] + "\n";
		}		
		return res;
	}
}