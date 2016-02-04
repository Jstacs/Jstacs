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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements;

import java.text.NumberFormat;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * Scaled transition element for an HMM with scaled transition matrices (SHMM).
 * 
 * @author Michael Seifert, Jens Keilwagen
 */
public class ScaledTransitionElement extends ReferenceBasedTransitionElement {
	
	/**
	 * The Graphviz options for distinguishing the arrows of the different transitions. 
	 */
	protected String[] arrowOptions;	
	
	/**
	 * The scaling factors of the individual transition classes.
	 */
	protected double[] scalingFactor;
	
	/**
	 * Creates an object representing the transition probabilities of a Hidden Markov TrainableStatisticalModel with scaled transition matrices (SHMM) for the given context. 
	 *  
	 * @param context context of states
	 * 
	 * @param states states that can be reached from the context
	 * 
	 * @param probabilities transition probabilities from the context to the states
	 * 
	 * @param ess ess
	 * 
	 * @param scalingFactor scaling factor
	 * 
	 * @param annotationID identifier for decoding transition classes 
	 */
	public ScaledTransitionElement( int[] context, int[] states, double[] probabilities, double ess, double[] scalingFactor, String annotationID ) {
		super(context, states, ess, probabilities, annotationID );
		this.scalingFactor = scalingFactor.clone();
		init();
	}
	
	protected void init() {
		super.init();
		if( scalingFactor != null ) {
			int numberOfTransitionClasses = this.scalingFactor.length;
			this.statisticsTransitionProb = new double[ numberOfTransitionClasses ][ states.length ];
			this.arrowOptions = new String[ numberOfTransitionClasses ];
			for( int i = 0; i < numberOfTransitionClasses; i++ )
			{
				this.arrowOptions[ i ] = "style=\"setlinewidth(" + (i + 1 )+ ")\"";
			}
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ScaledTransitionElement} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ScaledTransitionElement} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public ScaledTransitionElement(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	public ScaledTransitionElement clone() throws CloneNotSupportedException {
		ScaledTransitionElement clone = (ScaledTransitionElement) super.clone();
		clone.arrowOptions = arrowOptions.clone();
		clone.scalingFactor = scalingFactor.clone();
		return clone;
	}

	/**
	 * Returns the index of the transition matrix used for the transition from <code>pos - 1</code>  to <code>pos</code> in sequences <code>seq</code>.
	 * <ul>
	 * 	<li>corresponding transition class is assigned to <code>pos</code> 
	 * </ul>
	 * 
	 * @param pos position in seq
	 * 
	 * @param seq sequences
	 * 
	 * @return index of transition matrix
	 */
	protected int getIndex( int pos, Sequence seq )	{		
		return ((ReferenceSequenceAnnotation)seq.getSequenceAnnotationByTypeAndIdentifier( ReferenceSequenceAnnotation.TYPE, annotationID )).getReferenceSequence().discreteVal( pos );
	}
	
	public void addToStatistic( int childIdx, double weight, Sequence sequence, int sequencePosition ) {
		//context is already known becaues of transition element
		//statistics is summed up for the next state
		//transition class for the transition from context to the next state is stored at sequence position
		//weight is the epsilon in my world
		//[transition class][next state]
		this.statisticsTransitionProb[ this.getIndex( sequencePosition, sequence ) ][ childIdx ] += weight;	
	}
	
	public void resetStatistic() {
		for( int c = 0; c < statisticsTransitionProb.length; c++ )
		{
			for( int j = 0; j < this.states.length; j++ )
			{
				this.statisticsTransitionProb[ c ][ j ] = 0; 
			}
		}
		
		//Resets each element of the 'statistic' array to 'ess'  
		super.resetStatistic();
	}
	
	public void estimateFromStatistic()
	{
		//Estimate transition probabilities
		double sumOfGammas = 0;
		for( int i = 0; i < states.length; i++ ) {
			for( int c = 0; c < statisticsTransitionProb.length; c++ ) {
				statistic[ i ] += this.statisticsTransitionProb[ c ][ i ];
			}
			sumOfGammas += statistic[i];
		}
		logNorm = 0;
		
		//Self-transition probabilities
		parameters[diagElement] = this.determineDiagonalElement( sumOfGammas, 1E-3, true, false );
			
		//Non-self transition probabilities
		for( int j = 0; j < states.length; j++ )
		{
			if( j != diagElement )
			{
				//Lambda
				//this.probabilities[ 1 ][ i ][ j ] = ( 1 - this.probabilities[ 1 ][ i ][ i ] ) * sumOfEpsilons[ i ][ j ] / ( sumOfGammas[ i ] - sumOfEpsilons[ i ][ i ] );
				parameters[j] = ( 1 - this.parameters[ diagElement ] ) * statistic[ j ] / ( sumOfGammas - statistic[ diagElement ] );
			}
		}
	}
	
	/**
	 * Determines the diagonal element (self-transition probability) for each state in MAP training.
	 * 
	 * @param sumOfGammas summed gammas
	 * 
	 * @param state state (0 = '-'; 1 = '='; 2 = '+')
	 * 
	 * @param firstStep step in newton's method
	 * 
	 * @param output output
	 */
	private double determineDiagonalElement( double sumOfGammas, double currentDiagonalElement, boolean firstStep, boolean output )
	{
		String stateSymbol = "" + diagElement;
		
		if( firstStep && output ) {
			System.out.println( "<------ Start '" + stateSymbol + "' ------>" );
			System.out.println( "\tdiag = " + currentDiagonalElement );
		}		
		
		double nextDiagonalElement;
		
		
		//compute f_i^(1)(currentDiagonalElement) := L(currentDiagonalElement) and f_i^(1)'
		double variablePart = 0;
		double variablePartDerivative = 0;

		for( int c = 0; c < scalingFactor.length; c++ )
		{
			variablePart += this.scalingFactor[ c ] / ( currentDiagonalElement - 1 + this.scalingFactor[ c ] ) * statisticsTransitionProb[ c ][ diagElement ];
			variablePartDerivative -= this.scalingFactor[ c ] / Math.pow( currentDiagonalElement - 1 + this.scalingFactor[ c ], 2.0 ) * statisticsTransitionProb[ c ][ diagElement ];
		}
		
		//add transition prior
		//Lambda
		variablePart += this.hyperParameters[ diagElement ] / currentDiagonalElement;
		variablePartDerivative -= this.hyperParameters[ diagElement ] / Math.pow( currentDiagonalElement, 2.0 );

		
		//compute nextDiagonalElement
		//Like in Diss
		nextDiagonalElement = currentDiagonalElement - ( variablePart - sumOfGammas ) / variablePartDerivative;


		//Test for root: Ineffektiv kann verbessert werden
		//compute f_i^(1)(nextDiagonalElement)
		variablePart = 0;
		for( int c = 0; c < scalingFactor.length; c++ )
		{
			variablePart += this.scalingFactor[ c ] / ( nextDiagonalElement - 1 + this.scalingFactor[ c ] ) * statisticsTransitionProb[ c ][ diagElement ];
		}
		
		//add transition prior
		//Lambda
		variablePart += this.hyperParameters[ diagElement ] / nextDiagonalElement;	
		
		if( output )
		{
			System.out.println( "\tdiag = " + nextDiagonalElement + "\t\tf = " + ( variablePart - sumOfGammas ) );
		}
	
		
		//one more iteration?
		if( Math.abs( nextDiagonalElement - currentDiagonalElement ) > 1E-10 )
		{
			nextDiagonalElement = this.determineDiagonalElement( sumOfGammas, nextDiagonalElement, false, output );
		}
		
		if( firstStep )
		{			
			if( output )
			{
				System.out.println( "<------ End   '" + stateSymbol + "' ------>" );
			}
		}
		return nextDiagonalElement;
	}
	

	public double getLogScoreFor( int state, Sequence sequence, int sequencePosition ) {
		return( Math.log( getTransitionProb( state, this.getIndex( sequencePosition, sequence ) ) ) );
	}
	
	/**
	 * Returns the probability for a transition to state 'j' in transition class 'tc'.
	 * 
	 * @param j next state
	 * @param tc transition class
	 * @return transition probability
	 */
	private double getTransitionProb( int j, int tc ) {
		double transitionProb;
		if( j == diagElement ) {
			//Self-Transition
			transitionProb = ( this.parameters[ j ] - 1 + this.scalingFactor[ tc ] ) / this.scalingFactor[ tc ];
		} else {
			//Non-Self-Transition
			transitionProb = this.parameters[ j ] / this.scalingFactor[ tc ];
		}
			
		return( transitionProb );
	}
	
	protected void appendTransitions( StringBuffer representation, String contextNodeRepresentation, NumberFormat nf, String arrowOption, boolean graphical  ) {
		for( int s = 0; s < states.length; s++ ) {
			for( int tc = 0; tc < scalingFactor.length; tc++ )
			{
				representation.append( "\t" + contextNodeRepresentation + "->" + states[s] + getArrowOption( nf, this.getTransitionProb( s, tc ), getGraphvizEdgeWeight(s), this.arrowOptions[ tc ], graphical ) + "\n" );	
			}
		}
	}
	
	protected void appendFurtherInformation( StringBuffer xml ) {
		super.appendFurtherInformation( xml );
		XMLParser.appendObjectWithTags( xml, scalingFactor, "scalingFactor" );
	}

	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		super.extractFurtherInformation( xml );
		scalingFactor = (double[]) XMLParser.extractObjectForTags( xml, "scalingFactor" );
	}
	
	protected String getXMLTag() {
		return XML_TAG;
	}
	
	private static final String XML_TAG = "SCALED_TRANSITION_ELEMENT";
	
	/**
	 * Returns a string representation of the transition probabilities.
	 * @param stateNames the names of the states
	 * @return representation of transition probabilities
	 */
	public String toString( String[] stateNames ) {
		if( parameters.length > 0 ) {
			StringBuffer sb = new StringBuffer();
			String context = getContext( stateNames ), c;
			for( int tc = 0; tc < scalingFactor.length; tc++ ) {
				c = context;
				if( c.length() == 0 ) {
					c = "|";
				} else {
					c +=",";
				}
				c += "tc=" + tc;
				for( int i = 0; i < parameters.length; i++ ) {
					sb.append("P(" + getLabel( stateNames, states[i] ) + c + ") \t= " + getTransitionProb(i, tc) );
					sb.append("\t");
				}
				sb.append( "\n" );
			}
			return sb.toString();
		} else {
			return "";
		}
	}
}