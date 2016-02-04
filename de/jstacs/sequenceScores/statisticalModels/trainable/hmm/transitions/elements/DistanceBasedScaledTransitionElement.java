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
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * Distance-based scaled transition element for an HMM with distance-scaled transition matrices (DSHMM).
 * 
 * @author Michael Seifert
 *
 */
public class DistanceBasedScaledTransitionElement extends ReferenceBasedTransitionElement
{
	/**
	 * The Graphviz options for drawing the arrows. 
	 */
	protected String arrowOptions;	
	
	/**
	 * The maximal scaling factor.
	 */
	protected double scalingFactor;
	
	/**
	 * Represents the summarized epsilons required for estimating the transition probabilities from the <code>context</code>.
	 * <ul>
	 * 	<li> statisticsTransitionProb[j]: epsilon( context, states[j] )
	 *  <li> No pseudocounts included.
	 * </ul>
	 */
	protected double[] statisticsTransitionProb;//TODO: masks field in superclass. intended?
	
	/**
	 * Contains the single epsilons of the diagonal elements required for estimating the self-transition probability.
	 * <ul>
	 *  <li> reset should not be necessary after each training step
	 * </ul>
	 */
	protected Hashtable<Sequence,double []> diagonalWeights; 
	
	/**
	 * Creates an object representing the transition probabilities of a Hidden Markov TrainableStatisticalModel with scaled transition matrices (SHMM) for the given context. 
	 *  
	 * @param context context of states
	 * 
	 * @param states states that can be reached from the context
	 * 
	 * @param probabilities transition probabilities from the context to the states
	 * 
	 * @param ess the equivalent sample size (ess)
	 * 
	 * @param scalingFactor scaling factor
	 * 
	 * @param annotationID identifier for decoding transition classes 
	 * 
	 * @see #DistanceBasedScaledTransitionElement(int[], int[], double[], double, double, String, double[])
	 */
	public DistanceBasedScaledTransitionElement( int[] context, int[] states, double[] probabilities, double ess, double scalingFactor, String annotationID ) {
		this( context, states, probabilities, ess, scalingFactor, annotationID, null );
	}
	
	/**
	 * Creates an object representing the transition probabilities of a Hidden Markov TrainableStatisticalModel with scaled transition matrices (SHMM) for the given context. 
	 *  
	 * @param context context of states
	 * 
	 * @param states states that can be reached from the context
	 * 
	 * @param probabilities transition probabilities from the context to the states
	 * 
	 * @param ess the equivalent sample size (ess)
	 * 
	 * @param scalingFactor scaling factor
	 * 
	 * @param annotationID identifier for decoding transition classes 
	 * 
	 * @param weight the weight for plotting the edge in Graphviz, enables to modify the edge length, larger weights imply shorter edges (default: 1)
	 */
	public DistanceBasedScaledTransitionElement( int[] context, int[] states, double[] probabilities, double ess, double scalingFactor, String annotationID, double[] weight ) {
		super( context, states, ess, probabilities, annotationID, weight );
		this.scalingFactor = scalingFactor;
		init();
	}
	
	/**
	 * Basic initialization.
	 */
	protected void init()
	{
		super.init();
		this.statisticsTransitionProb = new double[ states.length ];
		this.diagonalWeights = new Hashtable<Sequence, double[]>();

		this.arrowOptions = "style=\"setlinewidth(1)\"";
	}

	/**
	 * Extracts a distance-base scaled transition element from XML.
	 * @param xml the XML-representation
	 * @throws NonParsableException if the representation could not be parsed
	 */
	public DistanceBasedScaledTransitionElement(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	/**
	 * Clones the object.
	 */
	public DistanceBasedScaledTransitionElement clone() throws CloneNotSupportedException
	{
		DistanceBasedScaledTransitionElement clone = (DistanceBasedScaledTransitionElement) super.clone();
		clone.scalingFactor = scalingFactor;
		return clone;
	}
	
	/**
	 * Returns the distance integrated into the transition from <code>pos - 1</code>  to <code>pos</code> in sequences <code>seq</code>.
	 * 
	 * @param pos position in seq
	 * 
	 * @param seq sequences
	 * 
	 * @return distance
	 */
	protected double getIndex( int pos, Sequence seq )
	{
		return ((ReferenceSequenceAnnotation)seq.getSequenceAnnotationByTypeAndIdentifier( ReferenceSequenceAnnotation.TYPE, annotationID )).getReferenceSequence().discreteVal( pos );
	}

	public void addToStatistic( int childIdx, double weight, Sequence sequence, int sequencePosition )
	{				
		//Summation: sum of epsilons
		//Was genau ist childIdx -> soll man den auf children mappen?
		this.statisticsTransitionProb[ childIdx ] += weight;
		
		//Minimaler Wert von sequencePosition ist 1???
		if( sequencePosition == 0 )
		{
			System.out.println( "Hier seqPos = 0" );
		}

		
		//Additional separate storage of individual epsilons for self-transitions 
		if( childIdx == this.diagElement )
		{
			boolean contained = this.diagonalWeights.containsKey( sequence );
			if( !contained )
			{				
				//this.diagonalWeights.put( sequence, new double[ sequence.getLength() - 1 ] );
				this.diagonalWeights.put( sequence, new double[ sequence.getLength() ] );
			}
			
			double[] epsilons = this.diagonalWeights.get( sequence );
			epsilons[ sequencePosition ] = weight;
		}
	}
	
	public void resetStatistic()
	{
		//Reset count statistics
		for( int j = 0; j < this.states.length; j++ )
		{
			this.statisticsTransitionProb[ j ] = 0; 
		}
	

		if( !this.diagonalWeights.isEmpty() )
		{
			Enumeration<Sequence> keyIterator = this.diagonalWeights.keys();	
			
			do
			{
				Sequence seq = keyIterator.nextElement();
				double[] diagEpsilons = this.diagonalWeights.get( seq );			
				Arrays.fill( diagEpsilons, 0 );		
				
			} while( keyIterator.hasMoreElements() );
		}

		
		//Resets each element of the 'statistic' array to 'ess'  
		super.resetStatistic();
	}
	
	public void estimateFromStatistic()
	{
		//Estimate transition probabilities
			//- 'statistic' already contains the value of the hyper-parameter
			//- 'sumOfGammas' contains the value of the hyper-parameter after summing
		double sumOfGammas = 0;
		for( int i = 0; i < states.length; i++ )
		{
			statistic[ i ] += this.statisticsTransitionProb[ i ];

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
		
		
		Enumeration<Sequence> keyIterator = this.diagonalWeights.keys();
				
		
		//compute f_i^(1)(currentDiagonalElement) := L(currentDiagonalElement) and f_i^(1)'
		double variablePart = 0;
		double variablePartDerivative = 0;

		
		do
		{
			Sequence seq = keyIterator.nextElement();
			int T = seq.getLength();
			double[] diagEpsilons = this.diagonalWeights.get( seq );
						
			
			/*
			for( int t = 0; t < T - 1; t++ )
			{
				//double tc = this.getIndex( t, seq );
				double tc = this.getIndex( t + 1, seq );
				double distanceBasedScalingFactor = 1 + ( this.scalingFactor - 1 ) * ( 1 - tc );				
				
				//variablePart += distanceBasedScalingFactor / ( currentDiagonalElement - 1 + distanceBasedScalingFactor ) * diagEpsilons[ t ];
				//variablePartDerivative -= distanceBasedScalingFactor / Math.pow( currentDiagonalElement - 1 + distanceBasedScalingFactor, 2.0 ) * diagEpsilons[ t ];
				
				//corrected because weights in diagEpsilons start at position 1
				variablePart += distanceBasedScalingFactor / ( currentDiagonalElement - 1 + distanceBasedScalingFactor ) * diagEpsilons[ t + 1 ];
				variablePartDerivative -= distanceBasedScalingFactor / Math.pow( currentDiagonalElement - 1 + distanceBasedScalingFactor, 2.0 ) * diagEpsilons[ t + 1 ];
			}
			*/
			
			
			for( int t = this.context.length; t < T; t++ )
			{
				double tc = this.getIndex( t, seq );				
				double distanceBasedScalingFactor = 1 + ( this.scalingFactor - 1 ) * ( 1 - tc );				
				
				//variablePart += distanceBasedScalingFactor / ( currentDiagonalElement - 1 + distanceBasedScalingFactor ) * diagEpsilons[ t ];
				//variablePartDerivative -= distanceBasedScalingFactor / Math.pow( currentDiagonalElement - 1 + distanceBasedScalingFactor, 2.0 ) * diagEpsilons[ t ];
				
				//corrected because weights in diagEpsilons start at position 1
				//variablePart += distanceBasedScalingFactor / ( currentDiagonalElement - 1 + distanceBasedScalingFactor ) * diagEpsilons[ t + 1 ];
				variablePart += distanceBasedScalingFactor / ( currentDiagonalElement - 1 + distanceBasedScalingFactor ) * diagEpsilons[ t ];
				//variablePartDerivative -= distanceBasedScalingFactor / Math.pow( currentDiagonalElement - 1 + distanceBasedScalingFactor, 2.0 ) * diagEpsilons[ t + 1 ];
				variablePartDerivative -= distanceBasedScalingFactor / Math.pow( currentDiagonalElement - 1 + distanceBasedScalingFactor, 2.0 ) * diagEpsilons[ t ];
			}
			
		} while( keyIterator.hasMoreElements() );
		
		
		//add transition prior
		//Lambda
		variablePart += this.hyperParameters[ diagElement ] / currentDiagonalElement;
		variablePartDerivative -= this.hyperParameters[ diagElement ] / Math.pow( currentDiagonalElement, 2.0 );

		
		//compute nextDiagonalElement
		//Like in Diss
		nextDiagonalElement = currentDiagonalElement - ( variablePart - sumOfGammas ) / variablePartDerivative;


		//Test for root
		//compute f_i^(1)(nextDiagonalElement)
		keyIterator = this.diagonalWeights.keys();
		
		variablePart = 0;
		do
		{
			Sequence seq = keyIterator.nextElement();
			int T = seq.getLength();
			double[] diagEpsilons = this.diagonalWeights.get( seq );

			for( int t = this.context.length; t < T; t++ )
			{
				double tc = this.getIndex( t, seq );
				//double tc = this.getIndex( t + 1, seq );
				double distanceBasedScalingFactor = 1 + ( this.scalingFactor - 1 ) * ( 1 - tc );				
				
				//variablePart += distanceBasedScalingFactor / ( nextDiagonalElement - 1 + distanceBasedScalingFactor ) * diagEpsilons[ t ];
				//corrected because weights in diagEpsilons start at position 1
				//variablePart += distanceBasedScalingFactor / ( nextDiagonalElement - 1 + distanceBasedScalingFactor ) * diagEpsilons[ t + 1 ];
				variablePart += distanceBasedScalingFactor / ( nextDiagonalElement - 1 + distanceBasedScalingFactor ) * diagEpsilons[ t ];

			}

			
			/*
			for( int t = 0; t < T - 1; t++ )
			{
				//double tc = this.getIndex( t, seq );
				double tc = this.getIndex( t + 1, seq );
				double distanceBasedScalingFactor = 1 + ( this.scalingFactor - 1 ) * ( 1 - tc );				
				
				//variablePart += distanceBasedScalingFactor / ( nextDiagonalElement - 1 + distanceBasedScalingFactor ) * diagEpsilons[ t ];
				//corrected because weights in diagEpsilons start at position 1
				variablePart += distanceBasedScalingFactor / ( nextDiagonalElement - 1 + distanceBasedScalingFactor ) * diagEpsilons[ t + 1 ];

			}
			*/			
			
		} while( keyIterator.hasMoreElements() );
		
		
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
	
	public double getLogScoreFor( int state, Sequence sequence, int sequencePosition )
	{
		return( Math.log( getTransitionProb( state, this.getIndex( sequencePosition, sequence ) ) ) );
	}
	
	/**
	 * Returns the probability for a transition to state 'j' in transition class 'tc'.
	 * 
	 * @param j next state
	 * @param tc transition class
	 * @return transition probability
	 */
	private double getTransitionProb( int j, double tc )
	{
		double transitionProb;
		double distanceBasedScalingFactor = 1 + ( this.scalingFactor - 1 ) * ( 1 - tc );
		if( j == diagElement )
		{
			//Self-Transition
			transitionProb = ( this.parameters[ j ] - 1 + distanceBasedScalingFactor ) / distanceBasedScalingFactor;
		}
		else
		{
			//Non-Self-Transition
			transitionProb = this.parameters[ j ] / distanceBasedScalingFactor;
		}
			
		return( transitionProb );
	}
	
	protected void appendTransitions( StringBuffer res, String contextNodeRepresentation, NumberFormat nf, String arrowOption, boolean graphical  )
	{
		//Just show the basic transition probabilities
		for( int s = 0; s < states.length; s++ )
		{
				res.append( "\t" + contextNodeRepresentation + "->" + states[s] + getArrowOption( nf, parameters[ s ], getGraphvizEdgeWeight(s), "\n", graphical ) );	
		}
	}
	
	protected void appendFurtherInformation( StringBuffer xml )
	{
		super.appendFurtherInformation( xml );
		XMLParser.appendObjectWithTags( xml, scalingFactor, "scalingFactor" );
	}

	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException
	{
		super.extractFurtherInformation( xml );
		scalingFactor = (Double) XMLParser.extractObjectForTags( xml, "scalingFactor" );
	}
	
	protected String getXMLTag() {
		return XML_TAG;
	}
	
	private static final String XML_TAG = "DISTANCE_BASED_SCALED_TRANSITION_ELEMENT";
	
	/**
	 * Returns a string representation of the transition probabilities.
	 * @param stateNames the names of the states
	 * @return representation of transition probabilities
	 */
	public String toString( String[] stateNames )
	{
		if( parameters.length > 0 )
		{
			StringBuffer sb = new StringBuffer();
			String context = getContext( stateNames ), c;
		
			c = context;
			if( c.length() == 0 )
			{
				c = "|";
			}
	
			for( int i = 0; i < parameters.length; i++ )
			{
				sb.append("P(" + getLabel( stateNames, states[i] ) + c + ") \t= " + getTransitionProb( i, 1 ) );
				sb.append("\t");
			}
			sb.append( "\n" );
		
			return sb.toString();
		}
		else
		{
			return "";
		}
	}
}