package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements;

import de.jstacs.NonParsableException;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;

/**
 * This class implements transition elements that utilize a reference sequence to determine the transition probability.
 * The reference sequence is given as annotation of the original sequence and can be utilized using {@link #annotationID}.
 * 
 * @author Jens Keilwagen, Michael Seifert
 * 
 * @see de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation
 */
public class ReferenceBasedTransitionElement extends BasicTransitionElement {

	/**
	 * The annotation id used for determining the transition matrix from the {@link de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation}s.
	 */
	protected String annotationID;
	
	/**
	 * Represents the initial the transition probabilities.   
	 */
	protected double[] probabilities;
	
	/**
	 * Represents the gammas required for estimating the transition probabilities not including pseudocounts.
	 */
	protected double[][] statisticsTransitionProb;
	
	/**
	 * The equivalent sample size (ess) used in the prior of this instance.
	 */
	protected double ess;
	
	/**
	 * The index of the self transition.
	 */
	protected int diagElement;
	
	private static double[] getHyperParameters( double ess, int states ) {
		double[] hyper = new double[states];
		ess /= states;
		for( int i = 0; i < states; i++ ) {
			hyper[i] = ess;
		}
		return hyper;
	}
	
	/**
	 * Creates an {@link ReferenceBasedTransitionElement} representing a conditional transition probability distribution of a Hidden Markov TrainableStatisticalModel for the given context. 
	 *  
	 * @param context context of states
	 * 
	 * @param states states that can be reached from the context
	 * 
	 * @param probabilities transition probabilities from the context to the states
	 * 
	 * @param ess the equivalent sample size (ess)
	 * 
	 * @param annotationID identifier for decoding transition classes 
	 * 
	 * @see #ReferenceBasedTransitionElement(int[], int[], double, double[], String, double[])
	 */
	public ReferenceBasedTransitionElement( int[] context, int[] states, double ess, double[] probabilities, String annotationID ) {
		this( context, states, ess, probabilities, annotationID, null );
	}
	
	/**
	 * Creates an {@link ReferenceBasedTransitionElement} representing a conditional transition probability distribution of a Hidden Markov TrainableStatisticalModel for the given context. 
	 *  
	 * @param context context of states
	 * 
	 * @param states states that can be reached from the context
	 * 
	 * @param probabilities transition probabilities from the context to the states
	 * 
	 * @param ess the equivalent sample size (ess)
	 * 
	 * @param annotationID identifier for decoding transition classes
	 * 
	 * @param weight the weight for plotting the edge in Graphviz, enables to modify the edge length, larger weights imply shorter edges (default: 1)
	 */
	public ReferenceBasedTransitionElement( int[] context, int[] states, double ess, double[] probabilities, String annotationID, double[] weight ) {
		super(context, states, getHyperParameters( ess, states.length ), weight );
		this.probabilities = probabilities.clone();
		this.initializeRandomly();
		this.annotationID  = annotationID;
		this.ess = ess;
		init();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ReferenceBasedTransitionElement} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ReferenceBasedTransitionElement} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public ReferenceBasedTransitionElement(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	protected void init() {
		determineDiagonalElement();
		super.init();
	}
	
	public ReferenceBasedTransitionElement clone() throws CloneNotSupportedException {
		ReferenceBasedTransitionElement clone = (ReferenceBasedTransitionElement) super.clone();
		clone.probabilities = probabilities.clone();
		clone.statisticsTransitionProb = ArrayHandler.clone( statisticsTransitionProb );
		return clone;
	}
	
	/**
	 * This method determines the self transition.
	 * 
	 * @see #diagElement
	 */
	protected void determineDiagonalElement() {
		diagElement = -1;
		if( context.length != 0 ) {
			int last = context[context.length-1];
			diagElement = 0;
			while( diagElement < states.length && states[diagElement] != last ) {
				diagElement++;
			}
		}
		if( diagElement < 0 || diagElement >= states.length ) {
			throw new IllegalArgumentException( "Could not determine the diagonal element." );
		}
		
		//System.out.println( "DiaEl = " + diagElement );
	}
	
	protected void precompute() {
	}

	public double getLogPriorTerm()
	{
		if( this.ess == 0 ) {
			//Maximum-Likelihood estimation
			return( 0 );
		} else {
			double logPriorVal = 0;
			//Transition prior
			for( int i = 0; i < this.states.length; i++ )
			{
				//Lambda
				logPriorVal += Math.log( this.parameters[ i ] ) * this.hyperParameters[ i ];				
			}
			return( logPriorVal );
		}
	}
	
	public double getLogGammaScoreFromStatistic() {
		return Double.NaN;//XXX
	}
	
	//XXX not random
	public void initializeRandomly() {
		System.arraycopy( probabilities, 0, parameters, 0, states.length);
	}
	
	protected void appendFurtherInformation( StringBuffer xml ) {
		XMLParser.appendObjectWithTags( xml, annotationID, "annotationID" );
		XMLParser.appendObjectWithTags( xml, probabilities, "probabilities" );
		XMLParser.appendObjectWithTags( xml, ess, "ess" );
	}

	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		annotationID = (String) XMLParser.extractObjectForTags( xml, "annotationID" );
		probabilities = (double[]) XMLParser.extractObjectForTags( xml, "probabilities" );
		ess = (Double) XMLParser.extractObjectForTags( xml, "ess" );
	}
}
