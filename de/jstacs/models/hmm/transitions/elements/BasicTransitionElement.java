package de.jstacs.models.hmm.transitions.elements;

import de.jstacs.NonParsableException;
import de.jstacs.models.hmm.transitions.BasicHigherOrderTransition.AbstractTransitionElement;

/**
 * This class implements the probability distribution for a given context, i.e. it contains all possible
 * transition and the corresponding probabilities for a given set of previously visited states.
 * This class just implements the basic functionalities.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class BasicTransitionElement extends AbstractTransitionElement {

	/**
	 * This is the main constructor creating a new instance with given context, descendant states, and hyper parameters.
	 * 
	 * @param context the context (=previously visited state indices); last entry corresponds to the last state visited
	 * @param states the transitions to all possible states; if <code>null</code> then no transition allowed
	 * @param hyperParameters the hyper parameters for the transitions; if <code>null</code> then no prior is used
	 * 
	 * @see #BasicTransitionElement(int[], int[], double[], double[])
	 */
	public BasicTransitionElement( int[] context, int[] states, double[] hyperParameters ){
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
	public BasicTransitionElement( int[] context, int[] states, double[] hyperParameters, double[] weight ){
		super( context, states, hyperParameters, weight );
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
	public BasicTransitionElement( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	protected void appendFurtherInformation( StringBuffer xml ) {
	}

	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
	}
}
