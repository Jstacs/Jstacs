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

import de.jstacs.io.NonParsableException;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.BasicHigherOrderTransition.AbstractTransitionElement;

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
		this( context, states, hyperParameters, weight, true );
	}

	/**
	 * This is the main constructor creating a new instance with given context, descendant states, and hyper parameters.
	 * 
	 * @param context the context (=previously visited state indices); last entry corresponds to the last state visited
	 * @param states the transitions to all possible states; if <code>null</code> then no transition allowed
	 * @param hyperParameters the hyper parameters for the transitions; if <code>null</code> then no prior is used
	 * @param weight the weight for plotting the edges in Graphviz, enables to modify the edge length, larger weights imply shorter edges (default: 1)
	 * @param norm whether a normalized or unnormalized variant should be created
	 */
	public BasicTransitionElement( int[] context, int[] states, double[] hyperParameters, double[] weight, boolean norm ){
		super( context, states, hyperParameters, weight, norm );
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
