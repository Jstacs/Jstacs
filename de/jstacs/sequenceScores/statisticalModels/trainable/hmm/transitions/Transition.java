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

package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions;

import java.text.NumberFormat;

import de.jstacs.Storable;
import de.jstacs.data.sequences.Sequence;

/**
 * This interface declares the methods of the transition used in a hidden Markov model.
 *  
 * <a id="sequence">
 * The interface declares the method {@link #getLogScoreFor(int, int, int, Sequence, int)}, which
 * surprisingly has two additional parameter, namely the actual position and the sequence. These
 * additional parameters allow to switch between different transition matrices. For more details,
 * we refer, for instance, to {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.ReferenceBasedTransitionElement}.
 * </a>
 * 
 * @author Jan Grau, Jens Keilwagen, Michael Scharfe
 * 
 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM
 */
public interface Transition extends Cloneable, Storable {	
	/**
	 * This method returns a deep clone of the current instance.
	 * 
	 * @return a deep clone of the current instance.
	 * 
	 * @throws CloneNotSupportedException
	 *             if the instance could not be cloned
	 */
	public Transition clone() throws CloneNotSupportedException;
	
	/**
	 * This method returns the number of children states for given index, i.e. context, and
	 * a given layer of the matrix.
	 * 
	 * @param layer
	 *            the layer of the matrix
	 * @param index
	 *            the index encoding the context
	 *            
	 * @return the number of children states
	 */
	public int getNumberOfChildren( int layer, int index );
	
	
	/**
	 * This method returns the maximal used Markov order.
	 * 
	 * @return maximal used Markov order
	 */
	public int getMaximalMarkovOrder();

	
	/**
	 * Returns a value that is proportional to the log of the prior. For maximum
	 * likelihood (ML) 0 should be returned.
	 * 
	 * @return a value that is proportional to the log of the prior
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorTerm();
	
	/**
	 * This method returns the number of states underlying this transition
	 * instance.
	 * 
	 * @return the number of states underlying this transition instance
	 */
	public int getNumberOfStates();

	/**
	 * This method randomly initializes the parameters of the transition.
	 */
	public void initializeRandomly();
	
	/**
	 * This method returns a {@link String} representation of the structure that
	 * can be used in <i>Graphviz</i> to create an image.
	 * 
	 * @param nf
	 *            the {@link NumberFormat} used for the probabilities, if
	 *            <code>null</code> no probabilities will we written
	 * @param arrowOption
	 *            this parameter gives the possibility to set some arrow option 
	 * @param graphical represent transition probabilities as thickness of edges
	 * 			instead of textual output
	 * 
	 * @return {@link String} representation of the state
	 */
	public String getGraphizNetworkRepresentation(NumberFormat nf, String arrowOption, boolean graphical);

	
	/**
	 * This method fills all relevant information for a specific edge in a given container.
	 * 
	 * @param layer the layer in the matrix
	 * @param index the index encoding the context
	 * @param childIdx the index of the child that will be visited
	 * @param container
	 * 			<code>container[0]</code> the index of the state;
	 * 			<code>container[1]</code> the index encoding the new context;
	 * 			<code>container[2]</code> the distance for the new layer, i.e. 0 for the same layer (=silent state), and 1 for the next layer (=non-silent state)
	 */
	public void fillTransitionInformation( int layer, int index, int childIdx, int[] container );
	
	/**
	 * This method returns the logarithm of the score for the transition.
	 * 
	 * @param layer the layer in the matrix
	 * @param index the index encoding the context
	 * @param childIdx the index of the child that will be visited
	 * @param sequence the sequence
	 * @param sequencePosition the position within the sequence
	 * 
	 * @return the logarithm of the score for the transition
	 */
	public double getLogScoreFor( int layer, int index, int childIdx, Sequence sequence, int sequencePosition);

	/**
	 * This method computes the number of different indexes for a given layer of the matrix.
	 * 
	 * @param layer the layer of the matrix
	 * 
	 * @return the number of different indexes for a given position
	 */
	public int getNumberOfIndexes(int layer);

	/**
	 * This method answers the question whether the instance models any self
	 * transitions. It returns <code>true</code> if the current instance has
	 * any self transitions, otherwise <code>false</code>
	 * 
	 * @return <code>true</code> if the current instance has any self
	 *         transitions, otherwise <code>false</code>
	 */
	public boolean hasAnySelfTransitions();

	/**
	 * This method returns the maximal out degree of any context used in this transition instance.
	 * 
	 * @return maximal out degree of any context
	 */
	public int getMaximalInDegree();
	
	/**
	 * This method returns the maximal number of children for any context used in this transition instance.
	 * 
	 * @return maximal number of children for any context
	 */
	public int getMaximalNumberOfChildren();

	/**
	 * The method returns the index of the state of the context, if there is no context -1 is returned.
	 * 
	 * @param layer the layer in the matrix
	 * @param index the index encoding the context
	 * 
	 * @return the index of the last state of the context, if there is no context -1 is returned 
	 */
	public int getLastContextState( int layer, int index );

	/**
	 * This method returns the child index of the state, if this state is no child of the context -1 is returned
	 * 
	 * @param layer the layer in the matrix
	 * @param index the index encoding the context
	 * @param state the index of the state
	 * 
	 * @return the child index of the state, if this state is no child of the context -1 is returned
	 */
	public int getChildIdx( int layer, int index, int state );
	
	/**
	 * This method returns for each state whether it is absorbing or not.
	 * 
	 * @return for each state whether it is absorbing or not
	 */
	public boolean[] isAbsorbing();
	
	/**
	 * This method returns a {@link String} representation of the {@link Transition} using the given names of the states.
	 * 
	 * @param stateNames the names of the states, can be <code>null</code>
	 * @param nf the {@link NumberFormat} for the {@link String} representation of probabilities
	 * 
	 * @return a {@link String} representation of the {@link Transition} using the given names of the states
	 */
	public String toString( String[] stateNames, NumberFormat nf );
	
	/**
	 * Set values of parameters of the instance to the value of the parameters of the given instance.
	 * It can be assumed that the given instance and the current instance are from the same class.
	 * 
	 * This method might be used for instance in a multi-threaded optimization to broadcast the parameters. 
	 * 
	 * @param t the transition with the parameters to be set 
	 * 
	 * @throws IllegalArgumentException if the assumption about the same class for given and current instance is wrong
	 */
	public void setParameters( Transition t ) throws IllegalArgumentException;
}