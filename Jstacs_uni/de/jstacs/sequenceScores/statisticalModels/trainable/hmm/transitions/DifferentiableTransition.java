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

package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions;

import java.util.LinkedList;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;


/**
 * This class declares methods that allow for optimizing the parameters numerically using the {@link de.jstacs.algorithms.optimization.Optimizer}.
 * 
 * @author Jens Keilwagen
 */
public interface DifferentiableTransition extends Transition {

	/**
	 * This method sets the internal offset of the parameter index.
	 * 
	 * @param offset the offset
	 * 
	 * @return the offset for the next parameters
	 */
	public int setParameterOffset( int offset );
	
	/**
	 * This method allows to set the parameters of the transition.
	 * 
	 * @param params the parameters
	 * @param start the (global) start position
	 */
	public void setParameters( double[] params, int start );
	
	/**
	 * This method allows to fill the parameters of the transition in a given array.
	 * 
	 * @param params the parameters
	 */
	public void fillParameters( double[] params );
	
	/**
	 * This method allows to compute the logarithm of the score and the gradient for a specific transition.
	 * 
	 * @param layer the layer of the matrix
	 * @param index the index encoding the context
	 * @param childIdx the index of the child
	 * @param indices a list for the parameter indices
	 * @param partDer a list for the partial derivations
	 * @param sequencePosition the position within the sequence
	 * @param sequence the sequence
	 * 
	 * @return the logarithm of the score
	 * 
	 * @see Transition#getLogScoreFor(int, int, int, Sequence, int)
	 */
	public double getLogScoreAndPartialDerivation(int layer, int index, int childIdx, IntList indices, DoubleList partDer, Sequence sequence, int sequencePosition );
	
	/**
	 * This method computes the gradient of {@link Transition#getLogPriorTerm()} for each
	 * parameter of this transition. The results are added to the array
	 * <code>gradient</code> beginning at index <code>start</code>.
	 * 
	 * @param gradient
	 *            the array of gradients
	 * @param start
	 *            the start index in the <code>gradient</code> array, where the
	 *            partial derivations for the parameters of this {@link Transition} shall be
	 *            entered
	 */
	public void addGradientForLogPriorTerm( double[] gradient, int start );

	/**
	 * Returns the size of the event space, i.e., the number of possible outcomes,
	 * for the random variable of parameter <code>index</code>
	 * @param index the index of the parameter
	 * @return the size of the event space
	 */
	public int getSizeOfEventSpace( int index );

	/**
	 * Adds the groups of indexes of those parameters of this transition that should be sampled
	 * together in one step of a grouped sampling procedure, each as an <code>int[]</code>, into <code>list</code>. 
	 * In most cases, one group should contain the
	 * parameters that are living on a common simplex, e.g. the parameters of one {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement}
	 * of this transition. The internal indexes of the parameters are incremeneted by an external <code>parameterOffset</code>
	 * @param parameterOffset the external parameter offset
	 * @param list the list of sampling groups
	 */
	public void fillSamplingGroups( int parameterOffset, LinkedList<int[]> list );
}