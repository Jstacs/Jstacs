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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions;

import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This interface declares all methods needed in an emission during a numerical optimization of HMM.
 * 
 * @author Jens Keilwagen
 */
public interface DifferentiableEmission extends Emission {
	
	/**
	 * Fills the current parameters in the global <code>params</code> array using the internal offset.
	 * 
	 * @param params the global parameter array of the HMM
	 * 
	 * @see #setParameterOffset(int)
	 */
	public void fillCurrentParameter( double[] params );

	/**
	 * This method sets the internal parameters using the given global parameter array, the global offset of the HMM and the internal offset.
	 * 
	 * @param params the global parameter array of the classifier
	 * @param offset the offset of the HMM
	 * 
	 * @see #setParameterOffset(int)
	 */
	public void setParameter( double[] params, int offset );

	/**
	 * This method sets the internal parameter offset and returns the new parameter offset for further use.
	 * 
	 * @param offset the offset to be set
	 * 
	 * @return the new parameter offset
	 */
	public int setParameterOffset( int offset );
	
	/** 
	 * Returns the logarithmic score for a {@link Sequence} beginning at
	 * position <code>start</code> in the {@link Sequence} and fills lists with
	 * the indices and the partial derivations.
	 * 
	 * @param forward
	 *            a switch whether to use the forward or the reverse complementary strand of the sequence
	 * @param seq
	 *            the {@link Sequence}
	 * @param startPos
	 *            the start position in the {@link Sequence}
	 * @param endPos
	 *            the end position in the {@link Sequence}
	 * @param indices
	 *            an {@link IntList} of indices, after method invocation the
	 *            list should contain the indices i where
	 *            {@latex.inline $\\frac{\\partial \\log score(seq)}{\\partial \\lambda_i}$}
	 *            is not zero
	 * @param partDer
	 *            a {@link DoubleList} of partial derivations, after method
	 *            invocation the list should contain the corresponding
	 *            {@latex.inline $\\frac{\\partial \\log score(seq)}{\\partial \\lambda_i}$}
	 *            that are not zero
	 * 
	 * @return the logarithmic score for the {@link Sequence}
	 * 
	 * @throws OperationNotSupportedException if <code>forward==false</code> and the reverse complement of the sequence can not be computed 
	 */
	public double getLogProbAndPartialDerivationFor( boolean forward, int startPos, int endPos, IntList indices, DoubleList partDer, Sequence seq ) throws OperationNotSupportedException;

	/**
	 * This method computes the gradient of {@link #getLogPriorTerm()} for each
	 * parameter of this model. The results are added to the array
	 * <code>grad</code> beginning at index (<code>offset</code> + internal offset).
	 * 
	 * @param grad
	 *            the array of gradients
	 * @param offset
	 *            the start index of the HMM in the <code>grad</code> array, where the
	 *            partial derivations for the parameters of the HMM shall be
	 *            entered
	 * 
	 * @see Emission#getLogPriorTerm()
	 * @see #setParameterOffset(int)
	 */
	public void addGradientOfLogPriorTerm( double[] grad, int offset );
	
	/**
	 * Adds the groups of indexes of those parameters of this emission that should be sampled
	 * together in one step of a grouped sampling procedure, each as an <code>int[]</code>, into <code>list</code>. 
	 * In most cases, one group should contain the parameters that are living on a common simplex.
	 * The internal indexes of the parameters are incremeneted by an external <code>parameterOffset</code>
	 * @param parameterOffset the external parameter offset
	 * @param list the list of sampling groups
	 */
	public void fillSamplingGroups( int parameterOffset, LinkedList<int[]> list );
	
	/**
	 * Returns the number of parameters of this emission.
	 * @return the number of parameters
	 */
	public int getNumberOfParameters();

	/**
	 * Returns the size of the event space, i.e., the number of possible outcomes,
	 * for the random variables of this emission
	 * @return the size of the event space
	 */
	public int getSizeOfEventSpace();

}