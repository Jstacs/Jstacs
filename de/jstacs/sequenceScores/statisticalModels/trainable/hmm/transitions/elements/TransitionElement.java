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

import java.util.Arrays;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;

/**
 * This class implements a transition element used
 * for training via sampling or gradient based optimization approach.
 * 
 * @author Jens Keilwagen
 */
public class TransitionElement extends BasicTransitionElement {

	/**
	 * The precomputed probabilities for each possible transition.
	 * 
	 * @see #getLogScoreAndPartialDerivation(int, IntList, DoubleList, Sequence, int)
	 * @see #addGradientForLogPriorTerm(double[], int)
	 */
	protected double[] probs;
	
	/**
	 * The internally used parameter offset.
	 * 
	 * @see #setParameterOffset(int)
	 * @see #getLogScoreAndPartialDerivation(int, IntList, DoubleList, Sequence, int)
	 * @see #addGradientForLogPriorTerm(double[], int)
	 */
	protected int offset;
	
	/**
	 * This is the main constructor creating a new instance with given context, descendant states, and hyper parameters.
	 * 
	 * @param context the context (=previously visited state indices); last entry corresponds to the last state visited
	 * @param states the transitions to all possible states; if <code>null</code> than no transition allowed
	 * @param hyperParameters the hyper parameters for the transitions; if <code>null</code> than no prior is used
	 * 
	 * @see #TransitionElement(int[], int[], double[], double[])
	 */
	public TransitionElement( int[] context, int[] states, double[] hyperParameters ) {
		this( context, states, hyperParameters, null );
	}
	
	/**
	 * This is the main constructor creating a new instance with given context, descendant states, and hyper parameters.
	 * 
	 * @param context the context (=previously visited state indices); last entry corresponds to the last state visited
	 * @param states the transitions to all possible states; if <code>null</code> than no transition allowed
	 * @param hyperParameters the hyper parameters for the transitions; if <code>null</code> than no prior is used
	 * @param weight the weight for plotting the edge in Graphviz, enables to modify the edge length, larger weights imply shorter edges (default: 1)
	 */
	public TransitionElement( int[] context, int[] states, double[] hyperParameters, double[] weight ) {
		this( context, states, hyperParameters, weight, true );
	}
	
	/**
	 * This is the main constructor creating a new instance with given context, descendant states, and hyper parameters.
	 * 
	 * @param context the context (=previously visited state indices); last entry corresponds to the last state visited
	 * @param states the transitions to all possible states; if <code>null</code> than no transition allowed
	 * @param hyperParameters the hyper parameters for the transitions; if <code>null</code> than no prior is used
	 * @param weight the weight for plotting the edge in Graphviz, enables to modify the edge length, larger weights imply shorter edges (default: 1)
	 * @param norm whether a normalized or unnormalized variant should be created
	 */
	public TransitionElement( int[] context, int[] states, double[] hyperParameters, double[] weight, boolean norm ) {
		super( context, states, hyperParameters, weight, norm );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link TransitionElement} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link TransitionElement} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 *             
	 */
	public TransitionElement( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	public TransitionElement clone() throws CloneNotSupportedException {
		TransitionElement clone = (TransitionElement) super.clone();
		clone.probs = probs.clone();
		return clone;
	}
	
	protected void init() {
		probs = new double[parameters.length];
		super.init();
	}
	
	protected void precompute() {
		if( norm ) {
			logNorm = Normalisation.logSumNormalisation( parameters, 0, parameters.length, probs, 0 );
		} else { 
			logNorm=0;
			Arrays.fill(probs, 0);
		}
	}
	
	/**
	 * Returns the logarithmic score and fills lists with
	 * the indices and the partial derivations.
	 * 
	 * @param childIdx the index of the child to be visited
	 * @param indices
	 *            an {@link IntList} of indices, after method invocation the
	 *            list should contain the indices i where
	 *            \( \frac{\partial \log score(seq)}{\partial \lambda_i} \)
	 *            is not zero
	 * @param partialDer
	 *            a {@link DoubleList} of partial derivations, after method
	 *            invocation the list should contain the corresponding
	 *            \( \frac{\partial \log score(seq)}{\partial \lambda_i} \)
	 *            that are not zero
	 * @param sequence the {@link Sequence} that can be used to extract {@link de.jstacs.data.sequences.annotation.SequenceAnnotation}
	 * @param sequencePosition the position within the sequence
	 * 
	 * @return the logarithmic score
	 */
	public double getLogScoreAndPartialDerivation( int childIdx, IntList indices, DoubleList partialDer, Sequence sequence, int sequencePosition ) {
		if( parameters.length>1 ) {
			if( norm ) {
				for( int i = 0; i < parameters.length; i++) {
					indices.add( offset + i );
					partialDer.add( - probs[i] );
				}
			}
			
			indices.add( offset + childIdx );
			partialDer.add( 1 );
		}
		return parameters[childIdx] - logNorm;
	}
	
	/**
	 * This method sets the internal {@link #offset} used  for several methods (cf. see tags).
	 * 
	 * @param o the offset to be set
	 * 
	 * @return the new offset that can be used for the next {@link TransitionElement}
	 * 
	 * @see #offset
	 * @see #getLogScoreAndPartialDerivation(int, IntList, DoubleList, Sequence, int)
	 * @see #addGradientForLogPriorTerm(double[], int)
	 */
	public int setParameterOffset( int o ) {
		offset = o;
		if( parameters.length>1 ) {
			return offset + parameters.length;
		} else {
			return offset;
		}
	}
	
	/**
	 * This method fills the current parameters of this {@link TransitionElement} into the given array <code>params</code>
	 * starting at position <code>offset</code>.
	 * 
	 * @param params the array for filling the parameters
	 * @param offset the start index
	 * 
	 * @return the next start index for further parameter filling
	 */
	public int fillParameters( double[] params, int offset ) {
		if( parameters.length>1 ) {
			for( int i = 0; i < parameters.length; i++) {
				params[offset+i] = parameters[i];
			}
			return offset + parameters.length;
		} else {
			return offset;
		}
	}
	
	/**
	 * This method sets the internal parameters to the values of
	 * <code>params</code> beginning at index <code>start</code>.
	 * 
	 * @param params
	 *            the new parameters
	 * @param start
	 *            the start index in <code>params</code>
	 *            
	 * @return the first index in <code>params</code> that has not been used
	 */
	public int setParameters( double[] params, int start ) {
		if( parameters.length > 1 ) {
			for( int i = 0; i < parameters.length; i++) {
				parameters[i] = params[start + i];
			}
			precompute();
			return start + parameters.length;
		} else {
			return start;
		}
	}
	
	/**
	 * This method computes the gradient of {@link TransitionElement#getLogPriorTerm()} for each
	 * parameter of this transition element. The results are added to the array
	 * <code>gradient</code> using the index <code>start</code>.
	 * 
	 * @param gradient
	 *            the array of gradients
	 * @param start
	 *            the start index in the <code>gradient</code> array used to add
	 *            partial derivations for the parameters
	 * 
	 * @see #offset
	 */
	public void addGradientForLogPriorTerm( double[] gradient, int start ) {
		if( hyperParameters.length>1 ) {
			double sum = 0;
			for( int i = 0; i < hyperParameters.length; i++ ) {
				sum += hyperParameters[i];
			}
			for( int i = 0; i < hyperParameters.length; i++ ) {
				gradient[start + offset + i] += hyperParameters[i] - sum*probs[i];
			}
		}
	}
	
	/**
	 * This method returns the minimal hyper parameters of this {@link TransitionElement}.
	 * 
	 * @return the minimal hyper parameters of this {@link TransitionElement}
	 */
	public double getMinimalHyperparameter() {
		return ToolBox.min( hyperParameters );
	}
	
	/**
	 * This method computes the log posterior from the internal sufficient statistic.
	 *  
	 * @return the log posterior from the internal sufficient statistic
	 */
	public double getLogPosteriorFromStatistic() {
		double logPost =0;
		 for( int i = 0; i < parameters.length; i++ ) {
			 logPost += statistic[i] * (parameters[i]-logNorm);
         }
		return logPost;
	}
}