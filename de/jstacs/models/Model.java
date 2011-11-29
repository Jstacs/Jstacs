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

package de.jstacs.models;

import de.jstacs.NotTrainedException;
import de.jstacs.SequenceScoringFunction;
import de.jstacs.StatisticalModel;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;

/**
 * This interface defines all methods for a probabilistic model.
 * 
 * @author Andre Gohr, Jan Grau, Jens Keilwagen
 */
public interface Model extends StatisticalModel {

	/**
	 * Creates a clone (deep copy) of the current {@link Model} instance.
	 * 
	 * @return the cloned instance
	 * 
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 */
	public Model clone() throws CloneNotSupportedException;

	/**
	 * Trains the {@link Model} object given the data as {@link Sample}. <br>
	 * This method should work non-incrementally. That means the result of the
	 * following series: <code>train(data1)</code>; <code>train(data2)</code>
	 * should be a fully trained model over <code>data2</code> and not over
	 * <code>data1+data2</code>. All parameters of the model were given by the
	 * call of the constructor.
	 * 
	 * @param data
	 *            the given sequences as {@link Sample}
	 * @throws Exception
	 *             if the training did not succeed
	 * 
	 * @see Sample#getElementAt(int)
	 * @see de.jstacs.data.Sample.ElementEnumerator
	 */
	public void train(Sample data) throws Exception;

	/**
	 * Trains the {@link Model} object given the data as {@link Sample} using
	 * the specified weights. The weight at position i belongs to the element at
	 * position i. So the array <code>weight</code> should have the number of
	 * sequences in the sample as dimension. (Optionally it is possible to use
	 * <code>weight == null</code> if all weights have the value one.)<br>
	 * This method should work non-incrementally. That means the result of the
	 * following series: <code>train(data1)</code>; <code>train(data2)</code>
	 * should be a fully trained model over <code>data2</code> and not over
	 * <code>data1+data2</code>. All parameters of the model were given by the
	 * call of the constructor.
	 * 
	 * @param data
	 *            the given sequences as {@link Sample}
	 * @param weights
	 *            the weights of the elements, each weight should be
	 *            non-negative
	 * @throws Exception
	 *             if the training did not succeed (e.g. the dimension of
	 *             <code>weights</code> and the number of sequences in the
	 *             sample do not match)
	 * 
	 * @see Sample#getElementAt(int)
	 * @see de.jstacs.data.Sample.ElementEnumerator
	 */
	public void train(Sample data, double[] weights) throws Exception;

	/**
	 * Returns the logarithm of the probability of (a part of) the given
	 * sequence given the model. If at least one random variable is continuous
	 * the value of density function is returned.
	 * 
	 * <br>
	 * <br>
	 * 
	 * For more details see {@link Model#getProbFor(Sequence, int, int)}
	 * 
	 * @param sequence
	 *            the given sequence
	 * @param startpos
	 *            the start position within the given sequence
	 * @param endpos
	 *            the last position to be taken into account
	 * 
	 * @return the logarithm of the probability or the value of the density
	 *         function of (the part of) the given sequence given the model
	 * 
	 * @throws Exception
	 *             if the sequence could not be handled (e.g.
	 *             <code>startpos &gt; </code>, <code>endpos
	 *             &gt; sequence.length</code>, ...) by the model
	 * @throws NotTrainedException
	 *             if the model is not trained yet
	 * 
	 * @see Model#getProbFor(Sequence, int, int)
	 */
	public double getLogProbFor(Sequence sequence, int startpos, int endpos)
			throws Exception;

	/**
	 * Returns the logarithm of the probability of (a part of) the given
	 * sequence given the model. If at least one random variable is continuous
	 * the value of density function is returned.
	 * 
	 * <br>
	 * <br>
	 * 
	 * For more details see {@link Model#getProbFor(Sequence, int)}
	 * 
	 * @param sequence
	 *            the given sequence
	 * @param startpos
	 *            the start position within the given sequence
	 * 
	 * @return the logarithm of the probability or the value of the density
	 *         function of (the part of) the given sequence given the model
	 * 
	 * @throws Exception
	 *             if the sequence could not be handled by the model
	 * @throws NotTrainedException
	 *             if the model is not trained yet
	 * 
	 * @see Model#getProbFor(Sequence, int)
	 */
	public double getLogProbFor(Sequence sequence, int startpos)
			throws Exception;

	/**
	 * Returns the logarithm of the probability of the given sequence given the
	 * model. If at least one random variable is continuous the value of density
	 * function is returned.
	 * 
	 * <br>
	 * <br>
	 * 
	 * For more details see {@link Model#getProbFor(Sequence)}
	 * 
	 * @param sequence
	 *            the given sequence for which the logarithm of the
	 *            probability/the value of the density function should be
	 *            returned
	 * 
	 * @return the logarithm of the probability or the value of the density
	 *         function of the part of the given sequence given the model
	 * 
	 * @throws Exception
	 *             if the sequence could not be handled by the model
	 * @throws NotTrainedException
	 *             if the model is not trained yet
	 * 
	 * @see Model#getProbFor(Sequence)
	 */
	public double getLogProbFor(Sequence sequence) throws Exception;

	/**
	 * This method computes the logarithm of the probabilities of all sequences
	 * in the given sample. The values are stored in an array according to the
	 * index of the respective sequence in the sample.
	 * 
	 * <br>
	 * <br>
	 * 
	 * The probability for any sequence shall be computed independent of all
	 * other sequences in the sample. So the result should be exactly the same
	 * as for the method {@link #getLogProbFor(Sequence)}.
	 * 
	 * @param data
	 *            the sample of sequences
	 * 
	 * @return an array containing the logarithm of the probabilities of all
	 *         sequences of the sample
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see Model#getLogProbFor(Sequence)
	 */
	public double[] getLogProbFor(Sample data) throws Exception;

	/**
	 * This method computes and stores the logarithm of the probabilities for
	 * any sequence in the sample in the given <code>double</code>-array.
	 * 
	 * <br>
	 * <br>
	 * 
	 * The probability for any sequence shall be computed independent of all
	 * other sequences in the sample. So the result should be exactly the same
	 * as for the method {@link #getLogProbFor(Sequence)}.
	 * 
	 * @param data
	 *            the sample of sequences
	 * @param res
	 *            the array for the results, has to have length
	 *            <code>data.getNumberOfElements()</code> (which returns the
	 *            number of sequences in the sample)
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see Model#getLogProbFor(Sample)
	 */
	public void getLogProbFor(Sample data, double[] res) throws Exception;

	

	

	/**
	 * Returns <code>true</code> if the model has been trained successfully,
	 * <code>false</code> otherwise.
	 * 
	 * @return <code>true </code>if the model has been trained successfully,
	 *         <code>false</code> otherwise
	 */
	public boolean isTrained();


	/**
	 * Should give a simple representation (text) of the model as {@link String}.
	 * 
	 * @return the representation as {@link String}
	 */
	public String toString();

	/**
	 * This method tries to set a new instance of an {@link AlphabetContainer}
	 * for the current model. <b>This instance has to be consistent with the
	 * underlying instance of an {@link AlphabetContainer}.</b>
	 * 
	 * <br>
	 * <br>
	 * 
	 * This method can be very useful to save time.
	 * 
	 * @param abc
	 *            the alphabets in an {@link AlphabetContainer}
	 * 
	 * @return <code>true</code> if the new instance could be set
	 * 
	 * @see Model#getAlphabetContainer()
	 * @see AlphabetContainer#checkConsistency(AlphabetContainer)
	 */
	public boolean setNewAlphabetContainerInstance(AlphabetContainer abc);

}