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

package de.jstacs.sequenceScores;

import de.jstacs.Storable;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.ResultSet;

/**
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public interface SequenceScore extends Cloneable, Storable {
	/**
	 * Creates a clone (deep copy) of the current {@link SequenceScore} instance.
	 * 
	 * @return the cloned instance
	 * 
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 */
	public SequenceScore clone() throws CloneNotSupportedException;
	
	/**
	 * Returns the container of alphabets that were used when constructing the instance.
	 * 
	 * @return the container of alphabets that were used when constructing the instance
	 */
	public AlphabetContainer getAlphabetContainer();

	/**
	 * Should return a <b>short</b> instance name such as iMM(0), BN(2), ...
	 * 
	 * @return a short instance name
	 */
	public String getInstanceName();

	/**
	 * Returns the length of sequences this instance can score. Instances that can
	 * only score sequences of defined length are e.g. PWM or inhomogeneous
	 * Markov models. If the instance can score sequences of arbitrary length,
	 * e.g. homogeneous Markov models, this method returns 0 (zero).
	 * 
	 * @return the length of sequences the instance can score
	 */
	public int getLength();
	
	/**
	 * Returns some information characterizing or describing the current
	 * instance. This could be e.g. the number of edges for a
	 * Bayesian network or an image showing some representation of the instance.
	 * The set of characteristics should always include the XML-representation
	 * of the instance. The corresponding result type is
	 * {@link de.jstacs.results.StorableResult}.
	 * 
	 * @return the characteristics of the current instance 
	 * 
	 * @throws Exception
	 *             if some of the characteristics could not be defined
	 * 
	 * @see de.jstacs.results.StorableResult
	 */
	public ResultSet getCharacteristics() throws Exception;

	/**
	 * Returns the subset of numerical values that are also returned by
	 * {@link #getCharacteristics()}.
	 * 
	 * @return the numerical characteristics of the current instance
	 * 
	 * @throws Exception
	 *             if some of the characteristics could not be defined
	 */
	public NumericalResultSet getNumericalCharacteristics() throws Exception;
	
	/**
	 * Returns the logarithmic score for the {@link Sequence} <code>seq</code>.
	 * 
	 * @param seq
	 *            the sequence
	 * 
	 * @return the logarithmic score for the sequence
	 */
	public double getLogScoreFor(Sequence seq);

	/**
	 * Returns the logarithmic score for the {@link Sequence} <code>seq</code>
	 * beginning at position <code>start</code> in the {@link Sequence}.
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * @param start
	 *            the start position in the {@link Sequence}
	 * 
	 * @return the logarithmic score for the {@link Sequence}
	 */
	public double getLogScoreFor(Sequence seq, int start);
	
	/**
	 * Returns the logarithmic score for the {@link Sequence} <code>seq</code>
	 * beginning at position <code>start</code> in the {@link Sequence}.
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * @param start
	 *            the start position in the {@link Sequence}
	 * @param end
	 *            the end position (inclusive) in the {@link Sequence}
	 * 
	 * @return the logarithmic score for the {@link Sequence}
	 * 
	 * @throws Exception if, for instance, the subsequence length can not be handled
	 */
	public double getLogScoreFor(Sequence seq, int start, int end) throws Exception;
	
	/**
	 * This method computes the logarithm of the scores of all sequences
	 * in the given sample. The values are stored in an array according to the
	 * index of the respective sequence in the sample.
	 * 
	 * <br>
	 * <br>
	 * 
	 * The score for any sequence shall be computed independent of all
	 * other sequences in the sample. So the result should be exactly the same
	 * as for the method {@link #getLogScoreFor(Sequence)}.
	 * 
	 * @param data
	 *            the sample of sequences
	 * 
	 * @return an array containing the logarithm of the score of all
	 *         sequences of the sample
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see #getLogScoreFor(Sequence)
	 */
	public double[] getLogScoreFor(DataSet data) throws Exception;

	/**
	 * This method computes and stores the logarithm of the scores for
	 * any sequence in the sample in the given <code>double</code>-array.
	 * 
	 * <br>
	 * <br>
	 * 
	 * The score for any sequence shall be computed independent of all
	 * other sequences in the sample. So the result should be exactly the same
	 * as for the method {@link #getLogScoreFor(Sequence)}.
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
	 * @see #getLogScoreFor(Sequence)
	 * @see #getLogScoreFor(DataSet)
	 */
	public void getLogScoreFor(DataSet data, double[] res) throws Exception;

	/**
	 * This method can be used to determine whether the instance is initialized. If
	 * the instance is initialized you should be able to invoke {@link #getLogScoreFor(Sequence)}.
	 * 
	 * @return <code>true</code> if the instance is initialized, <code>false</code>
	 *         otherwise
	 */
	public boolean isInitialized();
}
