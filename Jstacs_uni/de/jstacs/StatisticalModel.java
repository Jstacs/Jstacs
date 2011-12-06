package de.jstacs;

import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;


/**
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public interface StatisticalModel extends SequenceScoringFunction {

	/**
	 * Returns the logarithm of the probability of (a part of) the given
	 * sequence given the model. If at least one random variable is continuous
	 * the value of density function is returned.
	 * 
	 * <br>
	 * <br>
	 * 
	 * It extends the possibility given by the method
	 * {@link #getLogProbFor(Sequence, int)} by the fact, that the model could be
	 * e.g. homogeneous and therefore the length of the sequences, whose
	 * probability should be returned, is not fixed. Additionally, the end
	 * position of the part of the given sequence is given and the probability
	 * of the part from position <code>startpos</code> to <code>endpos</code>
	 * (inclusive) should be returned.
	 * 
	 * <br>
	 * 
	 * The <code>length</code> and the <code>alphabets</code> define the type of
	 * data that can be modeled and therefore both has to be checked.
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
	 */
	public double getLogProbFor(Sequence sequence, int startpos, int endpos) throws Exception;

	/**
	 * Returns the logarithm of the probability of (a part of) the given
	 * sequence given the model. If at least one random variable is continuous
	 * the value of density function is returned.
	 * 
	 * <br>
	 * <br>
	 * 
	 * If the length of the sequences, whose probability should be returned, is
	 * fixed (e.g. in a inhomogeneous model) and the given sequence is longer
	 * than their fixed length, the start position within the given sequence is
	 * given by <code>startpos</code>. E.g. the fixed length is 12. The length
	 * of the given sequence is 30 and the <code>startpos</code>=15 the logarithm
	 * of the probability of the part from position 15 to 26 (inclusive) given
	 * the model should be returned.
	 * 
	 * <br>
	 * 
	 * The <code>length</code> and the <code>alphabets</code> define the type of
	 * data that can be modeled and therefore both has to be checked.
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
	 * @see #getLogProbFor(Sequence, int, int)
	 */
	public double getLogProbFor(Sequence sequence, int startpos) throws Exception;

	/**
	 * Returns the logarithm of the probability of the given sequence given the
	 * model. If at least one random variable is continuous the value of density
	 * function is returned.
	 * 
	 * <br>
	 * <br>
	 * 
	 * The <code>length</code> and the <code>alphabets</code> define the type of
	 * data that can be modeled and therefore both has to be checked.
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
	 * @see #getLogProbFor(Sequence, int, int)
	 */
	public double getLogProbFor(Sequence sequence) throws Exception;
	
	/**
	 * Returns a value that is proportional to the log of the prior. For maximum likelihood (ML) 0
	 * should be returned.
	 * 
	 * @return a value that is proportional to the log of the prior
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public double getLogPriorTerm() throws Exception;
	
	/**
	 * This method returns a {@link DataSet} object containing artificial
	 * sequence(s).
	 * 
	 * <br>
	 * <br>
	 * 
	 * There are two different possibilities to create a sample for a model with
	 * length 0 (homogeneous models).
	 * <ol>
	 * <li> <code>emitDataSet( int n, int l )</code> should return a data set with
	 * <code>n</code> sequences of length <code>l</code>.
	 * <li> <code>emitDataSet( int n, int[] l )</code> should return a data set with
	 * <code>n</code> sequences which have a sequence length corresponding to
	 * the entry in the given array <code>l</code>.
	 * </ol>
	 * 
	 * <br>
	 * 
	 * There are two different possibilities to create a sample for a model with
	 * length greater than 0 (inhomogeneous models).<br>
	 * <code>emitDataSet( int n )</code> and
	 * <code>emitDataSet( int n, null )</code> should return a sample with
	 * <code>n</code> sequences of length of the model (
	 * {@link StatisticalModel#getLength()}).
	 * 
	 * <br>
	 * <br>
	 * 
	 * The standard implementation throws an {@link Exception}.
	 * 
	 * @param numberOfSequences
	 *            the number of sequences that should be contained in the
	 *            returned sample
	 * @param seqLength
	 *            the length of the sequences for a homogeneous model; for an
	 *            inhomogeneous model this parameter should be <code>null</code>
	 *            or an array of size 0.
	 * 
	 * @return a {@link DataSet} containing the artificial sequence(s)
	 * 
	 * @throws Exception
	 *             if the emission did not succeed
	 * @throws NotTrainedException
	 *             if the model is not trained yet
	 * 
	 * @see DataSet
	 */
	public DataSet emitDataSet(int numberOfSequences, int... seqLength)
			throws NotTrainedException, Exception;

	

	/**
	 * This method returns the maximal used Markov order, if possible.
	 * 
	 * @return maximal used Markov order
	 * 
	 * @throws UnsupportedOperationException
	 *             if the model can't give a proper answer
	 */
	public byte getMaximalMarkovOrder() throws UnsupportedOperationException;
	
}
