package de.jstacs;

import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.models.Model;


public interface StatisticalModel extends SequenceScoringFunction {

	public double getLogProbFor(Sequence sequence) throws Exception;
	
	public double getLogProbFor(Sequence sequence, int startpos) throws Exception;
	
	public double getLogProbFor(Sequence sequence, int startpos, int endpos) throws Exception;
	
	public double getESS();
	
	/**
	 * Returns a value that is proportional to the prior. For ML 1 should be
	 * returned.
	 * 
	 * @return a value that is proportional to the prior
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public double getPriorTerm() throws Exception;

	/**
	 * Returns a value that is proportional to the log of the prior. For maximum likelihood (ML) 0
	 * should be returned.
	 * 
	 * @return a value that is proportional to the log of the prior
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see Model#getPriorTerm()
	 */
	public double getLogPriorTerm() throws Exception;
	
	/**
	 * This method returns a {@link Sample} object containing artificial
	 * sequence(s).
	 * 
	 * <br>
	 * <br>
	 * 
	 * There are two different possibilities to create a sample for a model with
	 * length 0 (homogeneous models).
	 * <ol>
	 * <li> <code>emitSample( int n, int l )</code> should return a sample with
	 * <code>n</code> sequences of length <code>l</code>.
	 * <li> <code>emitSample( int n, int[] l )</code> should return a sample with
	 * <code>n</code> sequences which have a sequence length corresponding to
	 * the entry in the given array <code>l</code>.
	 * </ol>
	 * 
	 * <br>
	 * 
	 * There are two different possibilities to create a sample for a model with
	 * length greater than 0 (inhomogeneous models).<br>
	 * <code>emitSample( int n )</code> and
	 * <code>emitSample( int n, null )</code> should return a sample with
	 * <code>n</code> sequences of length of the model (
	 * {@link Model#getLength()}).
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
	 * @return a {@link Sample} containing the artificial sequence(s)
	 * 
	 * @throws Exception
	 *             if the emission did not succeed
	 * @throws NotTrainedException
	 *             if the model is not trained yet
	 * 
	 * @see Sample
	 */
	public Sample emitSample(int numberOfSequences, int... seqLength)
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
