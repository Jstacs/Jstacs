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

package de.jstacs.trainableStatisticalModels;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.Storable;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.data.DataSet.ElementEnumerator;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;

/**
 * Abstract class for a model for pattern recognition. <br>
 * For writing or reading a {@link StringBuffer} to or from a file (
 * {@link #fromXML(StringBuffer)}, {@link #toXML()}) you can use the class
 * {@link de.jstacs.io.FileManager}.
 * 
 * @see de.jstacs.io.FileManager
 * 
 * @author Andre Gohr, Jan Grau, Jens Keilwagen
 */
public abstract class AbstractTrainSM implements Cloneable, Storable, TrainableStatisticalModel {
	/**
	 * The length of the sequences the model can classify. For models that can
	 * take sequences of arbitrary length this value should be set to 0
	 */
	protected int length;

	/**
	 * The underlying alphabets
	 */
	protected AlphabetContainer alphabets;

	/**
	 * Constructor that sets the length of the model to <code>length</code> and
	 * the {@link AlphabetContainer} to <code>alphabets</code>.
	 * 
	 * <br>
	 * 
	 * The parameter <code>length</code> gives the length of the sequences the
	 * model can classify. Models that can only classify sequences of defined
	 * length are e.g. PWM or inhomogeneous Markov models. If the model can
	 * classify sequences of arbitrary length, e.g. homogeneous Markov models,
	 * this parameter must be set to 0 (zero).
	 * 
	 * <br>
	 * 
	 * The <code>length</code> and <code>alphabets</code> define the type of
	 * data that can be modeled and therefore both has to be checked before any
	 * evaluation (e.g. {@link #getLogScoreFor(Sequence)})
	 * 
	 * @param alphabets
	 *            the alphabets in an {@link AlphabetContainer}
	 * @param length
	 *            the length of the sequences a model can classify, 0 for
	 *            arbitrary length
	 */
	public AbstractTrainSM(AlphabetContainer alphabets, int length) {
		this.length = length;
		this.alphabets = alphabets;
		if (alphabets.getPossibleLength() > 0
				&& alphabets.getPossibleLength() != length) {
			throw new IllegalArgumentException(
					"The length and the alphabet container does not match.");
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link AbstractTrainSM} out of a {@link StringBuffer}.
	 * 
	 * @param stringBuff
	 *            the {@link StringBuffer} to be parsed
	 * 
	 * @throws NonParsableException
	 *             is thrown if the {@link StringBuffer} could not be parsed
	 */
	public AbstractTrainSM(StringBuffer stringBuff) throws NonParsableException {
		alphabets = null;
		length = -1;
		fromXML(stringBuff);
		if (alphabets == null) {
			throw new NonParsableException(
					"The alphabets were not set correctly.");
		}
		if (length < 0) {
			throw new NonParsableException("The length was not set correctly.");
		}
		if (alphabets.getPossibleLength() > 0
				&& alphabets.getPossibleLength() != length) {
			throw new IllegalArgumentException(
					"The length and the alphabet container doesnot not match.");
		}
	}

	/**
	 * Follows the conventions of {@link Object}'s <code>clone()</code>-method.
	 * 
	 * @return an object, that is a copy of the current {@link AbstractTrainSM}
	 *         (the member-{@link AlphabetContainer} isn't deeply cloned since
	 *         it is assumed to be immutable). The type of the returned object
	 *         is defined by the class <code>X</code> directly inherited from
	 *         {@link AbstractTrainSM}. Hence <code>X</code>'s
	 *         <code>clone()</code>-method should work as:<br>
	 *         1. <code>Object o = (X)super.clone();</code> <br>
	 *         2. all additional member variables of <code>o</code> defined by
	 *         <code>X</code> that are not of simple data-types like
	 *         <code>int</code>, <code>double</code>, ... have to be deeply
	 *         copied <br>
	 *         3. <code>return o</code>
	 */
	@Override
	public AbstractTrainSM clone() throws CloneNotSupportedException {
		return (AbstractTrainSM) super.clone();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#train(de.jstacs.data.Sample)
	 */
	public void train(DataSet data) throws Exception {
		train(data, null);
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getLogProbFor(de.jstacs.data.Sequence)
	 */
	public double getLogProbFor(Sequence sequence) throws Exception {
		return getLogProbFor(sequence, 0, sequence.getLength() - 1);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getLogProbFor(de.jstacs.data.Sequence, int)
	 */
	public double getLogProbFor(Sequence sequence, int startpos) throws Exception {
		if (length == 0) {
			return getLogProbFor(sequence, startpos, sequence.getLength() - 1);
		} else {
			return getLogProbFor(sequence, startpos, startpos + length - 1);
		}
	}
	
	/**
	 * This method checks all parameters before a probability can be computed for a sequence.
	 * Hence, should be used in {@link #getLogProbFor(Sequence, int, int)}.
	 * 
	 * @param sequence
	 *            the given sequence
	 * @param startpos
	 *            the start position within the given sequence
	 * @param endpos
	 *            the last position to be taken into account
	 * 
	 * @throws IllegalArgumentException
	 *             if the sequence could not be handled (e.g.
	 *             <code>startpos &gt; </code>, <code>endpos
	 *             &gt; sequence.length</code>, ...) by the model
	 * @throws NotTrainedException
	 *             if the model is not trained yet
	 */
	protected void check( Sequence sequence, int startpos, int endpos ) throws NotTrainedException, IllegalArgumentException {
		if( !isInitialized() ) {
			throw new NotTrainedException();
		} else if( !alphabets.checkConsistency( sequence.getAlphabetContainer().getSubContainer( startpos, endpos - startpos + 1 ) ) ) {
			throw new IllegalArgumentException( "This sequence is not possible with the given alphabet." );
		} else if( startpos < 0 ) {
			throw new IllegalArgumentException( "This startposition is impossible. Try: 0 <= startposition" );
		} else if( startpos > endpos || endpos >= sequence.getLength() ) {
			throw new IllegalArgumentException( "This endposition is impossible. Try: startposition <= endposition < sequence.length" );
		} else if( length != 0 &&  endpos - startpos + 1 != length ) {
			throw new IllegalArgumentException( "This sequence has not length " + length + "." );
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getLogScoreFor(de.jstacs.data.Sequence)
	 */
	public double getLogScoreFor(Sequence sequence) {
		return getLogScoreFor( sequence, 0 );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getLogScoreFor(de.jstacs.data.Sequence, int)
	 */
	public double getLogScoreFor(Sequence sequence, int startpos) {
		try {
			return getLogProbFor( sequence, startpos );
		} catch( Exception e ) {
			RuntimeException r = new RuntimeException();
			r.setStackTrace( e.getStackTrace() );
			throw r;
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getLogScoreFor(int, de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor( Sequence sequence, int startpos, int endpos) {
		try {
			return getLogProbFor( sequence, startpos, endpos );
		} catch( Exception e ) {
			RuntimeException r = new RuntimeException();
			r.setStackTrace( e.getStackTrace() );
			throw r;
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getLogScoreFor(de.jstacs.data.Sample)
	 */
	public double[] getLogScoreFor(DataSet data) throws Exception {
		double[] res = new double[data.getNumberOfElements()];
		getLogScoreFor(data, res);
		return res;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getLogScoreFor(de.jstacs.data.Sample, double[])
	 */
	public void getLogScoreFor(DataSet data, double[] res) throws Exception {
		if (res.length != data.getNumberOfElements()) {
			throw new IllegalArgumentException("The array has wrong dimension.");
		}
		ElementEnumerator ei = new ElementEnumerator(data);
		for (int i = 0; i < res.length; i++) {
			res[i] = getLogScoreFor(ei.nextElement());
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.StatisticalModel#emitSample(int, int[])
	 */
	public DataSet emitDataSet(int numberOfSequences, int... seqLength) throws NotTrainedException, Exception {
		throw new Exception( "Standard implementation of emitSample used for "
						+ getInstanceName()	+ ". You have to overwrite this method to use it in a proper way.");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getAlphabetContainer()
	 */
	public final AlphabetContainer getAlphabetContainer() {
		return alphabets;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getLength()
	 */
	public final int getLength() {
		return length;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.StatisticalModel#getMaximalMarkovOrder()
	 */
	public byte getMaximalMarkovOrder() throws UnsupportedOperationException {
		throw new UnsupportedOperationException( "The maximal markov order for this model in undefined.");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getCharacteristics()
	 */
	public ResultSet getCharacteristics() throws Exception {
		return new ResultSet(getNumericalCharacteristics().getResults(),
				new Result[] { new StorableResult("model", "the xml representation of the model", this) });
	}

	/**
	 * This method should only be used by the constructor that works on a
	 * {@link StringBuffer}. It is the counter part of {@link #toXML()}.
	 * 
	 * @param xml
	 *            the XML representation of the model
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} is not parsable or the
	 *             representation is conflicting
	 * 
	 * @see AbstractTrainSM#AbstractTrainSM(StringBuffer)
	 */
	protected abstract void fromXML(StringBuffer xml) throws NonParsableException;
}
