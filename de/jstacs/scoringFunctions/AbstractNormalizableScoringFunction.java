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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.scoringFunctions;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.data.DataSet.ElementEnumerator;

/**
 * This class is the main part of any {@link de.jstacs.classifier.scoringFunctionBased.ScoreClassifier}. It implements
 * many methods of the interface {@link NormalizableScoringFunction}.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public abstract class AbstractNormalizableScoringFunction extends AbstractScoringFunction implements NormalizableScoringFunction {

	/**
	 * The main constructor.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer} of this {@link ScoringFunction}
	 * @param length
	 *            the length of this {@link ScoringFunction}, i.e. the length of
	 *            the modeled sequences
	 *            
	 * @throws IllegalArgumentException            
	 *            if the length is negative or does not match with {@link AlphabetContainer#getPossibleLength()}
	 */
	public AbstractNormalizableScoringFunction( AlphabetContainer alphabets, int length ) throws IllegalArgumentException {
		super( alphabets, length );
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link AbstractNormalizableScoringFunction} out of a {@link StringBuffer}
	 * .
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public AbstractNormalizableScoringFunction( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.AbstractScoringFunction#clone()
	 */
	@Override
	public AbstractNormalizableScoringFunction clone() 	throws CloneNotSupportedException {
		return (AbstractNormalizableScoringFunction) super.clone();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.NormalizableScoringFunction#isNormalized()
	 */
	public boolean isNormalized() {
		return false;
	}

	/**
	 * This method checks whether all given {@link NormalizableScoringFunction}s
	 * are normalized.
	 * 
	 * @param function
	 *            the {@link NormalizableScoringFunction}s to be checked
	 * 
	 * @return <code>true</code> if all {@link NormalizableScoringFunction}s are
	 *         already normalized, otherwise <code>false</code>
	 * 
	 * @see NormalizableScoringFunction#isNormalized()
	 */
	public static boolean isNormalized(ScoringFunction... function) {
		int i = 0;
		while( i < function.length && (function[i] == null || (function[i] instanceof NormalizableScoringFunction && ((NormalizableScoringFunction)function[i]).isNormalized())) ) {
			i++;
		}
		return i == function.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.NormalizableScoringFunction#getInitialClassParam
	 * (double)
	 */
	public double getInitialClassParam( double classProb ) {
		return Math.log( classProb ) - getLogNormalizationConstant();
	}

	@Override
	public double getLogProbFor( Sequence sequence ) {
		return getLogScoreFor( sequence ) - getLogNormalizationConstant();
	}

	@Override
	public double getLogProbFor( Sequence sequence, int startpos ) throws Exception {
		return getLogScoreFor( sequence, startpos ) - getLogNormalizationConstant();
	}

	@Override
	public double getLogProbFor( Sequence sequence, int startpos, int endpos ) {
		return getLogScoreFor( sequence.getSubSequence( startpos, endpos-startpos+1 ) ) - getLogNormalizationConstant();
	}
	
	@Override
	public double[] getLogScoreFor( DataSet data ) throws Exception {
		double[] probs = new double[data.getNumberOfElements()];
		getLogScoreFor( data, probs );
		return probs;
	}

	@Override
	public void getLogScoreFor( DataSet data, double[] res ) throws Exception {
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
	 * @see de.jstacs.StatisticalModel#emitDataSet(int, int[])
	 */
	@Override
	public DataSet emitDataSet(int numberOfSequences, int... seqLength) throws NotTrainedException, Exception {
		throw new Exception( "Standard implementation of emitSample used for "
						+ getInstanceName()	+ ". You have to overwrite this method to use it in a proper way.");
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.StatisticalModel#getMaximalMarkovOrder()
	 */
	@Override
	public byte getMaximalMarkovOrder() throws UnsupportedOperationException {
		throw new UnsupportedOperationException( "The maximal markov order for this model in undefined.");
	}
}
