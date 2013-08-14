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

package de.jstacs.sequenceScores.differentiable;

import java.text.NumberFormat;
import java.util.Locale;
import java.util.Random;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.DataSet.ElementEnumerator;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class is the main part of any {@link de.jstacs.classifiers.differentiableSequenceScoreBased.ScoreClassifier}. It implements
 * many methods of the interface {@link DifferentiableSequenceScore}.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public abstract class AbstractDifferentiableSequenceScore implements DifferentiableSequenceScore {
	
	/**
	 * Returns the number of recommended starts in a numerical optimization.
	 * 
	 * @param score
	 *            the scoring functions
	 * 
	 * @return the number of recommended starts
	 * 
	 * @see DifferentiableSequenceScore#getNumberOfRecommendedStarts()
	 */
	public static final int getNumberOfStarts( DifferentiableSequenceScore[] score ) {
		int starts = score[0].getNumberOfRecommendedStarts();
		for( int i = 1; i < score.length; i++ ) {
			starts = Math.max( starts, score[i].getNumberOfRecommendedStarts() );
		}
		return starts;
	}
	
	/**
	 * This object can be used for drawing initial parameters.
	 * 
	 * @see AbstractDifferentiableSequenceScore#getCurrentParameterValues()
	 */
	protected static final Random r = new Random();

	/**
	 * The {@link AlphabetContainer} of this {@link AbstractDifferentiableSequenceScore}
	 * .
	 */
	protected AlphabetContainer alphabets;

	/**
	 * The length of the modeled sequences.
	 */
	protected int length;

	/**
	 * The main constructor.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer} of this {@link DifferentiableSequenceScore}
	 * @param length
	 *            the length of this {@link DifferentiableSequenceScore}, i.e. the length of
	 *            the modeled sequences
	 *            
	 * @throws IllegalArgumentException            
	 *            if the length is negative or does not match with {@link AlphabetContainer#getPossibleLength()}
	 */
	public AbstractDifferentiableSequenceScore( AlphabetContainer alphabets, int length ) throws IllegalArgumentException {
		this.alphabets = alphabets;
		int l = alphabets.getPossibleLength();
		if( length < 0 || (length != 0 && l != 0 && length != l ) ) {
			throw new IllegalArgumentException( "The given length could not be used. The length has to be not negative and has to match with the possible length defined by the AlphabetContainer." );
		} else {
			this.length = length;
		}
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link AbstractDifferentiableSequenceScore} out of a {@link StringBuffer}
	 * .
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public AbstractDifferentiableSequenceScore( StringBuffer xml ) throws NonParsableException {
		alphabets = null;
		length = -1;
		fromXML(xml);
		if (alphabets == null) {
			throw new NonParsableException( "AlphabetContainer could not be parsed." );
		}
		if (length < 0 || alphabets == null) {
			throw new NonParsableException("Length could not be parsed.");
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public AbstractDifferentiableSequenceScore clone() throws CloneNotSupportedException {
		return (AbstractDifferentiableSequenceScore) super.clone();
	}

	/**
	 * This method is called in the constructor for the {@link de.jstacs.Storable}
	 * interface to create a scoring function from a {@link StringBuffer}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} could not be parsed
	 * 
	 * @see AbstractDifferentiableSequenceScore#AbstractDifferentiableSequenceScore(StringBuffer)
	 */
	protected abstract void fromXML(StringBuffer xml) throws NonParsableException;

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getAlphabetContainer()
	 */
	public final AlphabetContainer getAlphabetContainer() {
		return alphabets;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLength()
	 */
	public int getLength() {
		return length;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScore(de.jstacs.data.Sequence)
	 */
	@Override
	public final double getLogScoreFor( Sequence seq ) {
		return getLogScoreFor( seq, 0 );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScore(de.jstacs.data.Sequence)
	 */
	@Override
	public double getLogScoreFor( Sequence seq, int startpos, int endpos ) throws WrongLengthException {
		if( endpos-startpos+1 == length ) {
			return getLogScoreFor( seq, startpos );
		} else {
			throw new WrongLengthException( "The sequence length of " + (endpos-startpos+1) + " is not supported by this instance." );
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScoreAndPartialDerivation
	 * (de.jstacs.data.Sequence, de.jstacs.utils.IntList,
	 * de.jstacs.utils.DoubleList)
	 */
	public final double getLogScoreAndPartialDerivation( Sequence seq, IntList indices, DoubleList partialDer ) {
		return getLogScoreAndPartialDerivation( seq, 0, indices, partialDer );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScoreAndPartialDerivation(int, de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public double getLogScoreAndPartialDerivation( Sequence seq, int startpos, int endpos, IntList indices, DoubleList partialDer ) throws WrongLengthException {
		if( endpos-startpos+1 == length ) {
			return getLogScoreAndPartialDerivation( seq, startpos, indices, partialDer );
		} else {
			throw new WrongLengthException( "The sequence length of " + (endpos-startpos+1) + " is not supported by this instance." );
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getNumberOfRecommendedStarts()
	 */
	public int getNumberOfRecommendedStarts() {
		return 1;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getInitialClassParam(double)
	 */
	public double getInitialClassParam( double classProb ) {
		return Math.log( classProb );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getLogScoreFor(de.jstacs.data.DataSet)
	 */
	public double[] getLogScoreFor(DataSet data) throws Exception {
		double[] res = new double[data.getNumberOfElements()];
		getLogScoreFor(data, res);
		return res;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getLogScoreFor(de.jstacs.data.DataSet, double[])
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
	 * @see de.jstacs.SequenceScore#getCharacteristics()
	 */
	public ResultSet getCharacteristics() throws Exception {
		return new ResultSet(getNumericalCharacteristics().getResults(),
				new Result[] { new StorableResult("model", "the xml representation of the model", this) });
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#getNumericalCharacteristics()
	 */
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		return new NumericalResultSet( new NumericalResult( "number of parameters", "the number of parameters used in this instance to score sequences", getNumberOfParameters() ) );
	}
	
	public final String toString() {
		NumberFormat nf = NumberFormat.getInstance(Locale.US);
		nf.setMaximumFractionDigits(3);
		return toString(nf);
	}
}
