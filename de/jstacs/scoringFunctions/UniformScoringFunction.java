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
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This {@link ScoringFunction} does nothing. So it is possible to save
 * parameters in an optimization.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class UniformScoringFunction extends BasicUniformScoringFunction implements SamplingScoringFunction {
	private double ess, logP;

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link UniformScoringFunction} that models each sequence uniformly.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param length
	 *            the length of the modeled sequences
	 * @param ess
	 *            the equivalent sample size (ess) of the class
	 */
	public UniformScoringFunction(AlphabetContainer alphabets, int length,
			double ess) {
		super(alphabets, length);
		if (!alphabets.isDiscrete()) {
			throw new IllegalArgumentException( "The given AlphabetContainer has to be discrete." );
		}
		if (ess < 0) {
			throw new IllegalArgumentException( "The given ess has to be non-negative." );
		}
		this.ess = ess;
		computeLogP();
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link UniformScoringFunction} out of its XML
	 * representation as returned by {@link #fromXML(StringBuffer)}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public UniformScoringFunction(StringBuffer xml) throws NonParsableException {
		super( xml );
		computeLogP();
	}

	private void computeLogP() {
		logP = 1;
		for (int i = 0; i < length; i++) {
			logP *= alphabets.getAlphabetLengthAt(i);
		}
		logP = -Math.log(logP);
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.BasicUniformScoringFunction#getLogScore(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		return logP;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.BasicUniformScoringFunction#getLogScoreAndPartialDerivation(de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList dList) {
		return logP;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.BasicUniformScoringFunction#getFurtherInformation()
	 */
	@Override
	protected StringBuffer getFurtherInformation() {
		StringBuffer b = new StringBuffer(1000);
		XMLParser.appendObjectWithTags(b, ess, "ess");
		return b;
	}

	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.BasicUniformScoringFunction#extractFurtherInformation(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		try {
			ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
		} catch( NonParsableException n ){
			ess = 0;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * getLogNormalizationConstant()
	 */
	public double getLogNormalizationConstant() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunction(int,
	 * boolean, de.jstacs.data.Sample[], double[][])
	 */
	public void initializeFunction(int index, boolean meila, DataSet[] data,
			double[][] weights) {
		// does nothing
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * getSizeOfEventSpaceForRandomVariablesOfParameter(int)
	 */
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * getLogPartialNormalizationConstant(int)
	 */
	public double getLogPartialNormalizationConstant(int parameterIndex) throws Exception {
		throw new IndexOutOfBoundsException(
				"Since a uniform scoring function has no parameters, this method can not be used");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.NormalizableScoringFunction#getEss()
	 */
	public double getESS() {
		return ess;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.BasicUniformScoringFunction#toString()
	 */
	@Override
	public String toString() {
		StringBuffer info = new StringBuffer(length * 100);
		double val;
		for (int j = 0; j < length; j++) {
			val = 1d / alphabets.getAlphabetLengthAt(0);
			info.append(j + "\t" + val + " for each element of "
					+ alphabets.getAlphabetAt(j).toString());
			if (j < length - 1) {
				info.append("\n");
			}
		}
		return info.toString();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.NormalizableScoringFunction#getLogPriorTerm()
	 */
	public double getLogPriorTerm() {
		// since the normalization constant does not depend on any parameter,
		// it is constant and therefore left out
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.NormalizableScoringFunction#addGradientOfLogPriorTerm(double[], int)
	 */
	public void addGradientOfLogPriorTerm(double[] grad, int start) {
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#isNormalized()
	 */
	@Override
	public boolean isNormalized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.SamplingScoringFunction#getSamplingGroups(int)
	 */
	@Override
	public int[][] getSamplingGroups(int parameterOffset) {
		return new int[0][];
	}
	
	@Override
	public double getLogProbFor( Sequence sequence, int startpos ) throws Exception {
		return getLogScoreFor( sequence, startpos ) - getLogNormalizationConstant();
	}

	@Override
	public double getLogProbFor( Sequence sequence ) throws Exception {
		return getLogScoreFor( sequence ) - getLogNormalizationConstant();
	}

	@Override
	public double getLogProbFor(Sequence sequence, int startpos, int endpos) throws Exception {
		if( endpos-startpos+1 != length ) {
			throw new IllegalArgumentException();
		} else {
			return getLogProbFor( sequence, startpos );
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