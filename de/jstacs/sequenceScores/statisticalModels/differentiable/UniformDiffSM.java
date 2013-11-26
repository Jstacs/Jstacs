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

package de.jstacs.sequenceScores.statisticalModels.differentiable;

import java.text.NumberFormat;

import de.jstacs.NotTrainedException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.differentiable.UniformDiffSS;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This {@link DifferentiableStatisticalModel} does nothing. So it is possible to save
 * parameters in an optimization.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class UniformDiffSM extends UniformDiffSS implements SamplingDifferentiableStatisticalModel {
	private double ess, logP;

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link UniformDiffSM} that models each sequence uniformly.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param length
	 *            the length of the modeled sequences
	 * @param ess
	 *            the equivalent sample size (ess) of the class
	 */
	public UniformDiffSM(AlphabetContainer alphabets, int length,
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
	 * Creates a new {@link UniformDiffSM} out of its XML
	 * representation as returned by {@link #fromXML(StringBuffer)}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public UniformDiffSM(StringBuffer xml) throws NonParsableException {
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
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.BasicUniformDiffSM#getLogScore(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		return logP;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.BasicUniformDiffSM#getLogScoreAndPartialDerivation(de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList dList) {
		return logP;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.BasicUniformDiffSM#getFurtherInformation()
	 */
	@Override
	protected StringBuffer getFurtherInformation() {
		StringBuffer b = new StringBuffer(1000);
		XMLParser.appendObjectWithTags(b, ess, "ess");
		return b;
	}

	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.BasicUniformDiffSM#extractFurtherInformation(java.lang.StringBuffer)
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
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogNormalizationConstant()
	 */
	@Override
	public double getLogNormalizationConstant() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#initializeFunction(int,
	 * boolean, de.jstacs.data.DataSet[], double[][])
	 */
	public void initializeFunction(int index, boolean meila, DataSet[] data,
			double[][] weights) {
		// does nothing
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getSizeOfEventSpaceForRandomVariablesOfParameter(int)
	 */
	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogPartialNormalizationConstant(int)
	 */
	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex) throws Exception {
		throw new IndexOutOfBoundsException(
				"Since a uniform scoring function has no parameters, this method can not be used");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getEss()
	 */
	public double getESS() {
		return ess;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.BasicUniformDiffSM#toString()
	 */
	@Override
	public String toString( NumberFormat nf ) {
		StringBuffer info = new StringBuffer(length * 100);
		double val;
		for (int j = 0; j < length; j++) {
			val = 1d / alphabets.getAlphabetLengthAt(0);
			info.append(j + "\t" + (nf==null?val:nf.format(val)) + " for each element of "
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
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorTerm() {
		// since the normalization constant does not depend on any parameter,
		// it is constant and therefore left out
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#addGradientOfLogPriorTerm(double[], int)
	 */
	public void addGradientOfLogPriorTerm(double[] grad, int start) {
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel#isNormalized()
	 */
	@Override
	public boolean isNormalized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.SamplingDifferentiableStatisticalModel#getSamplingGroups(int)
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
		throw new Exception( "Standard implementation of emitDataSet used for "
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