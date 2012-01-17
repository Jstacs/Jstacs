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

package de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This scoring function does nothing. So it is possible to save parameters in
 * an optimization.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class UniformHomogeneousDiffSM extends HomogeneousDiffSM {
	
	private double ess, p, logP;

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link UniformHomogeneousDiffSM} that models each sequence
	 * uniformly.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param ess
	 *            the equivalent sample size (ess) for the class
	 */
	public UniformHomogeneousDiffSM( AlphabetContainer alphabets, double ess ) {
		super(alphabets);
		if (ess < 0) {
			throw new IllegalArgumentException(
					"The given ess has to be non-negative.");
		}
		this.ess = ess;
		computeLogP();
	}

	private void computeLogP() {
		p = 1d / alphabets.getAlphabetLengthAt(0);
		logP = Math.log(p);
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link UniformHomogeneousDiffSM} out of its XML representation
	 * as returned by {@link #fromXML(StringBuffer)}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public UniformHomogeneousDiffSM(StringBuffer xml)
			throws NonParsableException {
		super(xml);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getInstanceName()
	 */
	public String getInstanceName() {
		return "uniform";
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractVariableLengthDiffSM#getLogScoreFor(int, de.jstacs.data.Sequence, int)
	 */
	public double getLogScoreFor( Sequence seq, int start, int end ) {
		return (end-start+1) * logP;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractVariableLengthDiffSM#getLogScoreAndPartialDerivation(int, de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public double getLogScoreAndPartialDerivation( Sequence seq, int start, int end, IntList indices, DoubleList dList) {
		return getLogScoreFor( seq, start, end );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getNumberOfParameters()
	 */
	public int getNumberOfParameters() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#setParameters(double[],
	 * int)
	 */
	public void setParameters(double[] params, int start) {
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer b = new StringBuffer(1000);
		XMLParser.appendObjectWithTags(b, length, "length");
		XMLParser.appendObjectWithTags(b, alphabets, "alphabets");
		XMLParser.appendObjectWithTags(b, ess, "ess");
		XMLParser.addTags(b, getClass().getSimpleName());
		return b;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM#
	 * getLogNormalizationConstant(int)
	 */
	public double getLogNormalizationConstant(int length) {
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
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel#fromXML
	 * (java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		StringBuffer b = XMLParser.extractForTag(xml, getClass().getSimpleName());
		length = XMLParser.extractObjectForTags(b, "length", int.class );
		alphabets = XMLParser.extractObjectForTags(b, "alphabets", AlphabetContainer.class );
		ess = XMLParser.extractObjectForTags(b, "ess", double.class );
		computeLogP();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#
	 * getSizeOfEventSpaceForRandomVariablesOfParameter(int)
	 */
	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index ) {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM#
	 * getLogPartialNormalizationConstant(int, int)
	 */
	public double getLogPartialNormalizationConstant(int parameterIndex, int length)
			throws Exception {
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
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return p + " for each element of "
				+ alphabets.getAlphabetAt(0).toString();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorTerm() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#addGradientOfLogPriorTerm(double[], int)
	 */
	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start) {
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getCurrentParameterValues()
	 */
	public double[] getCurrentParameterValues() throws Exception {
		return new double[0];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#isInitialized()
	 */
	public boolean isInitialized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel#isNormalized
	 * ()
	 */
	@Override
	public boolean isNormalized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM#getMaximalMarkovOrder()
	 */
	@Override
	public byte getMaximalMarkovOrder() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#initializeFunctionRandomly
	 * (boolean)
	 */
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM#setStatisticForHyperparameters(int[], double[])
	 */
	public void setStatisticForHyperparameters(int[] length, double[] weights)
			throws Exception {
		// does nothing
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM#initializeUniformly(boolean)
	 */
	@Override
	public void initializeUniformly(boolean freeParams) {
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.SamplingDifferentiableStatisticalModel#getSamplingGroups(int)
	 */
	@Override
	public int[][] getSamplingGroups(int parameterOffset) {
		return new int[0][];
	}
}
