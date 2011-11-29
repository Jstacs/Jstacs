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

package de.jstacs.scoringFunctions.homogeneous;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
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
public class UniformHomogeneousScoringFunction extends
		HomogeneousScoringFunction {
	private double ess, p, logP;

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link UniformHomogeneousScoringFunction} that models each sequence
	 * uniformly.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param ess
	 *            the equivalent sample size (ess) for the class
	 */
	public UniformHomogeneousScoringFunction( AlphabetContainer alphabets, double ess ) {
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
	 * {@link UniformHomogeneousScoringFunction} out of its XML representation
	 * as returned by {@link #fromXML(StringBuffer)}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public UniformHomogeneousScoringFunction(StringBuffer xml)
			throws NonParsableException {
		super(xml);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getInstanceName()
	 */
	public String getInstanceName() {
		return "uniform";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.VariableLengthScoringFunction#getLogScore(
	 * de.jstacs.data.Sequence, int, int)
	 */
	public double getLogScore(Sequence seq, int start, int length) {
		return length * logP;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.VariableLengthScoringFunction#
	 * getLogScoreAndPartialDerivation(de.jstacs.data.Sequence, int, int,
	 * de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	public double getLogScoreAndPartialDerivation(Sequence seq, int start,
			int length, IntList indices, DoubleList dList) {
		return length * logP;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getNumberOfParameters()
	 */
	public int getNumberOfParameters() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#setParameters(double[],
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
	 * @see de.jstacs.scoringFunctions.VariableLengthScoringFunction#
	 * getLogNormalizationConstant(int)
	 */
	public double getLogNormalizationConstant(int length) {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunction(int,
	 * boolean, de.jstacs.data.Sample[], double[][])
	 */
	public void initializeFunction(int index, boolean meila, Sample[] data,
			double[][] weights) {
		// does nothing
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#fromXML
	 * (java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		StringBuffer b = XMLParser.extractForTag(xml, getClass()
				.getSimpleName());
		length = XMLParser.extractObjectForTags(b, "length", int.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		alphabets = XMLParser.extractObjectForTags(b, "alphabets", AlphabetContainer.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		ess = XMLParser.extractObjectForTags(b, "ess", double.class );
		computeLogP();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * getSizeOfEventSpaceForRandomVariablesOfParameter(int)
	 */
	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index ) {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.VariableLengthScoringFunction#
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
	 * @see de.jstacs.scoringFunctions.NormalizableScoringFunction#getEss()
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
	 * de.jstacs.scoringFunctions.NormalizableScoringFunction#getLogPriorTerm()
	 */
	public double getLogPriorTerm() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * addGradientOfLogPriorTerm(double[], int)
	 */
	public void addGradientOfLogPriorTerm(double[] grad, int start) {
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#getCurrentParameterValues()
	 */
	public double[] getCurrentParameterValues() throws Exception {
		return new double[0];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#isInitialized()
	 */
	public boolean isInitialized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#isNormalized
	 * ()
	 */
	@Override
	public boolean isNormalized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.homogeneous.HomogeneousScoringFunction#
	 * getMaximalMarkovOrder()
	 */
	@Override
	public int getMaximalMarkovOrder() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#initializeFunctionRandomly
	 * (boolean)
	 */
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.VariableLengthScoringFunction#setStatisticForHyperparameters(int[], double[])
	 */
	public void setStatisticForHyperparameters(int[] length, double[] weights)
			throws Exception {
		// does nothing
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.homogeneous.HomogeneousScoringFunction#initializeUniformly(boolean)
	 */
	@Override
	public void initializeUniformly(boolean freeParams) {
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.SamplingScoringFunction#getSamplingGroups(int)
	 */
	@Override
	public int[][] getSamplingGroups(int parameterOffset) {
		return new int[0][];
	}
}
