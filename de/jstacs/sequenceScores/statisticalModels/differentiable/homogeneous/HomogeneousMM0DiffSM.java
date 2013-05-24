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

import java.text.NumberFormat;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.MEMConstraint;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.FastDirichletMRGParams;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * This scoring function implements a homogeneous Markov model of order zero
 * (hMM(0)) for a fixed sequence length.
 * 
 * @author Jens Keilwagen
 */
public class HomogeneousMM0DiffSM extends HomogeneousDiffSM {
	private double ess, norm, sumOfHyperParams, logGammaSum;
	private int[] counter;

	private boolean freeParams, plugIn, optimize;

	private MEMConstraint params;

	private int anz;

	/**
	 * The main constructor that creates an instance of a homogeneous Markov
	 * model of order 0.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer} of the model
	 * @param length
	 *            the length of sequences the model can handle
	 * @param ess
	 *            the equivalent sample size (ess)
	 * @param plugIn
	 *            indicates if a plug-in strategy to initialize the parameters
	 *            should be used
	 * @param optimize
	 *            indicates if the parameters should be optimized or not after
	 *            they have been initialized
	 */
	public HomogeneousMM0DiffSM(AlphabetContainer alphabets, int length,
			double ess, boolean plugIn, boolean optimize) {
		super(alphabets, length);
		if (ess < 0) {
			throw new IllegalArgumentException(
					"The ess has to be non-negative.");
		}
		this.ess = ess;
		sumOfHyperParams = ess * length;
		params = new MEMConstraint(new int[] { 0 }, new int[] { (int) alphabets
				.getAlphabetLengthAt(0) });
		this.plugIn = plugIn;
		this.optimize = optimize;
		setFreeParams(false);
		norm = 1;
		double d = -Math.log(alphabets.getAlphabetLengthAt(0));
		for (int i = 0; i < counter.length; i++) {
			params.setLambda(i, d);
		}
		computeConstantsOfLogPrior();
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link HomogeneousMM0DiffSM} out of its XML representation as returned by
	 * {@link #fromXML(StringBuffer)}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public HomogeneousMM0DiffSM(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel#clone()
	 */
	public HomogeneousMM0DiffSM clone() throws CloneNotSupportedException {
		HomogeneousMM0DiffSM clone = (HomogeneousMM0DiffSM) super.clone();
		clone.params = params.clone();
		clone.counter = counter.clone();
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getInstanceName()
	 */
	public String getInstanceName() {
		return "hMM(0)";
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractVariableLengthDiffSM#getLogScoreFor(int, de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor(Sequence seq, int start, int end) {
		double erg = 0;
		int length = end-start+1;
		for (int l = 0; l < length; l++) {
			erg += params.getLambda(params.satisfiesSpecificConstraint(seq,
					start + l));
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractVariableLengthDiffSM#getLogScoreAndPartialDerivation(int, de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, int end, IntList indices, DoubleList dList) {
		Arrays.fill(counter, 0);
		int l = 0, length = end-start+1;
		for (; l < length; l++) {
			counter[params.satisfiesSpecificConstraint(seq, start + l)]++;
		}
		double erg = 0;
		for (l = 0; l < counter.length; l++) {
			if (counter[l] > 0) {
				erg += counter[l] * params.getLambda(l);
				if (l < anz) {
					indices.add(l);
					dList.add(counter[l]);
				}
			}
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getNumberOfParameters()
	 */
	public int getNumberOfParameters() {
		return anz;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#setParameters(double[],
	 * int)
	 */
	public void setParameters(double[] params, int start) {
		if (optimize) {
			norm = 0;
			for (int i = 0; i < anz; i++) {
				this.params.setLambda(i, params[start + i]);
				norm += this.params.getExpLambda(i);
			}
			if (anz < counter.length) {
				norm += this.params.getExpLambda(anz);
			}
		}
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
		XMLParser.appendObjectWithTags(b, sumOfHyperParams, "sumOfHyperParams");
		XMLParser.appendObjectWithTags(b, params, "params");
		XMLParser.appendObjectWithTags(b, freeParams, "freeParams");
		XMLParser.appendObjectWithTags(b, plugIn, "plugIn");
		XMLParser.appendObjectWithTags(b, optimize, "optimize");
		XMLParser.addTags(b, getClass().getSimpleName());
		return b;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getCurrentParameterValues()
	 */
	public double[] getCurrentParameterValues() {
		double[] erg = new double[anz];
		for (int i = 0; i < anz; i++) {
			erg[i] = params.getLambda(i);
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#initializeFunction(int,
	 * boolean, de.jstacs.data.DataSet[], double[][])
	 */
	public void initializeFunction(int index, boolean freeParams,
			DataSet[] data, double[][] weights) {
		params.reset();
		if (plugIn) {
			if (data != null && data[index] != null) {
				Sequence seq;
				for (int k, l, i = 0; i < data[index].getNumberOfElements(); i++) {
					seq = data[index].getElementAt(i);
					l = seq.getLength();
					for (k = 0; k < l; k++) {
						params.add(seq.discreteVal(k), weights[index][i]);
					}
				}
			}
			params.estimate(sumOfHyperParams);
			for (int i = 0; i < counter.length; i++) {
				params.setExpLambda(i, params.getFreq(i));
			}
		} else {
			double d = -Math.log(alphabets.getAlphabetLengthAt(0));
			for (int i = 0; i < counter.length; i++) {
				params.setLambda(i, d);
			}
		}
		norm = 1;
		setFreeParams(freeParams);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#initializeFunctionRandomly
	 * (boolean)
	 */
	public void initializeFunctionRandomly(boolean freeParams) {
		int n = counter.length;
		double[] p = DirichletMRG.DEFAULT_INSTANCE.generate(n,
				new FastDirichletMRGParams(sumOfHyperParams == 0 ? 1
						: (sumOfHyperParams / (double) n)));
		for (int i = 0; i < n; i++) {
			params.setExpLambda(i, p[i]);
		}
		norm = 1;
		setFreeParams(freeParams);
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
		StringBuffer b = XMLParser.extractForTag(xml, getClass()
				.getSimpleName());
		length = XMLParser.extractObjectForTags(b, "length", int.class );
		alphabets = (AlphabetContainer) XMLParser.extractObjectForTags( b, "alphabets" );
		ess = XMLParser.extractObjectForTags(b, "ess", double.class );
		sumOfHyperParams = XMLParser.extractObjectForTags(b, "sumOfHyperParams", double.class );
		params = XMLParser.extractObjectForTags(b, "params", MEMConstraint.class );
		plugIn = XMLParser.extractObjectForTags(b, "plugIn", boolean.class );
		optimize = XMLParser.extractObjectForTags(b, "optimize", boolean.class );
		setFreeParams(XMLParser.extractObjectForTags(b, "freeParams", boolean.class ));
		for (int i = 0; i < params.getNumberOfSpecificConstraints(); i++) {
			norm += params.getExpLambda(i);
		}
		computeConstantsOfLogPrior();
	}

	private void setFreeParams(boolean freeParams) {
		this.freeParams = freeParams;
		counter = new int[params.getNumberOfSpecificConstraints()];
		if (optimize) {
			anz = counter.length - (freeParams ? 1 : 0);
		} else {
			anz = 0;
		}
		// TODO OK?
		if (freeParams) {
			double d = params.getLambda( params.getNumberOfSpecificConstraints() - 1 );
			for (int i = 0; i < params.getNumberOfSpecificConstraints(); i++) {
				params.setLambda(i, params.getLambda(i) - d );
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getSizeOfEventSpaceForRandomVariablesOfParameter(int)
	 */
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		if (index < anz) {
			return params.getNumberOfSpecificConstraints();
		} else {
			throw new IndexOutOfBoundsException();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM#
	 * getLogNormalizationConstant(int)
	 */
	public double getLogNormalizationConstant(int length) {
		if (length == 0) {
			throw new RuntimeException(
					"The normalization constant can not be computed for length 0.");
		} else {
			return norm * length;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM#
	 * getLogPartialNormalizationConstant(int, int)
	 */
	public double getLogPartialNormalizationConstant(int parameterIndex, int length)
			throws Exception {
		if (parameterIndex < anz) {
			return length + norm * (length - 1) + params.getLambda(parameterIndex);
		} else {
			throw new IndexOutOfBoundsException();
		}
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
	 * @see de.jstacs.sequenceScores.SequenceScore#toString(java.text.NumberFormat)
	 */
	@Override
	public String toString( NumberFormat nf ) {
		StringBuffer info = new StringBuffer(100);
		info.append(alphabets.getSymbol(0, 0) + ": "
				+ nf.format(params.getExpLambda(0) / norm));
		for (int i = 1; i < params.getNumberOfSpecificConstraints(); i++) {
			info.append("\t" + alphabets.getSymbol(0, i) + ": "
					+ nf.format(params.getExpLambda(i) / norm));
		}
		return info.toString();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorTerm() {
		if (optimize) {
			double val = 0;
			int n = params.getNumberOfSpecificConstraints(), i = 0;
			while (i < n) {
				val += params.getLambda(i++);
			}
			return (val * sumOfHyperParams / (double) n) + logGammaSum;
		}
		return 0;
	}

	private void computeConstantsOfLogPrior() {
		int anz = params.getNumberOfSpecificConstraints();
		logGammaSum = Gamma.logOfGamma(sumOfHyperParams) - anz
				* Gamma.logOfGamma(sumOfHyperParams / (double) anz);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#
	 * addGradientOfLogPriorTerm(double[], int)
	 */
	public void addGradientOfLogPriorTerm(double[] grad, int start) {
		double d = sumOfHyperParams
				/ (double) params.getNumberOfSpecificConstraints();
		for (int i = 0; i < anz; i++) {
			grad[start + i] += d;
		}
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
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM#
	 * getMaximalMarkovOrder()
	 */
	@Override
	public byte getMaximalMarkovOrder() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM#
	 * setStatisticForHyperparameters(int[], double[])
	 */
	public void setStatisticForHyperparameters(int[] length, double[] weight)
			throws Exception {
		if (weight.length != length.length) {
			throw new IllegalArgumentException(
					"The length of both arrays (length, weight) have to be identical.");
		}
		sumOfHyperParams = 0;
		for (int i = 0; i < length.length; i++) {
			if (weight[i] < 0 || length[i] < 0) {
				throw new IllegalArgumentException(
						"check length and weight for entry " + i);
			} else {
				sumOfHyperParams += length[i] * weight[i];
			}
		}
		computeConstantsOfLogPrior();
	}

	@Override
	public void initializeUniformly(boolean freeParams) {
		double p = 1d / (double) counter.length;
		for( int i = 0; i < counter.length; i++ ) {
			params.setExpLambda( i, p );
		}
		norm = 1;
		setFreeParams(freeParams);
	}

	@Override
	public int[][] getSamplingGroups( int parameterOffset ) {
		int[][] res = new int[1][params.getNumberOfSpecificConstraints()];
		for( int i = 0; i < res[0].length; i++ ) {
			res[0][i] = parameterOffset+i;
		}
		return res;
	}
}
