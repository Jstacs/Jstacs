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
import java.util.Random;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.ByteSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.StationaryDistribution;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.random.DiMRGParams;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;
import de.jstacs.utils.random.FastDirichletMRGParams;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * This scoring function implements a homogeneous Markov model of arbitrary
 * order for any sequence length. The scoring function uses the parameterization
 * of Meila if one uses the free parameters, which yields in a non-concave log
 * conditional likelihood.
 * 
 * @author Jens Keilwagen
 */
public class HomogeneousMMDiffSM extends HomogeneousDiffSM {
	private boolean freeParams, plugIn, optimize;

	private int order, starts;
	private int[] powers;
	private double classEss, logGammaSum;

	// the index can be computed using the powers array
	private double[][] params, probs, logNorm, hyperParams;//(order, (context, symbol))
	private int[] counter, distCounter, offset;

	/**
	 * This method returns an array that can be used in the constructor
	 * {@link #HomogeneousMMDiffSM(AlphabetContainer, int, double, double[], boolean, boolean, int)}
	 * containing the sums of the specific hyper-parameters.
	 * 
	 * @param order the order of the model
	 * @param length the sequence length
	 * @param ess the class ESS
	 * 
	 * @return an array containing the sums of the specific hyper-parameters
	 * 
	 * @see #HomogeneousMMDiffSM(AlphabetContainer, int, double, double[], boolean, boolean, int)
	 */
	public static double[] getSumOfHyperParameters( int order, int length, double ess ) {
		double[] sumOfHyperParams = new double[order+1];
		Arrays.fill( sumOfHyperParams, ess );
		sumOfHyperParams[order] = (length-order)*ess;
		return sumOfHyperParams;
	}
	
	/**
	 * This is a convenience constructor for creating an instance of a homogeneous
	 * Markov model of arbitrary order.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param order
	 *            the oder of the model (has to be non-negative)
	 * @param classEss
	 *            the equivalent sample size (ess) of the class
	 * @param length
	 *            the sequence length (only used for computing the hyper-parameters)
	 *            
	 * @see #getSumOfHyperParameters(int, int, double)
	 * @see #HomogeneousMMDiffSM(AlphabetContainer, int, double, double[], boolean, boolean, int)
	 */
	public HomogeneousMMDiffSM(AlphabetContainer alphabets, int order, double classEss, int length ) {
		this( alphabets, order, classEss, getSumOfHyperParameters( order, length, classEss ), true, true, 1 );
	}
	
	
	/**
	 * This is the main constructor that creates an instance of a homogeneous
	 * Markov model of arbitrary order.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param order
	 *            the oder of the model (has to be non-negative)
	 * @param classEss
	 *            the equivalent sample size (ess) of the class
	 * @param sumOfHyperParams
	 *            the sum of the hyper-parameters for each order (length has to
	 *            be <code>order</code>, each entry has to be non-negative)
	 * @param plugIn
	 *            a switch which enables to use the MAP-parameters as plug-in
	 *            parameters
	 * @param optimize
	 *            a switch which enables to optimize or fix the parameters
	 * @param starts
	 *            the number of recommended starts
	 */
	public HomogeneousMMDiffSM(AlphabetContainer alphabets, int order,
			double classEss, double[] sumOfHyperParams, boolean plugIn,
			boolean optimize, int starts) {
		this( alphabets, order, classEss, getHyperParams(alphabets, sumOfHyperParams), plugIn, optimize, starts );
	}
	
	private static double[][] getHyperParams( AlphabetContainer alphabets, double[] sumOfHyperParams ) {
		int al = (int) alphabets.getAlphabetLengthAt(0);
		int i, pow = 1;
		double[][] hyp = new double[sumOfHyperParams.length][];
		for (i = 0; i < sumOfHyperParams.length; i++) {
			pow *= al;
			hyp[i] = new double[pow];
			Arrays.fill( hyp[i], sumOfHyperParams[i]/(double)pow );
		}
		return hyp;
	}

	/**
	 * This is the main constructor that creates an instance of a homogeneous
	 * Markov model of arbitrary order with given hyper-parameters for the prior.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param order
	 *            the oder of the model (has to be non-negative)
	 * @param classEss
	 *            the equivalent sample size (ess) of the class
	 * @param hyperParams
	 *            the hyper-parameters for each order (length has to
	 *            be <code>order</code>, each entry has to be non-negative)
	 * @param plugIn
	 *            a switch which enables to use the MAP-parameters as plug-in
	 *            parameters
	 * @param optimize
	 *            a switch which enables to optimize or fix the parameters
	 * @param starts
	 *            the number of recommended starts
	 */
	public HomogeneousMMDiffSM(AlphabetContainer alphabets, int order,
			double classEss, double[][] hyperParams, boolean plugIn,
			boolean optimize, int starts) {
		super(alphabets);
		if (order < 0) {
			throw new IllegalArgumentException(
					"The order has to be non-negative.");
		}
		this.order = order;
		createArrays();
		if (classEss < 0) {
			throw new IllegalArgumentException(
					"The ess for the class has to be non-negative.");
		}
		this.classEss = classEss;
		if (hyperParams == null) {
			this.hyperParams = new double[order + 1][];
		} else {
			if (hyperParams.length != order + 1) {
				throw new IllegalArgumentException(
						"Wrong dimension of the ess array.");
			} else {
				this.hyperParams = new double[order + 1][];
				double[] sumOfHyperParams = new double[ order + 1 ];
				for (int i = 0; i <= order; i++) {
					if(hyperParams[i].length != powers[ i+1 ]){
						throw new IllegalArgumentException( "Wrong length of hyper-parameters for order "+i );
					}
					sumOfHyperParams[i] = 0;
					this.hyperParams[i] = new double[hyperParams[i].length];
					for(int j=0;j<powers[ i+1 ];j++){
						if (hyperParams[i][j] < 0) {
							throw new IllegalArgumentException(
									"The ess has to be non-negative. Violated at position "
									+ i + ".");
						}
						this.hyperParams[i][j] = hyperParams[i][j];
						sumOfHyperParams[i] += hyperParams[i][j];
					}
					if (i > 0 && i < order
							&& sumOfHyperParams[i] > sumOfHyperParams[i - 1]) {//TODO Double precision!!!
						throw new IllegalArgumentException(
								"The ess for start probabilities of order "
								+ i
								+ " is inconsistent with the ess for the probabilities of the previous order.");
					}
				}
			}
		}
		params = new double[order + 1][];
		double uniform = 1d / (double) powers[1], logUniform = Math
				.log(uniform);
		for (int i = 0; i <= order; i++) {
			params[i] = new double[powers[i + 1]];
			probs[i] = new double[powers[i + 1]];
			logNorm[i] = new double[powers[i]];
			Arrays.fill(params[i], logUniform);
			Arrays.fill(probs[i], uniform);
		}
		this.plugIn = plugIn;
		this.optimize = optimize;
		if (starts <= 0) {
			throw new IllegalArgumentException(
					"The number of starts has to be positive.");
		}
		this.starts = starts;
		setFreeParams(false);
		computeConstantsOfLogPrior();
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link HomogeneousMMDiffSM} out of its XML representation as returned by
	 * {@link #fromXML(StringBuffer)}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public HomogeneousMMDiffSM(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	private void createArrays() {
		powers = new int[order + 2];
		powers[0] = 1;
		powers[1] = (int) alphabets.getAlphabetLengthAt(0);
		int i;
		for (i = 2; i < powers.length; i++) {
			powers[i] = powers[i - 1] * powers[1];
		}
		probs = new double[order + 1][];
		logNorm = new double[order + 1][];
		counter = new int[powers[order + 1]];
		distCounter = new int[powers[order]];
		offset = new int[order + 2];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel#clone()
	 */
	@Override
	public HomogeneousMMDiffSM clone() throws CloneNotSupportedException {
		HomogeneousMMDiffSM clone = (HomogeneousMMDiffSM) super.clone();
		// the powers do not have to be cloned
		clone.params = new double[params.length][];
		clone.probs = new double[probs.length][];
		clone.logNorm = new double[logNorm.length][];
		for (int i = 0; i <= order; i++) {
			clone.params[i] = params[i].clone();
			clone.probs[i] = probs[i].clone();
			clone.logNorm[i] = logNorm[i].clone();
		}
		clone.counter = counter.clone();
		clone.distCounter = distCounter.clone();
		clone.offset = offset.clone();
		clone.hyperParams = ArrayHandler.clone(hyperParams);
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getInstanceName()
	 */
	public String getInstanceName() {
		return "hMM(" + order + ")";
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractVariableLengthDiffSM#getLogScoreFor(int, de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor( Sequence seq, int start, int end ) {
		double erg = 0;
		int length = end-start+1, l = 0, indexOld, indexNew = 0, o = Math.min(order, length);
		for (; l < o; l++) {
			indexOld = indexNew;
			indexNew = indexOld * powers[1] + seq.discreteVal(start++);
			erg += params[l][indexNew] - logNorm[l][indexOld];
		}
		for (; l < length; l++) {
			indexOld = indexNew % powers[order];
			indexNew = indexOld * powers[1] + seq.discreteVal(start++);
			erg += params[order][indexNew] - logNorm[order][indexOld];
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractVariableLengthDiffSM#getLogScoreAndPartialDerivation(int, de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public double getLogScoreAndPartialDerivation( Sequence seq, int start, int end, IntList indices, DoubleList dList) {
		if (optimize) {
			Arrays.fill(counter, 0);
			Arrays.fill(distCounter, 0);
			double erg = 0;
			int stop = powers[1] - (freeParams ? 1 : 0);
			int length = end-start+1, l = 0, indexOld, indexNew = 0, h, o = Math.min(order, length), index, z;
			// start probabilities
			for (; l < o; l++) {
				indexOld = indexNew;
				z = indexOld * powers[1];
				indexNew = z + seq.discreteVal(start++);
				erg += params[l][indexNew] - logNorm[l][indexOld];
				h = z - (freeParams ? indexOld : 0);
				for (index = 0; index < stop; index++) {
					indices.add(offset[l] + h + index);
					if (z + index == indexNew) {
						dList.add(1 - probs[l][z + index]);
					} else {
						dList.add(-probs[l][z + index]);
					}
				}
			}
			// counting the usage of transition probability parameters for the
			// sequence
			for (; l < length; l++) {
				indexOld = indexNew % powers[order];
				indexNew = indexOld * powers[1] + seq.discreteVal(start++);
				distCounter[indexOld]++;
				counter[indexNew]++;
			}
			// computing the gradient and the score
			for (l = 0; l < distCounter.length; l++) {
				if (distCounter[l] > 0) {
					h = l * (powers[1] - (freeParams ? 1 : 0));
					o = l * powers[1];
					for (index = 0; index < stop; index++, h++, o++) {
						indices.add(offset[order] + h);
						dList.add(counter[o] - distCounter[l]
										* probs[order][o]);
						erg += counter[o] * params[order][o];
					}
					if (stop < powers[1]) {
						erg += counter[o] * params[order][o];
					}
					erg -= distCounter[l] * logNorm[order][l];
				}
			}
			return erg;
		} else {
			return getLogScoreFor(seq, start, end);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getNumberOfParameters()
	 */
	public int getNumberOfParameters() {
		if (optimize) {
			return offset[order + 1];
		} else {
			return 0;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#setParameters(double[],
	 * int)
	 */
	public void setParameters(double[] params, int start) {
		if (optimize) {
			int stop = powers[1] - (freeParams ? 1 : 0);
			for (int j, n, index, o = 0; o <= order; o++) {
				for (index = n = 0; n < logNorm[o].length; n++) {
					logNorm[o][n] = 0d;
					for (j = 0; j < stop; j++, start++) {
						this.params[o][index + j] = params[start];
						probs[o][index + j] = Math
								.exp(this.params[o][index + j]);
						logNorm[o][n] += probs[o][index + j];
					}
					if (j < powers[1]) {
						probs[o][index + j] = Math
								.exp(this.params[o][index + j]);
						logNorm[o][n] += probs[o][index + j];
					}
					for (j = 0; j < powers[1]; j++, index++) {
						this.probs[o][index] /= logNorm[o][n];
					}
					logNorm[o][n] = Math.log(logNorm[o][n]);
				}
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
		XMLParser.appendObjectWithTags(b, order, "order");
		XMLParser.appendObjectWithTags(b, classEss, "classEss");
		XMLParser.appendObjectWithTags( b, hyperParams, "hyperParams" );	
		XMLParser.appendObjectWithTags(b, params, "params");
		XMLParser.appendObjectWithTags(b, plugIn, "plugIn");
		XMLParser.appendObjectWithTags(b, optimize, "optimize");
		XMLParser.appendObjectWithTags(b, starts, "starts");
		XMLParser.appendObjectWithTags(b, freeParams, "freeParams");
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
		int l = optimize ? offset[order + 1] : 0;
		double[] erg = new double[l];
		if (optimize) {
			int stop = powers[1] - (freeParams ? 1 : 0);
			for (int j, i = 0, index, o = 0; o <= order; o++) {
				for (index = 0; index < params[o].length; index += powers[1]) {
					for (j = 0; j < stop; j++, i++) {
						erg[i] = params[o][index + j];
					}
				}
			}
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
		int len, o, indexOld, indexNew;
		if (optimize && plugIn && data != null && data[index] != null) {
			double hyper;
			for (int i = 0; i <= order; i++) {
				System.arraycopy(probs[i], 0, hyperParams[i], 0, probs[i].length );
				for( int k = 0, j = 0; j < logNorm[i].length; j++ ) {
					logNorm[i][j] = 0;
					for( int l = 0; l < powers[1]; l++, k++ ) {
						logNorm[i][j] += probs[i][k];
					}
				}
			}
			Sequence seq;
			int anz = data[index].getNumberOfElements();
			double w = 1;
			boolean externalWeights = weights != null && weights[index] != null;
			// counting
			for (int l, i = 0; i < anz; i++) {
				seq = data[index].getElementAt(i);
				len = seq.getLength();
				o = Math.min(len, order);
				indexNew = 0;
				if (externalWeights) {
					w = weights[index][i];
				}
				for (l = 0; l < o; l++) {
					indexOld = indexNew;
					indexNew = indexOld * powers[1] + seq.discreteVal(l);
					probs[l][indexNew] += w;
					logNorm[l][indexOld] += w;
				}
				for (; l < len; l++) {
					indexOld = indexNew % powers[order];
					indexNew = indexOld * powers[1] + seq.discreteVal(l);
					probs[order][indexNew] += w;
					logNorm[order][indexOld] += w;
				}
			}
			// computing freqs and parameters
			for (o = 0; o <= order; o++) {
				for (indexOld = indexNew = 0; indexOld < logNorm[o].length; indexOld++) {
					if( logNorm[o][indexOld] > 0 ) {
						for (len = 0; len < powers[1]; len++, indexNew++) {
							probs[o][indexNew] /= logNorm[o][indexOld];
							params[o][indexNew] = Math.log(probs[o][indexNew]);
						}
						if(freeParams){
							int last = indexNew-1;
							indexNew -= powers[1];
							for (len = 0; len < powers[1]; len++, indexNew++) {
								params[o][indexNew] -= params[o][last];
							}
						}
					} else {
						hyper = 1d/(double) powers[1];
						for (len = 0; len < powers[1]; len++, indexNew++) {
							probs[o][indexNew] = hyper;
							params[o][indexNew] = 0;
						}
					}
					logNorm[o][indexOld] = 0;
				}
			}
		} else {
			initializeFunctionRandomly(freeParams);
		}
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
		if (optimize) {
			int o, normCounter, paramCounter, len;
			DiMRGParams hyper;
			double[] p = new double[powers[1]];
			double offset = 0, sum;
			for (o = 0; o <= order; o++) {
				for (normCounter = paramCounter = 0; normCounter < logNorm[o].length; normCounter++) {
					sum = 0;
					for (len = 0; len < powers[1]; len++) {
						p[len] = hyperParams[o][len+paramCounter];
						sum +=p[len];
					}
					hyper = sum == 0
							? new FastDirichletMRGParams(1)
							: new DirichletMRGParams(p);
					
					logNorm[o][normCounter] = 0;
					DirichletMRG.DEFAULT_INSTANCE.generate(p, 0, powers[1], hyper);
					if( freeParams ) {
						offset = Math.log(p[powers[1]-1]);
					}
					for (len = 0; len < powers[1]; len++, paramCounter++) {
						probs[o][paramCounter] = p[len];
						params[o][paramCounter] = Math.log(p[len]) - offset;
					}
				}
			}
			setFreeParams(freeParams);
		}
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
		alphabets = (AlphabetContainer) XMLParser.extractObjectForTags( b, "alphabets" );
		order = XMLParser.extractObjectForTags(b, "order", int.class );
		createArrays();
		classEss = XMLParser.extractObjectForTags(b, "classEss", double.class );
		if( XMLParser.hasTag(b, "hyperParams", null, null) ) {
			hyperParams = XMLParser.extractObjectForTags( b, "hyperParams", double[][].class );
		} else {
			hyperParams = getHyperParams(alphabets, XMLParser.extractObjectForTags(b, "sumOfHyperParams", double[].class ));
		}
		params = XMLParser.extractObjectForTags(b, "params", double[][].class );
		plugIn = XMLParser.extractObjectForTags(b, "plugIn", boolean.class );
		optimize = XMLParser.extractObjectForTags(b, "optimize", boolean.class );
		starts = XMLParser.extractObjectForTags(b, "starts", int.class );
		setFreeParams(XMLParser.extractObjectForTags(b, "freeParams", boolean.class ));
		for (int j, n, index, o = 0; o <= order; o++) {
			probs[o] = new double[params[o].length];
			logNorm[o] = new double[powers[o]];
			for (n = index = 0; n < logNorm[o].length; n++) {
				logNorm[o][n] = 0d;
				for (j = 0; j < powers[1]; j++) {
					probs[o][index + j] = Math.exp(this.params[o][index + j]);
					logNorm[o][n] += probs[o][index + j];
				}
				for (j = 0; j < powers[1]; j++, index++) {
					this.probs[o][index] /= logNorm[o][n];
				}
				logNorm[o][n] = Math.log(logNorm[o][n]);
			}
		}
		computeConstantsOfLogPrior();
	}

	private void setFreeParams(boolean freeParams) {
		this.freeParams = freeParams;
		if (optimize) {
			offset[0] = 0;
			for (int i = 0; i <= order; i++) {
				offset[i + 1] = offset[i] + params[i].length
						- (freeParams ? powers[i] : 0);
			}
		} else {
			offset[order + 1] = 0;
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getSizeOfEventSpaceForRandomVariablesOfParameter(int)
	 */
	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		int i = 0;
		while( index >= offset[i] ) {
			i++;
		}
		return powers[i];
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
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM#
	 * getLogPartialNormalizationConstant(int, int)
	 */
	public double getLogPartialNormalizationConstant(int parameterIndex, int length)
			throws Exception {
		if (parameterIndex < offset[order + 1]) {
			return Double.NEGATIVE_INFINITY;
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
		return classEss;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.SequenceScore#toString(java.text.NumberFormat)
	 */
	@Override
	public String toString( NumberFormat nf ) {
		StringBuffer info = new StringBuffer(100);
		DiscreteAlphabet abc = (DiscreteAlphabet) alphabets.getAlphabetAt(0);
		int i = 0, o, index, l = (int) abc.length();
		String[] sym = new String[l];
		l--;
		for (; i <= l; i++) {
			sym[i] = abc.getSymbolAt(i);
			info.append("\t" + sym[i]);
		}
		info.append("\n");
		int[] context = new int[order + 1];
		for (o = 0; o <= order; o++) {
			info.append("P(X_" + o);
			for (i = 0; i < o; i++) {
				if (i == 0) {
					info.append("|");
				} else {
					info.append(" ");
				}
				info.append("X_" + i);
			}
			info.append(")\n");
			Arrays.fill(context, 0);
			for (index = 0; index < probs[o].length;) {
				for (i = 0; i < o; i++) {
					info.append(sym[context[i]]);
				}
				for (i = 0; i <= l; i++, index++) {
					info.append("\t" + nf.format(probs[o][index]));
				}
				info.append("\n");

				i = o - 1;
				while (i >= 0 && context[i] == l) {
					context[i] = 0;
					i--;
				}
				if (i >= 0) {
					context[i]++;
				}
			}
			info.append("\n");
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
			for (int j, n, index, o = 0; o <= order; o++) {
				for (n = index = 0; n < logNorm[o].length; n++) {
					val -= ToolBox.sum( index, index+powers[1], hyperParams[o] ) * logNorm[o][n];
					for (j = 0; j < powers[1]; j++, index++) {
						val += hyperParams[o][index] * params[o][index];
					}
				}

			}
			return val + logGammaSum;
		} else {
			return 0;
		}
	}

	private void computeConstantsOfLogPrior() {
		double hyper;
		logGammaSum = 0;
		for (int o = 0; o <= order; o++) {
			for(int j=0;j<powers[ o + 1 ]; j+= powers[1]){
				double temp = ToolBox.sum( j, j+powers[1], hyperParams[o] );
				logGammaSum += Gamma.logOfGamma( temp );
				for(int k=0;k<powers[1];k++){
					hyper = hyperParams[o][j+k];
					if( hyper > 0 ){
						logGammaSum -= Gamma.logOfGamma(hyper);
					}
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#addGradientOfLogPriorTerm(double[], int)
	 */
	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start) {
		if (optimize) {
			double hyperSum;
			int stop = powers[1] - (freeParams ? 1 : 0);
			for (int j, index, o = 0; o <= order; o++) {
				
				for (index = 0; index < params[o].length; index += powers[1]) {
					hyperSum = ToolBox.sum( index, index+powers[1], hyperParams[o] );
					for (j = 0; j < stop; j++, start++) {
						grad[start] += hyperParams[o][index + j] - hyperSum * probs[o][index + j];
					}
				}
			}
			// System.out.println( start );System.exit(0);
		}
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
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#isInitialized()
	 */
	public boolean isInitialized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM#getMaximalMarkovOrder()
	 */
	@Override
	public byte getMaximalMarkovOrder() {
		return (byte)order;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.differentiable.AbstractDifferentiableSequenceScore#getNumberOfRecommendedStarts()
	 */
	@Override
	public int getNumberOfRecommendedStarts() {
		return starts;
	}

	/**
	 * This method allows the user to specify whether the parameters should be
	 * optimized or not.
	 * 
	 * @param optimize
	 *            indicates if the parameters should be optimized or not
	 */
	public void setParameterOptimization(boolean optimize) {
		this.optimize = optimize;
	}

	/**
	 * This method returns the stationary conditional distributions. The first
	 * dimension of the result is used for the order, the second is used for
	 * encoding the context, and the third is used for the different values of
	 * the random variable.
	 * 
	 * <br>
	 * <br>
	 * 
	 * For an homogeneous Markov model of order 2 it returns an array containing
	 * the stationary symbol distribution as first entry, the conditional
	 * stationary distribution of order 1 as second entry and the conditional
	 * distribution of order 2 as third entry.
	 *  
	 * @return all conditional stationary distributions
	 */
	public double[][][] getAllConditionalStationaryDistributions() {
		return StationaryDistribution.getAllConditionalStationaryDistributions(
				probs[order], powers[1]);
	}
	
	/**
	 * Sets the start parameters of this homogeneous Markov model to
	 * the corresponding stationary distributions of the transition probabilities.
	 */
	public void setStartParamsToConditionalStationaryDistributions() {
		double[][][] p = getAllConditionalStationaryDistributions();
		for( int h, o = 0; o < order; o++ ) {
			h = 0;
			for( int i = 0; i < p[o].length; i++ ) {
				for( int j = 0; j < p[o][i].length; j++, h++ ) {
					probs[o][h] = p[o][i][j];
					params[o][h] = Math.log(p[o][i][j]);
				}
				logNorm[o][i] = 0;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM#setStatisticForHyperparameters(int[], double[])
	 */
	@Override
	public void setStatisticForHyperparameters(int[] length, double[] weight)
			throws Exception {
		if (weight.length != length.length) {
			throw new IllegalArgumentException(
					"The length of both arrays (length, weight) have to be identical.");
		}

		double[] sumOfHyperParams = new double[order+1];
		Arrays.fill(sumOfHyperParams, 0);
		for (int l, i = 0; i < length.length; i++) {
			if (weight[i] < 0 || length[i] < 0) {
				throw new IllegalArgumentException(
						"check length and weight for entry " + i);
			} else {
				for (l = 0; l < length[i] && l < order; l++) {
					sumOfHyperParams[l] += weight[i];
				}
				if (order < length[i]) {
					sumOfHyperParams[order] += (length[i] - order) * weight[i];
				}
			}
		}
		hyperParams = getHyperParams( alphabets, sumOfHyperParams );
		computeConstantsOfLogPrior();
	}

	/**
	 * This method returns a {@link DataSet} object containing artificial
	 * sequence(s).
	 * 
	 * <br>
	 * <br>
	 * 
	 * There are 2 different possibilities to create a {@link DataSet}:
	 * <ol>
	 * <li> <code>emitDataSet( int n, int l )</code> returns a {@link DataSet} with
	 * <code>n</code> sequences of length <code>l</code>.
	 * <li> <code>emitDataSet( int n, int[] l )</code> should return a
	 * {@link DataSet} with <code>n</code> sequences which have a sequence length
	 * corresponding to the entry in the array.
	 * </ol>
	 * 
	 * @param numberOfSequences
	 *            the number of sequences that should be contained in the
	 *            returned {@link DataSet}
	 * @param seqLength
	 *            the length of the sequences
	 * 
	 * @return a {@link DataSet} containing <code>numberOfSequences</code>
	 *         artificial sequence(s)
	 * 
	 * @throws Exception
	 *             if the emission of the artificial {@link DataSet} did not
	 *             succeed
	 * 
	 * @see DataSet
	 */
	@Override
	public DataSet emitDataSet(int numberOfSequences, int... seqLength)
			throws Exception {
		Random r = new Random();
		Sequence[] seqs = new Sequence[numberOfSequences];
		byte[] bytes;
		byte a;
		double p;
		for (int l = seqLength[0], parent, o, j, i = 0; i < numberOfSequences; i++) {
			// System.out.println( "sequence " + i );
			if (seqLength.length > 1) {
				l = seqLength[i];
			}
			bytes = new byte[l];
			parent = 0;
			o = 0;
			for (j = 0; j < l; j++) {
				p = r.nextDouble();
				// System.out.print( j + "\t" + p );
				a = 0;
				while (a < powers[1] && probs[o][parent + a] < p) {
					// System.out.print( "\t" + probs[o][parent+a] );
					p -= probs[o][parent + a];
					a++;
				}
				bytes[j] = a;
				// System.out.print( "\t=> " + a );

				parent += a;
				parent *= powers[1];
				parent %= powers[order + 1];
				// System.out.println( "\t-> " + parent );
				if (o < order) {
					o++;
				}
			}
			seqs[i] = new ByteSequence(alphabets, bytes);
			// System.out.println(seqs[i]);
		}
		return new DataSet("generated from " + getInstanceName(), seqs);
	}

	@Override
	public void initializeUniformly(boolean freeParams) {
		double p = 1d / (double) powers[1], v = Math.log(powers[1]);
		for( int o = 0; o <= order; o++) {
			Arrays.fill( params[o], 0 );
			Arrays.fill( logNorm[o], v );
			Arrays.fill( probs[o], p );
		}
		setFreeParams(freeParams);
	}

	@Override
	public int[][] getSamplingGroups(int parameterOffset) {
		int[][] res = new int[getNumberOfParameters() % powers[1]][powers[1]];
		for( int i = 0; i < res.length; i++ ) {
			for( int j = 0; j < res[i].length; j++, parameterOffset++ ) {
				res[i][j] = parameterOffset;
			}
		}
		return res;
	}
}
