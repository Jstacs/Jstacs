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

import java.util.AbstractList;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.ConstraintManager;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.MEMConstraint;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.SequenceIterator;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;

/**
 * This class implements the scoring function for any MRF (Markov Random Field).
 * 
 * @author Jens Keilwagen
 */
public final class MarkovRandomFieldDiffSM extends
		AbstractDifferentiableStatisticalModel {
	private MEMConstraint[] constr;

	private String name;

	private boolean freeParams;

	private int[] offset, help;

	private double ess, norm;
	private double[][] partNorm;
	private SequenceIterator seqIt;

	/**
	 * This constructor creates an instance of a {@link MarkovRandomFieldDiffSM} with
	 * equivalent sample size (ess) 0.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param length
	 *            the length of the sequences and accordingly the model
	 * @param constr
	 *            the constraints that are used for the model, see
	 *            {@link ConstraintManager#extract(int, String)}
	 * 
	 * @see MarkovRandomFieldDiffSM#MarkovRandomFieldDiffSM(AlphabetContainer, int,
	 *      double, String)
	 */
	public MarkovRandomFieldDiffSM(AlphabetContainer alphabets, int length,
			String constr) {
		this(alphabets, length, 0, constr);
	}

	/**
	 * This is the main constructor that creates an instance of a
	 * {@link MarkovRandomFieldDiffSM}.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param length
	 *            the length of the sequences and accordingly the model
	 * @param ess
	 *            the equivalent sample size (ess)
	 * @param constr
	 *            the constraints that are used for the model, see
	 *            {@link ConstraintManager#extract(int, String)}
	 */
	public MarkovRandomFieldDiffSM(AlphabetContainer alphabets, int length,
			double ess, String constr ) {
		super(alphabets, length);
		if (!alphabets.isDiscrete()) {
			throw new IllegalArgumentException(
					"The AlphabetContainer has to be discrete.");
		}
		if (ess < 0) {
			throw new IllegalArgumentException(
					"The ess has to be non-negative.");
		}
		this.ess = ess;
		int[] aLength = new int[length];
		for (int i = 0; i < length; aLength[i] = (int) alphabets
				.getAlphabetLengthAt(i++))
			;
		AbstractList<int[]> list = ConstraintManager.extract(length, constr);
		ConstraintManager.reduce(list);
		this.constr = ConstraintManager.createConstraints(list, aLength);
		this.name = constr;
		freeParams = false;
		getNumberOfParameters();
		init(Double.NaN);
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link MarkovRandomFieldDiffSM} out of a {@link StringBuffer} as
	 * returned by {@link #toXML()}.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public MarkovRandomFieldDiffSM(StringBuffer source) throws NonParsableException {
		super(source);
	}

	private void init(double n) {
		norm = n;
		if (partNorm == null) {
			partNorm = new double[constr.length][];
			for (int i = 0; i < partNorm.length; i++) {
				partNorm[i] = new double[constr[i]
						.getNumberOfSpecificConstraints()];
			}
			help = new int[2];
			int[] aLength = new int[length];
			for (int i = 0; i < length; i++ ) {
				aLength[i] = (int) alphabets.getAlphabetLengthAt(i);
			}
			seqIt = new SequenceIterator(length);
			seqIt.setBounds(aLength);
		} else {
			for (int i = 0; i < partNorm.length; i++) {
				Arrays.fill(partNorm[i], n);
			}
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
	protected void fromXML(StringBuffer representation)
			throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag(representation, XML_TAG);
		length = XMLParser.extractObjectForTags(xml, "length", int.class );
		alphabets = XMLParser.extractObjectForTags(xml, "alphabets", AlphabetContainer.class );
		ess = XMLParser.extractObjectForTags(xml, "ess", double.class );
		name = XMLParser.extractObjectForTags(xml, "name", String.class );
		constr = XMLParser.extractObjectForTags(xml, "constr", MEMConstraint[].class );
		freeParams = XMLParser.extractObjectForTags(xml, "freeParams", boolean.class );
		getNumberOfParameters();
		init(Double.NaN);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel#clone()
	 */
	@Override
	public MarkovRandomFieldDiffSM clone() throws CloneNotSupportedException {
		MarkovRandomFieldDiffSM clone = (MarkovRandomFieldDiffSM) super.clone();
		clone.constr = ArrayHandler.clone(constr);
		clone.partNorm = new double[partNorm.length][];
		for (int i = 0; i < partNorm.length; i++) {
			clone.partNorm[i] = partNorm[i].clone();
		}
		clone.norm = norm;

		clone.help = help.clone();
		clone.offset = null;
		clone.getNumberOfParameters();
		clone.seqIt = seqIt.clone();
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScore(de.jstacs.data
	 * .Sequence, int)
	 */
	public double getLogScoreFor(Sequence seq, int start) {
		double erg = 0;
		for (int i = 0; i < constr.length; i++) {
			erg += constr[i].getLambda( constr[i].satisfiesSpecificConstraint( seq, start) );
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScoreAndPartialDerivation
	 * (de.jstacs.data.Sequence, int, de.jstacs.utils.IntList,
	 * de.jstacs.utils.DoubleList)
	 */
	public double getLogScoreAndPartialDerivation(Sequence seq, int start,
			IntList indices, DoubleList partialDer) {
		double erg = 0;
		int i = 0, j, z;
		for (; i < constr.length; i++) {
			j = constr[i].satisfiesSpecificConstraint(seq, start);
			if ((z = offset[i] + j) < offset[i + 1]) {
				indices.add(z);
				partialDer.add(1);
			}
			erg += constr[i].getLambda(j);
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getNumberOfParameters()
	 */
	public int getNumberOfParameters() {
		if (offset == null) {
			int i = 0, anz = 0;
			offset = new int[constr.length + 1];
			while (i < constr.length) {
				anz += constr[i++].getNumberOfSpecificConstraints();
				if (freeParams) {
					anz -= 1;
				}
				offset[i] = anz;
			}
		}

		return offset[constr.length];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getInstanceName()
	 */
	public String getInstanceName() {
		return "MRF(" + name + ")";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#setParameters(double[],
	 * int)
	 */
	public void setParameters(double[] params, int start) {
		norm = Double.NaN;
		int i = 0, j, s = offset[0];
		for (; i < constr.length; i++) {
			//TODO constr[i].setLambda(constr[i].getNumberOfSpecificConstraints() - 1, 0);
			for (j = 0; s < offset[i + 1]; s++, j++) {
				constr[i].setLambda(j, params[start + s]);
			}
		}
	}

	
	public String toString() {
		int i = 0, j;
		StringBuffer res = new StringBuffer();
		res.append( getInstanceName() + "\n" );
		for( ; i < constr.length; i++ ) {
			res.append( constr[i] );
			for( j = 0; j < constr[i].getNumberOfSpecificConstraints(); j++ ) {
				res.append( "\t" + constr[i].getLambda( j ) );
			}
			res.append( "\n" );
		}
		return res.toString();
	}/**/

	private static final String XML_TAG = "MarkovRandomFieldDiffSM";

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer b = new StringBuffer(10000);
		XMLParser.appendObjectWithTags(b, length, "length");
		XMLParser.appendObjectWithTags(b, alphabets, "alphabets");
		XMLParser.appendObjectWithTags(b, ess, "ess");
		XMLParser.appendObjectWithTags(b, name, "name");
		XMLParser.appendObjectWithTags(b, constr, "constr");
		XMLParser.appendObjectWithTags(b, freeParams, "freeParams");
		XMLParser.addTags(b, XML_TAG);
		return b;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#initializeFunction(int,
	 * boolean, de.jstacs.data.DataSet[], double[][])
	 */
	public void initializeFunction(int index, boolean freeParams,
			DataSet[] data, double[][] weights) throws Exception {
		if (this.freeParams != freeParams) {
			offset = null;
			this.freeParams = freeParams;
			getNumberOfParameters();
		}

		// uniform
		double d = 0;
		for (int i = 0; i < length; i++) {
			d -= Math.log(alphabets.getAlphabetLengthAt(i));
		}
		d /= constr.length;
		for (int k, j, i = 0; i < constr.length; i++) {
			k = constr[i].getNumberOfSpecificConstraints();
			for (j = 0; j < k; j++) {
				constr[i].setLambda(j, d);
			}
		}
		
		/*TODO
		//somehow data specific
		ConstraintManager.countInhomogeneous( alphabets, length, data[index], weights == null ? null : weights[index], true, constr );
		for( int k, j, i = 0; i < constr.length; i++ ) {
			constr[i].estimate( ess );
			k = constr[i].getNumberOfSpecificConstraints();
			for (j = 0; j < k; j++) {
				constr[i].setExpLambda( j, constr[i].getFreq( j ) );
			}
		}/**/

		norm = Double.NaN;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#initializeFunctionRandomly
	 * (boolean)
	 */
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		if (this.freeParams != freeParams) {
			offset = null;
			this.freeParams = freeParams;
			getNumberOfParameters();
		}
		for (int k, j, i = 0; i < constr.length; i++) {
			k = constr[i].getNumberOfSpecificConstraints();
			for (j = 0; j < k; j++) {
				constr[i].setLambda(j, r.nextGaussian()/(double)k );
			}
		}
		norm = Double.NaN;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#
	 * getLogNormalizationConstant()
	 */
	public double getLogNormalizationConstant() {
		if ( Double.isNaN( norm ) ) {
			precompute();
		}
		return norm;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#
	 * getLogPartialNormalizationConstant(int)
	 */
	public double getLogPartialNormalizationConstant(int parameterIndex)
			throws Exception {
		if ( Double.isNaN( norm ) ) {
			precompute();
		}
		computeIndices( parameterIndex );
		return partNorm[help[0]][help[1]];
	}

	private void precompute() {
		// TODO current implementation is only the naive approach, so try to
		// make this faster?!?
		seqIt.reset();
		int i;
		int[] fulfilled = new int[constr.length];
		seqIt.reset();
		double s;
		init(Double.NEGATIVE_INFINITY);
		do {
			s = getLogScore(fulfilled, seqIt);
			for (i = 0; i < constr.length; i++) {
				partNorm[i][fulfilled[i]] = Normalisation.getLogSum( partNorm[i][fulfilled[i]], s );
			}
			norm = Normalisation.getLogSum( norm, s );
		} while (seqIt.next());
	}

	private double getLogScore(int[] fulfilled, SequenceIterator sequence) {
		double s = 0;
		for (int counter = 0; counter < constr.length; counter++) {
			fulfilled[counter] = constr[counter].satisfiesSpecificConstraint(sequence);
			s += constr[counter].getLambda(fulfilled[counter]);
		}
		return s;
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
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#
	 * getSizeOfEventSpaceForRandomVariablesOfParameter(int)
	 */
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		computeIndices(index);
		return constr[help[0]].getNumberOfSpecificConstraints();
	}

	private void computeIndices(int index) {
		help[0] = 0;
		while (index >= offset[help[0]]) {
			help[0]++;
		}
		help[0]--;
		help[1] = index - offset[help[0]];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorTerm() {
		double logPriorTerm = 0, d;
		int i = 0, j, s = 0;
		for (; i < constr.length; i++) {
			s = constr[i].getNumberOfSpecificConstraints();
			d = ess / (double) s;
			for (j = 0; j < s; j++) {
				logPriorTerm += constr[i].getLambda(j) * d;
			}
		}
		return logPriorTerm;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#
	 * addGradientOfLogPriorTerm(double[], int)
	 */
	public void addGradientOfLogPriorTerm(double[] grad, int start) {
		double d;
		int i = 0, s = offset[0];
		while( i < constr.length ) {
			d = ess / (double) constr[i].getNumberOfSpecificConstraints();
			i++;
			for(; s < offset[i]; s++, start++) {
				grad[start] += d;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getCurrentParameterValues()
	 */
	public double[] getCurrentParameterValues() throws Exception {
		double[] start = new double[offset[constr.length]];
		for( int j, i = 0, n = 0; i < constr.length; i++ ) {
			for( j = 0; n < offset[i + 1]; n++, j++ ) {
				start[n] = constr[i].getLambda(j);
			}
		}
		return start;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#isInitialized()
	 */
	public boolean isInitialized() {
		return true;
	}
}
