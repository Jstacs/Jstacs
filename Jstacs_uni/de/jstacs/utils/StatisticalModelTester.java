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

package de.jstacs.utils;

import java.io.IOException;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.DataSet.ElementEnumerator;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.sequenceScores.SequenceScore;
import de.jstacs.sequenceScores.statisticalModels.StatisticalModel;

/**
 * This class is useful for some test for any (discrete) models. It implements
 * several statistics (log-likelihood, Shannon entropy, AIC, BIC, ...) to
 * compare models.
 * 
 * @see de.jstacs.sequenceScores.statisticalModels.StatisticalModel
 * 
 * @author Jens Keilwagen
 */
public class StatisticalModelTester {
	/**
	 * Returns the Kullback-Leibler-divergence <code>D(p_m1||p_m2)</code>.
	 * 
	 * <br>
	 * <br>
	 * 
	 * Computes <code>\sum_x p(x|m1) * \log \frac{p(x|m1)}{p(x|m2)}</code>.
	 * 
	 * @param m1
	 *            one discrete model
	 * @param m2
	 *            another discrete model
	 * @param length
	 *            the length of the sequence (for inhomogeneous models length
	 *            has to be {@link StatisticalModel#getLength()})
	 * 
	 * @return the Kullback-Leibler-divergence
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getKLDivergence(StatisticalModel m1, StatisticalModel m2, int length)
			throws Exception {
		SeqIterator s = new SeqIterator(m1.getAlphabetContainer(), length);
		Sequence seq;
		double kl = 0, v;
		do {
			seq = s.getSequence();
			v = m1.getLogProbFor(seq);
			kl += Math.exp(v) * (v - m2.getLogProbFor(seq));
		} while (s.next());
		return kl;
	}

	/**
	 * Returns the difference of the Kullback-Leibler-divergences, i.e.
	 * <code>D(p_m1||p_m2) - D(p_m2||p_m1)</code>.
	 * 
	 * <br>
	 * <br>
	 * 
	 * Computes
	 * <code>\sum_x (p(x|m1)-p(x|m2)) * \log \frac{p(x|m1)}{p(x|m2)}</code>.
	 * 
	 * @param m1
	 *            one discrete model
	 * @param m2
	 *            another discrete model
	 * @param length
	 *            the length of the sequence (for inhomogeneous models length
	 *            has to be {@link StatisticalModel#getLength()})
	 * 
	 * @return the difference of the Kullback-Leibler-divergence
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getSymKLDivergence(StatisticalModel m1, StatisticalModel m2, int length)
			throws Exception {
		SeqIterator s = new SeqIterator(m1.getAlphabetContainer(), length);
		Sequence seq;
		double kl = 0, logP1, logP2;
		do {
			seq = s.getSequence();
			logP1 = m1.getLogProbFor(seq);
			logP2 = m2.getLogProbFor(seq);
			kl += (Math.exp(logP1) - Math.exp(logP2)) * (logP1 - logP2);
		} while (s.next());
		return kl;
	}

	/**
	 * Returns the log-likelihood of a {@link DataSet} <code>data</code> for a
	 * given model <code>m</code>.
	 * 
	 * @param m
	 *            the given model
	 * @param data
	 *            the {@link DataSet}
	 * 
	 * @return the log-likelihood of <code>data</code>
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getLogLikelihood(StatisticalModel m, DataSet data)
			throws Exception {
		return getLogLikelihood(m, data, null);
	}

	/**
	 * Returns the log-likelihood of a {@link DataSet} <code>data</code> for a
	 * given model <code>m</code>.
	 * 
	 * @param m
	 *            the given model
	 * @param data
	 *            the {@link DataSet}
	 * @param weights
	 *            the weight for each element of the {@link DataSet}
	 * 
	 * @return the log-likelihood of <code>data</code>
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getLogLikelihood(StatisticalModel m, DataSet data, double[] weights)
			throws Exception {
		int counter, d = data.getNumberOfElements();
		double erg = 0;
		ElementEnumerator ei = new ElementEnumerator(data);
		if (weights == null) {
			for (counter = 0; counter < d; counter++) {
				erg += m.getLogProbFor(ei.nextElement());
			}
		} else if (d != weights.length) {
			throw new IllegalArgumentException(
					"The weights and the data set does not match.");
		} else {
			for (counter = 0; counter < d; counter++) {
				erg += weights[counter] * m.getLogProbFor(ei.nextElement());
			}
		}
		return erg;
	}

	/**
	 * This method computes the marginal distribution for any discrete model
	 * <code>m</code> and all sequences that fulfill the <code>constraint</code>
	 * , if possible.
	 * 
	 * @param m
	 *            a discrete model
	 * @param constraint
	 *            <code>constraint[j][i] < 0</code> stands for an irrelevant
	 *            position, <code>constraint[j][i] = c</code> with
	 *            <code>0 <= c < m.getAlphabets()[(m.getLength==0)?0:i].getAlphabetLength()</code>
	 *            is the encoded character of position <code>i</code>
	 * 
	 * @return the marginal distributions of the given constraints for the discrete model
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double[] getMarginalDistribution(StatisticalModel m, int[]... constraint)
			throws Exception {
		int l = constraint[0].length, len = m.getLength();
		for( int i = 0; i < constraint.length; i++ ) {
			if (l!= constraint[i].length || (len != 0 && len != constraint[i].length)) {
				throw new IOException(
						"This model can only classify sequences of length "
								+ m.getLength() + ".");
			}
		}
		double[] erg = new double[constraint.length];
		Arrays.fill(erg, Double.NEGATIVE_INFINITY);
		SeqIterator s = new SeqIterator(m.getAlphabetContainer(),l);
		do {
			for( int i = 0; i < constraint.length; i++ ) {
				if (s.isSatisfied(constraint[i])) {
					erg[i] = Normalisation.getLogSum( erg[i], m.getLogProbFor(s.getSequence()) );
				}
			}
		} while (s.next());
		for( int i = 0; i < constraint.length; i++ ) {
			erg[i] = Math.exp( erg[i] );
		}
		return erg;
	}

	/**
	 * This method computes the maximum deviation between the probabilities for
	 * all sequences of <code>length</code> for discrete models <code>m1</code>
	 * and <code>m2</code>.
	 * 
	 * @param m1
	 *            one discrete model
	 * @param m2
	 *            another discrete model
	 * @param length
	 *            the length of the sequence (for inhomogeneous models length
	 *            has to be {@link StatisticalModel#getLength()})
	 * 
	 * @return the maximum deviation between the probabilities
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getMaxOfDeviation(StatisticalModel m1, StatisticalModel m2, int length)
			throws Exception {
		if (m1.getLength() != 0 && m1.getLength() != length) {
			throw new IOException(
					"The model m1 can only classify sequences of length "
							+ m1.getLength() + ".");
		}
		if (m2.getLength() != 0 && m2.getLength() != length) {
			throw new IOException(
					"This model m2 can only classify sequences of length "
							+ m2.getLength() + ".");
		}
		if (!m1.getAlphabetContainer().checkConsistency(
				m2.getAlphabetContainer())) {
			throw new IOException(
					"The models are training on different alphabets.");
		}
		double p, max = 0;
		SeqIterator s = new SeqIterator(m1.getAlphabetContainer(), length);
		Sequence seq;
		do {
			seq = s.getSequence();
			p = Math.abs(Math.exp(m1.getLogProbFor(seq, 0, s.last)) - Math.exp(m2.getLogProbFor(seq, 0, s.last)));
			if (p > max) {
				max = p;
			}
		} while (s.next());
		return max;
	}

	/**
	 * Returns one most probable sequence for the discrete model <code>m</code>.
	 * (Maybe there are more than one most probable sequences. In this case only
	 * one of them is returned.)
	 * 
	 * <br>
	 * <br>
	 * 
	 * This is only a standard implementation. For some special models like
	 * Markov models it is possible to compute the probabilities of the
	 * sequences much faster by using a dynamic-programming-algorithm.
	 * 
	 * @param m
	 *            the discrete model
	 * @param length
	 *            the length of the sequence (for inhomogeneous models length
	 *            has to be {@link SequenceScore#getLength()})
	 * 
	 * @return one most probable sequence
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static Sequence getMostProbableSequence(SequenceScore m, int length)
			throws Exception {
		SeqIterator s = new SeqIterator(m.getAlphabetContainer(), length);
		Sequence current, seq = s.getSequence();
		double pmax = m.getLogScoreFor(seq), p;
		while (s.next()) {
			current = s.getSequence();
			p = m.getLogScoreFor(current);
			if (p > pmax) {
				pmax = p;
				seq = current;
			}
		}
		return seq;
	}

	/**
	 * This method computes the Shannon entropy for any discrete model
	 * <code>m</code> and all sequences of <code>length</code>, if possible.
	 * 
	 * @param m
	 *            the discrete model
	 * @param length
	 *            the length of the sequence (for inhomogeneous models length
	 *            has to be {@link StatisticalModel#getLength()})
	 * 
	 * @return the Shannon entropy for a discrete model
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getShannonEntropy(StatisticalModel m, int length)
			throws Exception {
		if (m.getLength() != 0 && m.getLength() != length) {
			throw new IOException(
					"This model can only classify sequences of length "
							+ m.getLength() + ".");
		}
		double logP, erg = 0;
		SeqIterator s = new SeqIterator(m.getAlphabetContainer(), length);
		do {
			logP = m.getLogProbFor(s.getSequence());
			if ( !Double.isInfinite( logP ) ) {
				erg -= Math.exp(logP) * logP;
			}
			if ( logP > 0 ) {
				throw new IOException("The probability of sequence "
						+ s.getSequence() + " is not correct (" + Math.exp(logP) + ").");
			}
		} while (s.next());
		return erg;
	}

	/**
	 * This method computes the Shannon entropy in bits for any discrete model
	 * <code>m</code> and all sequences of <code>length</code>, if possible.
	 * 
	 * @param m
	 *            the discrete model
	 * @param length
	 *            the length of the sequence (for inhomogeneous models length
	 *            has to be {@link StatisticalModel#getLength()})
	 * 
	 * @return the Shannon entropy in bits for a discrete model
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getShannonEntropyInBits(StatisticalModel m, int length)
			throws Exception {
		return (getShannonEntropy(m, length) / Math.log(2));
	}

	/**
	 * This method computes the sum of deviations between the probabilities for
	 * all sequences of <code>length</code> for discrete models <code>m1</code>
	 * and <code>m2</code>.
	 * 
	 * @param m1
	 *            one discrete model
	 * @param m2
	 *            another discrete model
	 * @param length
	 *            the length of the sequence (for inhomogeneous models length
	 *            has to be {@link StatisticalModel#getLength()})
	 * 
	 * @return the sum of deviations between the probabilities
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getSumOfDeviation(StatisticalModel m1, StatisticalModel m2, int length)
			throws Exception {
		if (m1.getLength() != 0 && m1.getLength() != length) {
			throw new IOException(
					"The model m1 can only classify sequences of length "
							+ m1.getLength() + ".");
		}
		if (m2.getLength() != 0 && m2.getLength() != length) {
			throw new IOException(
					"This model m2 can only classify sequences of length "
							+ m2.getLength() + ".");
		}
		if (!m1.getAlphabetContainer().checkConsistency(
				m2.getAlphabetContainer())) {
			throw new IOException(
					"The models are training on different alphabets.");
		}
		double sum = 0;
		SeqIterator s = new SeqIterator(m1.getAlphabetContainer(), length);
		Sequence seq;
		do {
			seq = s.getSequence();
			sum += Math.abs(Math.exp(m1.getLogProbFor(seq, 0, s.last)) - Math.exp(m2.getLogProbFor(seq, 0, s.last)));
		} while (s.next());
		return sum;
	}

	/**
	 * This method computes the marginal distribution for any discrete model
	 * <code>m</code> and all sequences of <code>length</code>, if possible. So
	 * this method can be used to give a hint whether a model is a distribution
	 * or if some mistakes are in the implementation.
	 * 
	 * <br>
	 * <br>
	 * 
	 * It is expected that this method delivers the value 1.0, but because of
	 * the limited precision in Java the value 1.0 is unrealistic.
	 * 
	 * <br>
	 * <br>
	 * 
	 * <code>Math.abs( 1.0d - getSumOfDistribution( m, length )</code> should be
	 * smaller than <code>1E-10</code>.
	 * 
	 * @param m
	 *            the discrete model
	 * @param length
	 *            the length of the sequence (for inhomogeneous models length
	 *            has to be {@link StatisticalModel#getLength()})
	 * 
	 * @return the marginal distribution for a discrete model
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getSumOfDistribution(StatisticalModel m, int length)
			throws Exception {
		if (m.getLength() != 0 && m.getLength() != length) {
			throw new IOException(
					"This model can only classify sequences of length "
							+ m.getLength() + ".");
		}
		double erg = Double.NEGATIVE_INFINITY, p;
		SeqIterator s = new SeqIterator(m.getAlphabetContainer(), length);
		do {
			p = m.getLogProbFor(s.getSequence());
			if (p > 0) {
				throw new IOException("The probability (" + Math.exp(p)
						+ ") for sequence \"" + s.getSequence()
						+ "\" is not in [0,1].");
			}
			erg = Normalisation.getLogSum( erg, p );
		} while (s.next());
		return Math.exp(erg);
	}

	/**
	 * This method computes the value of Akaikes Information Criterion (AIC). It
	 * uses the formula: AIC = <code>2 * log L(t,x) - 2*k</code>, where
	 * <code>L(t,x)</code> is the likelihood of the {@link DataSet} and
	 * <code>k</code> is the number of parameters in the model.
	 * 
	 * <br>
	 * <br>
	 * 
	 * The value of the AIC can be used for model selection.
	 * 
	 * @param m
	 *            a trained model
	 * @param s
	 *            the {@link DataSet} for the test
	 * @param k
	 *            the number of parameters of the model <code>m</code>
	 * 
	 * @return the value of AIC
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getValueOfAIC(StatisticalModel m, DataSet s, int k)
			throws Exception {
		return 2 * getLogLikelihood(m, s) - 2 * k;
	}

	/**
	 * This method computes the value of the Bayesian Information Criterion
	 * (BIC). It uses the formula: BIC = <code>2 * log L(t,x) - k *
	 * log n</code>, where <code>L(t,x)</code> is the likelihood of the
	 * {@link DataSet}, <code>k</code> is the number of parameters in the model
	 * and <code>n</code> is the number of sequences in the {@link DataSet}.
	 * 
	 * <br>
	 * <br>
	 * 
	 * The value of the BIC can be used for model selection.
	 * 
	 * @param m
	 *            a trained model
	 * @param s
	 *            the {@link DataSet} for the test
	 * @param k
	 *            the number of parameters of the model <code>m</code>
	 * 
	 * @return value of AIC
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static double getValueOfBIC(StatisticalModel m, DataSet s, int k)
			throws Exception {
		return 2 * getLogLikelihood(m, s) - k * StrictMath.log(s.getNumberOfElements());
	}

	private static class SeqIterator {
		// the sequence
		private int[] seq;

		// simple or not
		private boolean simple;

		// the alphabet.length
		private int[] a;

		// length of the seq
		private int l, last;

		// the alphabets
		private AlphabetContainer abc;

		private SeqIterator(AlphabetContainer abc, int length)
				throws IllegalArgumentException {
			if (!abc.isDiscrete()) {
				throw new IllegalArgumentException("The model is not discrete.");
			}
			this.abc = abc;
			simple = abc.isSimple();
			a = new int[(simple ? 1 : length) + 1];
			int i = 0;
			for (; i < a.length - 1; i++) {
				a[i] = (int) abc.getAlphabetLengthAt(i) - 1;
			}
			a[i] = 1;
			l = length;
			last = l - 1;
			seq = new int[length + 1];
		}

		private boolean next() {
			int s_index = 0;
			while (seq[s_index] == a[simple ? 0 : s_index]) {
				seq[s_index++] = 0;
			}
			seq[s_index]++;
			return seq[l] == 0;
		}

		private boolean isSatisfied(int[] constr) {
			int i = 0;
			while (i < constr.length
					&& (constr[i] == -1 || constr[i] == seq[i])) {
				i++;
			}
			return (i == constr.length);
		}

		private Sequence getSequence() throws WrongAlphabetException,
				WrongSequenceTypeException {
			return new IntSequence(abc, seq, 0, l);
		}
	}
}