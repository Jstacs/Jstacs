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

import java.util.Random;

import de.jstacs.NotTrainedException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.StatisticalModel;

/**
 * Emits {@link DataSet}s for discrete inhomogeneous models by a naive implementation.
 * 
 * @see StatisticalModel#emitDataSet(int, int[])
 * 
 * @author Jens Keilwagen
 */
public class DiscreteInhomogenousDataSetEmitter {
	/**
	 * This method emits a data set with
	 * <code>n<code> sequences from the discrete inhomogeneous model <code>m</code>
	 * .
	 * 
	 * @param m
	 *            the model
	 * @param n
	 *            the number of sequences
	 * 
	 * @return the synthetic data set
	 * 
	 * @throws NotTrainedException
	 *             if the model is not trained
	 * @throws Exception
	 *             if something went wrong
	 */
	public static DataSet emitDataSet(StatisticalModel m, int n) throws NotTrainedException,
			Exception {
		if (!m.isInitialized()) {
			throw new NotTrainedException();
		}
		AlphabetContainer alphabet = m.getAlphabetContainer();
		if (!alphabet.isDiscrete()) {
			throw new IllegalArgumentException("The models has to be discrete.");
		}
		int length = m.getLength(), counter1, counter2, counter3, il1 = 1, il2 = 1, l = 0;
		if (length == 0) {
			throw new IllegalArgumentException(
					"The models has to be inhomogenous.");
		}

		int[] alphabetLength = new int[length + 1];
		boolean firstPart = true;
		for (counter1 = 0; counter1 < length; counter1++) {
			alphabetLength[counter1] = (int) alphabet
					.getAlphabetLengthAt(counter1);
			if (firstPart && il1 < Integer.MAX_VALUE / alphabetLength[counter1]) {
				il1 *= alphabetLength[counter1];
				l++;
			} else {
				firstPart = false;
				if (il2 < Integer.MAX_VALUE / alphabetLength[counter1]) {
					il2 *= alphabetLength[counter1];
				} else {
					throw new IllegalArgumentException(
							"It is not possible to emit a data set of sequences with this length and alphabets by this implementation. (needs to much memory)");
				}
			}
		}
		alphabetLength[length] = 2;

		double p = 0;
		// fill the array for all discrete values
		double[][] cumulative_p = new double[il2][il1];
		int[] sequence = new int[length + 1];
		for (counter2 = 0; counter2 < il2; counter2++) {
			for (counter1 = 0; counter1 < il1; counter1++) {
				p += Math.exp( m.getLogProbFor(new IntSequence(alphabet, sequence, 0, length)) );
				cumulative_p[counter2][counter1] = p;
				counter3 = 0;
				while (sequence[counter3] == alphabetLength[counter3] - 1) {
					sequence[counter3++] = 0;
				}
				sequence[counter3]++;
			}
		}

		// draw sequences (using binary search)
		Random r = new Random();
		int max, min;
		Sequence[] seqs = new Sequence[n];
		sequence = new int[length];
		for (int i = 0; i < n; i++) {
			// first (rough) search
			max = il2 - 1;
			min = -1;
			p = r.nextDouble();
			while (max - min > 1) {
				counter1 = (max + min) / 2;
				if (p >= cumulative_p[counter1][il1 - 1]) {
					min = counter1;
				} else {
					max = counter1;
				}
			}
			// cumulative_p[min][il-1] <= p < cumulative_p[max][il-1];
			counter1 = max;

			// second (fine) search
			// "cumulative_p[counter1][-1]" <= p < cumulative_p[max][il-1];
			max = il1 - 1;
			min = -1;
			while (max - min > 1) {
				counter2 = (max + min) / 2;
				if (p >= cumulative_p[counter1][counter2]) {
					min = counter2;
				} else {
					max = counter2;
				}
			}
			// cumulative_p[counter1][min] <= p < cumulative_p[counter1][max];
			counter2 = max;

			// decode
			counter3 = counter2;
			for (min = 0; min < l; min++) {
				sequence[min] = counter3 % alphabetLength[min];
				counter3 /= alphabetLength[min];
			}
			counter3 = counter1;
			for (; min < length; min++) {
				sequence[min] = counter3 % alphabetLength[min];
				counter3 /= alphabetLength[min];
			}
			seqs[i] = new IntSequence(alphabet, sequence);
		}

		return new DataSet("sampled from " + m.getInstanceName(), seqs);
	}
}
