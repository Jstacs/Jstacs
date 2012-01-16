/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous;

import java.util.Arrays;

/**
 * This class is used to iterate over a discrete sequence.
 * 
 * @author Jens Keilwagen
 */
public class SequenceIterator implements Cloneable {

	// the sequence
	int[] seq;

	// the length of the alphabet (maximal number of symbols) at each position
	private int[] maxSymbol;

	private int l;

	/**
	 * Creates a new {@link SequenceIterator} with maximal <code>length</code>.
	 * 
	 * @param length
	 *            the maximal length of the sequence
	 */
	public SequenceIterator( int length ) {
		seq = new int[length + 1];
		maxSymbol = null;
		l = -1;
	}
	
	public SequenceIterator clone() throws CloneNotSupportedException {
		SequenceIterator clone = (SequenceIterator) super.clone();
		clone.seq = seq.clone();
		clone.maxSymbol = maxSymbol.clone();
		return clone;
	}

	///**
	// * This method returns the value of the internal sequence at a given index.
	// * 
	// * @param index
	// *            the index
	// * 
	// * @return the value of the internal sequence at position index
	// */
	//public int getValueAt( int index )
	//{
	//	return seq[index];
	//}

	/**
	 * Changes the internal sequence representation to the next sequence.
	 * 
	 * @return <code>true</code> if the new sequence is correct, otherwise
	 *         <code>false</code>
	 */
	public boolean next() {
		int s_index = 0;
		while( seq[s_index] == maxSymbol[s_index] ) {
			seq[s_index++] = 0;
		}
		seq[s_index]++;
		return seq[l] == 0;
	}

	/**
	 * Resets the internal sequence representation. So the
	 * {@link SequenceIterator} starts again.
	 */
	public void reset() {
		Arrays.fill( seq, 0 );
	}

	/**
	 * This method sets the bounds for each position. It does not copy the array
	 * and it does not proof <code>bounds.length
	 * &lt;= length</code>. This has to be ensured by the user.
	 * 
	 * @param bounds
	 *            the array with bounds for each position
	 */
	public void setBounds( int[] bounds ) {
		maxSymbol = new int[bounds.length + 1];
		for( int counter = 0; counter < bounds.length; counter++ ) {
			maxSymbol[counter] = bounds[counter] - 1;
		}
		l = bounds.length;
		maxSymbol[l] = 1;
		reset();
	}

	/**
	 * This method skips some position.
	 * 
	 * @param firstPos
	 *            the first position that is interesting
	 * 
	 * @return <code>true</code> if the internal sequence is from the scope,
	 *         <code>false</code> otherwise
	 */
	public boolean skip( int firstPos ) {
		int counter1 = 0;
		while( counter1 < firstPos ) {
			seq[counter1++] = 0;
		}
		while( seq[counter1] == maxSymbol[counter1] ) {
			seq[counter1++] = 0;
		}
		seq[counter1]++;
		return seq[l] == 0;
	}
	
	/**
	 * This method returns the discrete value for a specific position.
	 * 
	 * @param pos the position
	 * 
	 * @return the value for the given position
	 */
	public int discreteValAt( int pos ) {
		if( pos < l ) {
			return seq[pos];
		} else {
			throw new IndexOutOfBoundsException();
		}
	}
	
	/**
	 * This method returns the number of sequences in this iterator,
	 * i.e., the number of times {@link #next()} returns <code>true</code> after using {@link #reset()}.
	 * 
	 * @return the number of sequences in this iterator
	 */
	public int getNumberOfSequences() {
		int res = maxSymbol[0]+1;
		for( int i = 1; i < maxSymbol.length-1; i++ ) {
			res *= (maxSymbol[i]+1);
		}
		return res;
	}
}
