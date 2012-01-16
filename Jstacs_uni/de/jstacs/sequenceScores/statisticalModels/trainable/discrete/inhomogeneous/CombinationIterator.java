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

/**
 * This class can be used for iterating over all possible combinations (in the
 * sense of combinatorics).
 * 
 * @author Jens Keilwagen
 */
public class CombinationIterator {

	private long[][] pascal;

	private int max, anz;

	private int[] currentCombi;

	/**
	 * Creates a new {@link CombinationIterator} with <code>n</code> elements
	 * and at most <code>max</code> selected elements.
	 * 
	 * @param n
	 *            the number of elements
	 * @param max
	 *            the maximal number of selected elements
	 */
	public CombinationIterator( int n, int max ) {
		int j = n - 1, i = 0, l = max + 1;
		this.max = max;
		this.anz = n;
		pascal = new long[n][l];
		for( i = 0; i < n; i++ ) {
			pascal[i][0] = 1;
		}

		for( i = n - 1; i >= 0; i-- ) {
			l = Math.min( max, n - i - 1 );
			for( j = 1; j <= l; j++ ) {
				pascal[i][j] = pascal[i + 1][j - 1] + pascal[i + 1][j];
			}
		}

		/*
		 * for( i = 0; i < n; i++ ) System.out.println( i + "\t" + Arrays.toString(pascal[i]) ); System.out.println();
		 */
	}

	/**
	 * Returns the number of possible combinations.
	 * 
	 * @param elements
	 *            the number of selected elements
	 * 
	 * @return the number of possible combinations
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>elements &lt; 0</code>
	 */
	public long getNumberOfCombinations( int elements ) throws IllegalArgumentException {
		if( elements <= 0 ) {
			if( elements == 0 ) {
				return 0;
			} else {
				throw new IllegalArgumentException( "out of range." );
			}
		} else {
			return pascal[0][elements] + pascal[0][elements - 1];
		}
	}

	/**
	 * This method sets the current used number of selected elements.
	 * 
	 * @param current
	 *            the current used number of selected elements
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>current &lt; 0</code> or
	 *             <code>current &gt; max</code>
	 * 
	 * @see CombinationIterator#reset()
	 */
	public void setCurrentLength( int current ) throws IllegalArgumentException {
		if( current < 0 || current > max ) {
			throw new IllegalArgumentException( "The value has to be in [0," + max + "]." );
		}
		currentCombi = new int[current];
		reset();
	}

	/**
	 * Steps to the next combination.
	 * 
	 * @return <code>true</code> if the method does not run out of range,
	 *         <code>false</code> otherwise
	 * 
	 * @throws IllegalArgumentException
	 *             if the method {@link #setCurrentLength(int)}was not used to
	 *             set the current used number of selected elements
	 */
	public boolean next() throws IllegalArgumentException {
		if( currentCombi == null ) {
			throw new IllegalArgumentException( "There is no length set (use method setCurrentLength(int))." );
		}
		int counter1 = currentCombi.length - 1, counter2 = 1;
		while( counter1 >= 0 && currentCombi[counter1] == anz - counter2 ) {
			counter1--;
			counter2++;
		}
		if( counter1 < 0 ) {
			return false;
		} else {
			currentCombi[counter1++]++;
			while( counter1 < currentCombi.length ) {
				currentCombi[counter1] = currentCombi[counter1 - 1];
				currentCombi[counter1++]++;
			}
			return true;
		}
	}

	/**
	 * Returns a clone of the internal combination.
	 * 
	 * @return a clone of the internal combination
	 */
	public int[] getCombination() {
		return currentCombi.clone();
	}

	/**
	 * This method returns an index for the sorted entries of a combination
	 * <code>combi</code>.
	 * 
	 * @param combi
	 *            the combination
	 * 
	 * @return the index
	 */
	public long getIndex( int[] combi ) {
		long erg = pascal[combi[0]][combi.length];
		for( int i = 1; i < combi.length; i++ ) {
			erg += pascal[combi[i]][combi.length - i];
		}
		return erg;
	}

	/**
	 * This method resets the internal combination.
	 * 
	 * @throws IllegalArgumentException
	 *             if the method {@link #setCurrentLength(int)} was not used to
	 *             set the current used number of selected elements
	 */
	public void reset() throws IllegalArgumentException {
		if( currentCombi == null ) {
			throw new IllegalArgumentException( "There is no length set (use method setCurrentLength(int))." );
		}
		for( byte i = 0; i < currentCombi.length; i++ ) {
			currentCombi[i] = i;
		}
	}
}
