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

package de.jstacs.data;

import java.util.Arrays;

import de.jstacs.data.alphabets.ComplementableDiscreteAlphabet;
import de.jstacs.data.sequences.IntSequence;

/**
 * This class enumerates over all {@link Sequence}s of a specific
 * {@link AlphabetContainer} and length.
 * 
 * @author Jens Keilwagen
 */
public class DiscreteSequenceEnumerator implements RecyclableSequenceEnumerator {

	private AlphabetContainer con;

	private int[] symbol, stop, compl;

	private int length;

	boolean sparse;

	/**
	 * Creates a new {@link DiscreteSequenceEnumerator} from a given
	 * {@link AlphabetContainer} and a length. <code>sparse</code> determines,
	 * if the enumerator delivers only {@link Sequence}s for which the reverse
	 * complement has not been returned before
	 * 
	 * @param con
	 *            the {@link AlphabetContainer} for the discrete
	 *            {@link Sequence}s
	 * @param length
	 *            the length of the discrete {@link Sequence}s
	 * @param sparse
	 *            indicates if the enumerator delivers only sequences for which
	 *            the reverse complement has not been returned before
	 */
	public DiscreteSequenceEnumerator( AlphabetContainer con, int length, boolean sparse ) {
		if( !con.isDiscrete() ) {
			throw new IllegalArgumentException( "The AlphabetContainer has to be discrete." );
		}
		int l = con.getPossibleLength(), i;
		if( l != 0 && l != length ) {
			throw new IllegalArgumentException( "Please check the length" );
		}
		this.con = con;
		this.length = length;
		symbol = new int[length + 1];
		stop = new int[length + 1];
		for( l = 0; l < length; l++ ) {
			stop[l] = (int)con.getAlphabetLengthAt( l ) - 1;
		}
		stop[length] = 1;
		if( sparse ) {
			if( !con.isSimple() ) {
				throw new IllegalArgumentException( "The enumerator can not be sparse, since it has a non-simple AlphabetContainer." );
			}
			ComplementableDiscreteAlphabet abc = (ComplementableDiscreteAlphabet)con.getAlphabetAt( 0 );
			compl = new int[stop[0] + 1];
			for( i = 0; i <= stop[0]; i++ ) {
				compl[i] = abc.getComplementaryCode( i );
			}
			System.out.println( Arrays.toString( compl ) );
		}
		this.sparse = sparse;
	}

	/* (non-Javadoc)
	 * @see java.util.Enumeration#hasMoreElements()
	 */
	public boolean hasMoreElements() {
		return symbol[length] == 0;
	}

	/* (non-Javadoc)
	 * @see java.util.Enumeration#nextElement()
	 */
	public Sequence<int[]> nextElement() {
		Sequence<int[]> seq;
		try {
			seq = new IntSequence( con, symbol, 0, length );
		} catch ( Exception e ) {
			RuntimeException r = new RuntimeException( e.getMessage() );
			r.setStackTrace( e.getStackTrace() );
			throw r;
		}
		int l, k;
		do {
			// next
			l = 0;
			while( l < length && symbol[l] == stop[l] ) {
				symbol[l] = 0;
				l++;
			}
			symbol[l]++;

			// has been seen ?
			if( sparse ) {
				l = 0;
				k = length - 1;
				while( l < length && symbol[l] == compl[symbol[k]] ) {
					l++;
					k--;
				}
				if( l < length && symbol[l] < compl[symbol[k]] ) {
					l = length;
				}
			}
		} while( symbol[length] == 0 && sparse && l != length );
		return seq;
	}

	public void reset() {
		Arrays.fill( symbol, 0 );
	}
}
