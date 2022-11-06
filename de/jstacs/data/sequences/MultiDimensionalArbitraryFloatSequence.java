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
package de.jstacs.data.sequences;

import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.ArrayHandler;

/**
 * This class is for multidimensional arbitrary float sequences.
 * 
 * @author Jens Keilwagen
 */
public class MultiDimensionalArbitraryFloatSequence extends MultiDimensionalSequence<float[]> {

	/**
	 * This constructor creates an {@link MultiDimensionalArbitraryFloatSequence} from a set of individual {@link Sequence}s.
	 * 
	 * @param seqAn the annotations for the aligned sequences
	 * @param sequence the individual sequences that have been aligned
	 * 
	 * @throws WrongLengthException if the sequences have different lengths
	 * @throws WrongAlphabetException if the sequences have different {@link de.jstacs.data.AlphabetContainer}s 
	 */
	public MultiDimensionalArbitraryFloatSequence( SequenceAnnotation[] seqAn, ArbitraryFloatSequence... sequence ) throws WrongLengthException, WrongAlphabetException {
		super( seqAn, sequence );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#fillContainer(java.lang.Object, int)
	 */
	@Override
	public void fillContainer(float[] container, int pos) {
		for( int s = 0; s < content.length; s++ ) {
			container[s] = (float)content[s].continuousVal( pos );
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getEmptyContainer()
	 */
	@Override
	public float[] getEmptyContainer() {
		return new float[content.length];
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.sequences.MultiDimensionalSequence#getInstance(de.jstacs.data.sequences.annotation.SequenceAnnotation[], de.jstacs.data.Sequence[])
	 */
	@Override
	protected MultiDimensionalArbitraryFloatSequence getInstance( SequenceAnnotation[] seqAn, Sequence... seqs ) throws WrongLengthException, WrongAlphabetException {
		return new MultiDimensionalArbitraryFloatSequence( seqAn, ArrayHandler.cast( ArbitraryFloatSequence.class, seqs ) );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#compareTo(java.lang.Object, java.lang.Object)
	 */
	@Override
	public int compareTo( float[] t1, float[] t2 ) {
		if( t1.length == t2.length ) {
			for( int i = 0; i < t1.length; i++ ) {
				if( t1[i] != t2[i] ) {
					return (int)Math.signum( t1[i] - t2[i] );
				}
			}
			return 0;
		} else {
			return t1.length - t2.length;
		}
	}
}
