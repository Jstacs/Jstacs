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
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */
package de.jstacs.data;

import java.util.HashSet;
import java.util.Iterator;

import javax.naming.OperationNotSupportedException;

/**
 * Class for an {@link RecyclableSequenceEnumerator} of {@link Sequence}s that enumerates all k-mers that exist in a given {@link DataSet}, optionally ignoring reverse complements.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class DataSetKMerEnumerator implements RecyclableSequenceEnumerator {

	private Iterator<Sequence> it;
	private HashSet<Sequence> current;
	
	/**
	 * Constructs a new SampleKMerEnumerator from a {@link DataSet} <code>data</code> by extracting all k-mers.
	 * @param data the data to extract the k-mers from
	 * @param k the length of the k-mers
	 * @param eliminateRevComp whether to ignore the reverse complement of already existing k-mers
	 * @throws OperationNotSupportedException if a {@link Sequence} in <code>data</code> does not support to compute the reverse complement but <code>eliminateRevComp</code> is true
	 */
	public DataSetKMerEnumerator(DataSet data, int k, boolean eliminateRevComp) throws OperationNotSupportedException{

		int f = eliminateRevComp ? 2 : 1;
		
		current = new HashSet<Sequence>(f * data.getMaximalElementLength());
		
		for( int j = 0; j < data.getNumberOfElements(); j++ ) {
			Sequence seq = data.getElementAt( j );
			for( int l = 0; l<=seq.getLength() - k; l++ ){
				Sequence temp = seq.getSubSequence( l, k );
				if( ! ( current.contains( temp ) || ( eliminateRevComp && current.contains( temp.reverseComplement() ) ) ) ){
					current.add( temp );
				}
			}
		}
		reset();
	}
	
	public void reset(){
		it = current.iterator();
	}
		
	public boolean hasMoreElements() {
		return it.hasNext();
	}

	public Sequence nextElement() {
		return it.next();
	}	
}
