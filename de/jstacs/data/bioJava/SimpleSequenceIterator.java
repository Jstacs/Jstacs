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

package de.jstacs.data.bioJava;

import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;

/**
 * Class that implements the {@link SequenceIterator} interface of BioJava in a
 * simple way, backed by an array of {@link Sequence}s.
 * 
 * @author Jan Grau
 */
public class SimpleSequenceIterator implements SequenceIterator {

	private Sequence[] seqs;

	private int nextIdx;

	/**
	 * Creates a new {@link SimpleSequenceIterator} from an array of
	 * {@link Sequence}s.
	 * 
	 * @param seqs
	 *            the array of {@link Sequence}s
	 */
	public SimpleSequenceIterator( Sequence... seqs ) {
		//TODO clone array reference?
		this.seqs = seqs;
		this.nextIdx = 0;
	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.seq.SequenceIterator#hasNext()
	 */
	public boolean hasNext() {
		return nextIdx < seqs.length;
	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.seq.SequenceIterator#nextSequence()
	 */
	public org.biojava.bio.seq.Sequence nextSequence() throws NoSuchElementException, BioException {
		return seqs[nextIdx++];
	}

	/**
	 * Provides the possibility to reset and afterwards reuse this iterator.
	 */
	public void reset() {
		nextIdx = 0;
	}
}