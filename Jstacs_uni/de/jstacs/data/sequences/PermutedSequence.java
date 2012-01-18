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

import java.util.Random;

import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;

/**
 * This class is for permuted sequences.
 * 
 * @param <T> the type of each position
 * 
 * @author Jan Grau
 */
public class PermutedSequence<T> extends Sequence.RecursiveSequence<T> {

	private static Random r = new Random();

	private int[] permutation;

	/**
	 * Creates a new {@link PermutedSequence} by shuffling the symbols of a
	 * given {@link Sequence}.
	 * 
	 * @param seq
	 *            the initial sequence
	 * 
	 * @throws WrongAlphabetException
	 *             if the {@link de.jstacs.data.AlphabetContainer} is not simple
	 * 
	 * @see de.jstacs.data.sequences.Sequence.RecursiveSequence#Sequence.RecursiveSequence(de.jstacs.data.AlphabetContainer, de.jstacs.data.sequences.Sequence)
	 */
	public PermutedSequence( Sequence<T> seq ) throws WrongAlphabetException {
		super( seq.getAlphabetContainer(), seq );
		if( !seq.getAlphabetContainer().isSimple() ) {
			throw new WrongAlphabetException( "Alphabet must be simple" );
		}
		permutation = new int[seq.getLength()];
		for( int i = 0; i < permutation.length; i++ ) {
			permutation[i] = i;
		}
		int temp;
		int temp2;
		for( int i = permutation.length - 1; i >= 0; i-- ) {
			temp = r.nextInt( i + 1 );
			temp2 = permutation[i];
			permutation[i] = permutation[temp];
			permutation[temp] = temp2;
		}
	}

	/**
	 * Creates a new {@link PermutedSequence} for a given permutation
	 * 
	 * @param seq
	 *            the original sequence
	 * @param permutation the permutation of the sequence
	 * 
	 * @throws WrongAlphabetException
	 *             if the {@link de.jstacs.data.AlphabetContainer} is not simple
	 * @throws Exception if the length of the permutation does not match that of the sequence
	 * 
	 * @see de.jstacs.data.sequences.Sequence.RecursiveSequence#Sequence.RecursiveSequence(de.jstacs.data.AlphabetContainer, de.jstacs.data.sequences.Sequence)
	 */
	public PermutedSequence( Sequence<T> seq, int[] permutation ) throws WrongAlphabetException, Exception {
		super( seq.getAlphabetContainer(), seq );
		if( !seq.getAlphabetContainer().isSimple() ) {
			throw new WrongAlphabetException( "Alphabet must be simple" );
		}
		if(permutation.length != seq.getLength()){
			throw new Exception("Length of permutation does not match length of sequence");
		}
		this.permutation = permutation.clone();
	}
	
	private PermutedSequence( Sequence<T> seq, SequenceAnnotation[] annotation, int[] perm ) {
		super( seq.getAlphabetContainer(), annotation, seq );
		this.permutation = perm;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.sequences.RecursiveSequence#getIndex(int)
	 */
	@Override
	protected int getIndex( int pos ) {
		return permutation[pos];
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getLength()
	 */
	@Override
	public int getLength() {
		return permutation.length;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#flatCloneWithoutAnnotation()
	 */
	@Override
	protected PermutedSequence<T> flatCloneWithoutAnnotation() {
		try {
			return new PermutedSequence<T>( this.content, null, permutation );
		} catch ( Exception doesnothappen ) {
			throw new RuntimeException( doesnothappen );
		}
	}

}
