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

import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.WrongAlphabetException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.ComplementableDiscreteAlphabet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.AbstractStringExtractor;
import de.jstacs.io.SymbolExtractor;

/**
 * This class is an implementation for sequences on one alphabet with length 4.
 * This implementation can be used, for instance, for DNA sequences.
 * 
 * <br>
 * <br>
 * 
 * The symbols are encoded in the bits of the primitive type <code>long</code>,
 * which allows to save 32 symbols in one <code>long</code>. On the one hand an
 * instance of this class is more memory efficient than any other
 * {@link DiscreteSequence}, e.g. {@link ByteSequence}. But on the other hand
 * this class will be a little bit slower when accessing single positions.
 * 
 * @author Jens Keilwagen
 */
public final class SparseSequence extends SimpleDiscreteSequence {

	/**
	 * The maximal number of symbols that can be encoded in one long.
	 */
	private static final int MAX_NUMBER_OF_SYMBOLS_IN_ONE_LONG = 32;

	private static final int LOG_MAX_NUMBER = 5;

	private static final int MAX_NUMBER_MINUS_ONE = 31;

	/**
	 * The length of the sequence.
	 */
	private int length;

	/**
	 * The content of the sequence.
	 */
	private long[] content;

	/**
	 * Creates a new {@link SparseSequence} from a {@link String}
	 * representation.
	 * 
	 * @param alphCon
	 *            the {@link AlphabetContainer}
	 * @param seq
	 *            the sequence as {@link String}
	 * 
	 * @throws WrongSequenceTypeException
	 *             if the {@link AlphabetContainer} is not simple or the
	 *             internal {@link de.jstacs.data.Alphabet} has more than 4
	 *             symbols
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} is not discrete
	 * 
	 * @see SparseSequence#SparseSequence(AlphabetContainer, SymbolExtractor)
	 */
	public SparseSequence( AlphabetContainer alphCon, String seq ) throws WrongSequenceTypeException, WrongAlphabetException {
		this( alphCon, new SymbolExtractor( seq, alphCon.getDelim() ) );
	}

	/**
	 * Creates a new {@link SparseSequence} from a {@link SymbolExtractor}.
	 * 
	 * @param alphCon
	 *            the {@link AlphabetContainer}
	 * @param se
	 *            the {@link SymbolExtractor}
	 * 
	 * @throws WrongSequenceTypeException
	 *             if the {@link AlphabetContainer} is not simple or the
	 *             internal {@link de.jstacs.data.Alphabet} has more than 4
	 *             symbols
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} is not discrete
	 * 
	 * @see SparseSequence#SparseSequence(AlphabetContainer, int, SequenceAnnotation[])
	 */
	public SparseSequence( AlphabetContainer alphCon, SymbolExtractor se ) throws WrongSequenceTypeException, WrongAlphabetException {
		this( alphCon, se.countElements(), null );
		if( !( alphCon.isSimple() && alphCon.getAlphabetLengthAt( 0 ) == 4 ) ) {
			throw new WrongSequenceTypeException();
		}
		DiscreteAlphabet abc = (DiscreteAlphabet)alphCon.getAlphabetAt( 0 );
		for( int j, index = 0, l = 0; l < length; index++ ) {
			content[index] = 0;
			for( j = 0; j < MAX_NUMBER_OF_SYMBOLS_IN_ONE_LONG && l < length; j++, l++ ) {
				content[index] += ( (long)abc.getCode( se.nextElement() ) << ( 2 * j ) );
			}
		}
	}

	/**
	 * This constructor creates an empty sequence of a given length (where empty
	 * means poly-(first-symbol)). This constructor will <b>ONLY</b> be used by
	 * other constructors or for {@link #reverse(int, int)},
	 * {@link #reverseComplement(int, int)} and {@link #complement(int, int)}.
	 * 
	 * @param con
	 *            the {@link AlphabetContainer}
	 * @param length
	 *            the sequence length
	 * @param annotation
	 *            the annotation of the sequence
	 * 
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} is not discrete (forwarded
	 *             from {@link DiscreteSequence})
	 * 
	 * @see DiscreteSequence#DiscreteSequence(AlphabetContainer,
	 *      SequenceAnnotation[])
	 * @see SparseSequence#reverse(int, int)
	 * @see SparseSequence#reverseComplement(int, int)
	 * @see SparseSequence#complement(int, int)
	 */
	private SparseSequence( AlphabetContainer con, int length, SequenceAnnotation[] annotation ) throws WrongAlphabetException {
		super( con, annotation );
		this.length = length;
		content = new long[length / MAX_NUMBER_OF_SYMBOLS_IN_ONE_LONG + ( length % MAX_NUMBER_OF_SYMBOLS_IN_ONE_LONG > 0 ? 1 : 0 )];
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#discreteVal(int)
	 */
	@Override
	public int discreteVal( int pos ) {
		if( pos < 0 || pos >= length ) {
			throw new IndexOutOfBoundsException( "The index " + pos + " is out of bounds [0," + length + "]" );
		} else {
			// TODO make this faster
			return (int)( ( content[pos >> LOG_MAX_NUMBER] >> ( ( pos & MAX_NUMBER_MINUS_ONE ) << 1 ) ) & 3 );
			// slower: return (int) ((content[pos / MAX_NUMBER_OF_SYMBOLS_IN_ONE_LONG] >> (2 * (pos %
			// MAX_NUMBER_OF_SYMBOLS_IN_ONE_LONG))) & 3);
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getLength()
	 */
	@Override
	public int getLength() {
		return length;
	}

	// the next 3 methods are not very time efficient, but they create memory efficient SparseSequences

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#complement(int, int)
	 */
	@Override
	public SparseSequence complement( int start, int end ) throws OperationNotSupportedException {
		if( alphabetCon.isReverseComplementable() ) {
			SparseSequence erg;
			try {
				erg = new SparseSequence( alphabetCon, end - start, null );
			} catch ( WrongAlphabetException e ) {
				// does not happen
				throw new RuntimeException( e.getMessage() );
			}
			ComplementableDiscreteAlphabet compAbc = (ComplementableDiscreteAlphabet)alphabetCon.getAlphabetAt( 0 );
			for( int j, l = start, index = 0; l < end; index++ ) {
				erg.content[index] = 0;
				for( j = 0; j < MAX_NUMBER_OF_SYMBOLS_IN_ONE_LONG && l < end; j++, l++ ) {
					erg.content[index] += ( (long)compAbc.getComplementaryCode( discreteVal( l ) ) << ( 2 * j ) );
				}
			}
			return erg;
		} else {
			throw new OperationNotSupportedException( "The alphabet of sequence has to be complementable." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#reverse(int, int)
	 */
	@Override
	public SparseSequence reverse( int start, int end ) throws OperationNotSupportedException {
		SparseSequence erg;
		try {
			erg = new SparseSequence( alphabetCon, end - start, null );
		} catch ( WrongAlphabetException e ) {
			// does not happen
			throw new RuntimeException( e.getMessage() );
		}
		for( int j, l = end - 1, index = 0; l >= start; index++ ) {
			erg.content[index] = 0;
			for( j = 0; j < MAX_NUMBER_OF_SYMBOLS_IN_ONE_LONG && l >= start; j++, l-- ) {
				erg.content[index] += ( (long)discreteVal( l ) << ( 2 * j ) );
			}
		}
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#reverseComplement(int, int)
	 */
	@Override
	public SparseSequence reverseComplement( int start, int end ) throws OperationNotSupportedException {
		if( rc != null && start == 0 && end == getLength() ) {
			return (SparseSequence)rc;
		} else if( alphabetCon.isReverseComplementable() ) {
			SparseSequence erg;
			try {
				erg = new SparseSequence( alphabetCon, end - start, null );
			} catch ( WrongAlphabetException e ) {
				// does not happen
				throw new RuntimeException( e.getMessage() );
			}
			ComplementableDiscreteAlphabet compAbc = (ComplementableDiscreteAlphabet)alphabetCon.getAlphabetAt( 0 );
			for( int j, l = end - 1, index = 0; l >= start; index++ ) {
				erg.content[index] = 0;
				for( j = 0; j < MAX_NUMBER_OF_SYMBOLS_IN_ONE_LONG && l >= start; j++, l-- ) {
					erg.content[index] += ( (long)compAbc.getComplementaryCode( discreteVal( l ) ) << ( 2 * j ) );
				}
			}
			if( start == 0 && end == getLength() ) {
				rc = erg;
				erg.rc = this;
			}
			return erg;
		} else {
			throw new OperationNotSupportedException( "The alphabet of sequence has to be complementable." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#flatCloneWithoutAnnotation()
	 */
	@Override
	protected SparseSequence flatCloneWithoutAnnotation() {
		try {
			SparseSequence erg = new SparseSequence( alphabetCon, length, null );
			erg.content = content;
			return erg;
		} catch ( Exception doesnothappen ) {
			throw new RuntimeException( doesnothappen );
		}
	}

	/**
	 * This method allows to create a {@link DataSet} containing {@link SparseSequence}s.
	 * 
	 * @param con the {@link AlphabetContainer} for the {@link DataSet} and {@link Sequence}s
	 * @param se the {@link AbstractStringExtractor}s that handle the {@link DataSet} as {@link String} 
	 * 
	 * @return a {@link DataSet} containing {@link SparseSequence}s
	 *  
	 * @throws WrongSequenceTypeException
	 *             if the {@link AlphabetContainer} is not simple or the
	 *             internal {@link de.jstacs.data.Alphabet} has more than 4
	 *             symbols
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} is not discrete
	 * @throws EmptyDataSetException if a {@link DataSet} with 0 (zero) {@link Sequence} should be created
	 */
	public static final DataSet getDataSet( AlphabetContainer con, AbstractStringExtractor... se ) throws WrongSequenceTypeException, WrongAlphabetException, EmptyDataSetException {
		LinkedList<SparseSequence> list = new LinkedList<SparseSequence>();
		SymbolExtractor symEx = new SymbolExtractor( con.getDelim() );
		String s, annot = null;
		for( int i = 0; i < se.length; i++ ) { 
			while( se[i].hasMoreElements() ) {
				s = se[i].nextElement();
				symEx.setStringToBeParsed( s );
				list.add(  new SparseSequence( con, symEx ) );
			}
			if( i == 0 ) {
				annot = se[i].getAnnotation();
			} else {
				annot += ", " + se[i].getAnnotation();
			}
		}
		return new DataSet( annot, list.toArray( new Sequence[0] ) );
	}
}
