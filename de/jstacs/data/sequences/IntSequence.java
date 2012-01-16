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

import de.jstacs.WrongAlphabetException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.SymbolExtractor;

/**
 * This class is for sequences with the alphabet symbols encoded as
 * <code>int</code>s and can therefore be used for discrete
 * {@link AlphabetContainer}s with alphabets that use a huge number of symbols.
 * 
 * @author Jens Keilwagen
 */
public class IntSequence extends SimpleDiscreteSequence {

	private int[] content;

	/**
	 * Creates a new {@link IntSequence} from an array of <code>int</code>-
	 * encoded alphabet symbols. This constructor is designed for the method
	 * {@link de.jstacs.sequenceScores.StatisticalModel#emitDataSet(int, int...)}.
	 * 
	 * @param alphabetContainer
	 *            the {@link AlphabetContainer} for the sequence
	 * @param content
	 *            an array containing the encoded symbols
	 * 
	 * @throws WrongAlphabetException
	 *             if the sequence is not defined over
	 *             <code>alphabetContainer</code>
	 * @throws WrongSequenceTypeException
	 *             if <code>alphabetContainer</code> contains alphabets that can
	 *             not be encoded with <code>int</code>s
	 * 
	 * @see de.jstacs.sequenceScores.StatisticalModel#emitDataSet(int, int...)
	 * @see IntSequence#IntSequence(AlphabetContainer, int[], int, int)
	 */
	public IntSequence( AlphabetContainer alphabetContainer, int... content ) throws WrongAlphabetException, WrongSequenceTypeException {
		this( alphabetContainer, content, 0, content.length );
	}

	/**
	 * Creates a new {@link IntSequence} from a part of the array of
	 * <code>int</code>- encoded alphabet symbols.
	 * 
	 * @param alphabetContainer
	 *            the {@link AlphabetContainer} for the sequence
	 * @param content
	 *            an array containing the encoded symbols
	 * @param start
	 *            the start index in <code>content</code>
	 * @param length
	 *            the length of the part of the array
	 * 
	 * @throws WrongAlphabetException
	 *             if the sequence is not defined over
	 *             <code>alphabetContainer</code>
	 * @throws WrongSequenceTypeException
	 *             if <code>alphabetContainer</code> contains alphabets that can
	 *             not be encoded with <code>int</code>s
	 * 
	 * @see SimpleDiscreteSequence#SimpleDiscreteSequence(AlphabetContainer, SequenceAnnotation[])
	 */
	public IntSequence( AlphabetContainer alphabetContainer, int[] content, int start, int length ) throws WrongAlphabetException,
																									WrongSequenceTypeException {
		super( alphabetContainer.getSubContainer( start, length ), null );
		if( alphabetCon.getMaximalAlphabetLength() > Integer.MAX_VALUE ) {
			throw new WrongSequenceTypeException();
		}
		this.content = new int[length];
		for( int i = 0; i < length; i++ ) {
			if( alphabetCon.isEncodedSymbol( i, content[start + i] ) ) {
				this.content[i] = content[start + i];
			} else {
				throw new WrongAlphabetException();
			}
		}
	}

	/**
	 * Creates a new {@link IntSequence} from a {@link String} representation
	 * using the default delimiter.
	 * 
	 * @param alphabetContainer
	 *            the {@link AlphabetContainer} for the sequence
	 * @param sequence
	 *            a {@link String} representation of the sequence
	 * 
	 * @throws WrongAlphabetException
	 *             if the sequence is not defined over
	 *             <code>alphabetContainer</code>
	 * @throws WrongSequenceTypeException
	 *             if <code>alphabetContainer</code> contains alphabets that can
	 *             not be encoded with <code>int</code>s
	 * 
	 * @see IntSequence#IntSequence(AlphabetContainer, SequenceAnnotation[],
	 *      String, String)
	 */
	public IntSequence( AlphabetContainer alphabetContainer, String sequence ) throws WrongAlphabetException, WrongSequenceTypeException {
		this( alphabetContainer, null, sequence, alphabetContainer.getDelim() );
	}

	/**
	 * Creates a new {@link IntSequence} from a {@link String} representation
	 * using the delimiter <code>delim</code>. Annotations for this sequence are
	 * considered by <code>annotation</code>.
	 * 
	 * @param alphabetContainer
	 *            the {@link AlphabetContainer} for the sequence
	 * @param annotation
	 *            the annotation for this sequence
	 * @param sequence
	 *            a {@link String} representation of the sequence
	 * @param delim
	 *            the delimiter, a {@link String} that separates the symbols
	 * 
	 * @throws WrongAlphabetException
	 *             if the sequence is not defined over
	 *             <code>alphabetContainer</code>
	 * @throws WrongSequenceTypeException
	 *             if <code>alphabetContainer</code> contains alphabets that can
	 *             not be encoded with <code>int</code>s
	 * 
	 * @see IntSequence#IntSequence(AlphabetContainer, SequenceAnnotation[],
	 *      SymbolExtractor)
	 */
	public IntSequence( AlphabetContainer alphabetContainer, SequenceAnnotation[] annotation, String sequence, String delim )
																																throws WrongAlphabetException,
																																WrongSequenceTypeException {
		this( alphabetContainer, annotation, new SymbolExtractor( sequence, delim ) );
	}

	/**
	 * Creates a new {@link IntSequence} from a {@link SymbolExtractor}.
	 * Annotations for this sequence are considered by <code>annotation</code>.
	 * 
	 * @param alphabetContainer
	 *            the {@link AlphabetContainer} for the sequence
	 * @param annotation
	 *            the annotation for this sequence
	 * @param extractor
	 *            the {@link SymbolExtractor}
	 * 
	 * @throws WrongAlphabetException
	 *             if the sequence is not defined over
	 *             <code>alphabetContainer</code>
	 * @throws WrongSequenceTypeException
	 *             if <code>alphabetContainer</code> contains alphabets that can
	 *             not be encoded with <code>int</code>s
	 * 
	 * @see SimpleDiscreteSequence#SimpleDiscreteSequence(AlphabetContainer, SequenceAnnotation[])
	 */
	public IntSequence( AlphabetContainer alphabetContainer, SequenceAnnotation[] annotation, SymbolExtractor extractor )
																															throws WrongAlphabetException,
																															WrongSequenceTypeException {
		super( alphabetContainer, annotation );
		if( alphabetContainer.getMaximalAlphabetLength() > Integer.MAX_VALUE ) {
			throw new WrongSequenceTypeException();
		}
		content = new int[extractor.countElements()];
		for( int k = 0; k < content.length; k++ ) {
			content[k] = (int)alphabetContainer.getCode( k, extractor.nextElement() );
		}
	}

	private IntSequence( AlphabetContainer cont, SequenceAnnotation[] annotation, int[] content ) throws WrongAlphabetException {
		super( cont, annotation );
		this.content = content;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#discreteVal(int)
	 */
	@Override
	public int discreteVal( int pos ) {
		return content[pos];
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getLength()
	 */
	@Override
	public int getLength() {
		return content.length;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#flatCloneWithoutAnnotation()
	 */
	@Override
	protected IntSequence flatCloneWithoutAnnotation() {
		try {
			return new IntSequence( getAlphabetContainer(), null, content );
		} catch ( Exception doesnothappen ) {
			throw new RuntimeException( doesnothappen );
		}
	}
}
