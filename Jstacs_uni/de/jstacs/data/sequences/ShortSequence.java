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

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.SymbolExtractor;

/**
 * This class is for sequences with the alphabet symbols encoded as
 * <code>shorts</code>s and can therefore be used for discrete
 * {@link AlphabetContainer}s with alphabets that use many different symbols.
 * 
 * @author Jens Keilwagen
 */
public class ShortSequence extends SimpleDiscreteSequence {

	private short[] content;

	/**
	 * Creates a new {@link ShortSequence} from an array of <code>short</code>-
	 * encoded alphabet symbols. This constructor is designed for the method
	 * {@link de.jstacs.sequenceScores.statisticalModels.StatisticalModel#emitDataSet(int, int...)}.
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
	 *             not be encoded with <code>short</code>s
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.StatisticalModel#emitDataSet(int, int...)
	 * @see SimpleDiscreteSequence#SimpleDiscreteSequence(AlphabetContainer, SequenceAnnotation[])
	 */
	public ShortSequence( AlphabetContainer alphabetContainer, short[] content ) throws WrongAlphabetException, WrongSequenceTypeException {
		super( alphabetContainer, null );
		if( alphabetContainer.getMaximalAlphabetLength() > Short.MAX_VALUE ) {
			throw new WrongSequenceTypeException();
		}
		this.content = new short[content.length];
		for( int i = 0; i < content.length; i++ ) {
			if( alphabetContainer.isEncodedSymbol( i, content[i] ) ) {
				this.content[i] = content[i];
			} else {
				throw new WrongAlphabetException();
			}
		}
	}

	private ShortSequence( AlphabetContainer cont, SequenceAnnotation[] annotation, short[] content ) throws WrongAlphabetException {
		super( cont, annotation );
		this.content = content;
	}

	/**
	 * Creates a new {@link ShortSequence} from a {@link String} representation
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
	 *             not be encoded with <code>short</code>s
	 * 
	 * @see ShortSequence#ShortSequence(AlphabetContainer, SequenceAnnotation[],
	 *      String, String)
	 */
	public ShortSequence( AlphabetContainer alphabetContainer, String sequence ) throws WrongAlphabetException, WrongSequenceTypeException {
		this( alphabetContainer, null, sequence, alphabetContainer.getDelim() );
	}

	/**
	 * Creates a new {@link ShortSequence} from a {@link String} representation
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
	 *             not be encoded with <code>short</code>s
	 * 
	 * @see ShortSequence#ShortSequence(AlphabetContainer, SequenceAnnotation[],
	 *      SymbolExtractor)
	 */
	public ShortSequence( AlphabetContainer alphabetContainer, SequenceAnnotation[] annotation, String sequence, String delim )
																																throws WrongAlphabetException,
																																WrongSequenceTypeException {
		this( alphabetContainer, annotation, new SymbolExtractor( sequence, delim ) );
	}

	/**
	 * Creates a new {@link ShortSequence} from a {@link SymbolExtractor}.
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
	 *             not be encoded with <code>short</code>s
	 * 
	 * @see SimpleDiscreteSequence#SimpleDiscreteSequence(AlphabetContainer, SequenceAnnotation[])
	 */
	public ShortSequence( AlphabetContainer alphabetContainer, SequenceAnnotation[] annotation, SymbolExtractor extractor )
																															throws WrongAlphabetException,
																															WrongSequenceTypeException {
		super( alphabetContainer, annotation );
		if( alphabetContainer.getMaximalAlphabetLength() > Short.MAX_VALUE ) {
			throw new WrongSequenceTypeException();
		}
		content = new short[extractor.countElements()];
		for( int k = 0; k < content.length; k++ ) {
			content[k] = (short)alphabetContainer.getCode( k, extractor.nextElement() );
		}
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
	protected ShortSequence flatCloneWithoutAnnotation() {
		try {
			return new ShortSequence( alphabetCon, null, content );
		} catch ( Exception doesnothappen ) {
			throw new RuntimeException( doesnothappen );
		}
	}
}
