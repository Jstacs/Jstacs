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
import de.jstacs.data.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.SymbolExtractor;

/**
 * This class is for any continuous or hybrid sequence.
 * 
 * @author Jens Keilwagen
 */
public class ArbitrarySequence extends Sequence<double[]> {

	private double[] content;

	/**
	 * Creates a new {@link ArbitrarySequence} from an array of
	 * <code>double</code>-encoded alphabet symbols. This constructor is
	 * designed for the method
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
	 *             not be encoded with <code>double</code>s
	 * 
	 * @see de.jstacs.sequenceScores.StatisticalModel#emitDataSet(int, int...)
	 * @see Sequence#Sequence(AlphabetContainer, SequenceAnnotation[])
	 */
	public ArbitrarySequence( AlphabetContainer alphabetContainer, double[] content ) throws WrongAlphabetException,
																						WrongSequenceTypeException {
		super( alphabetContainer, null );
		this.content = new double[content.length];
		for( int i = 0; i < content.length; i++ ) {
			if( alphabetContainer.isEncodedSymbol( i, content[i] ) ) {
				this.content[i] = content[i];
			} else {
				throw new WrongAlphabetException();
			}
		}
	}

	private ArbitrarySequence( AlphabetContainer cont, SequenceAnnotation[] annotation, double[] content ) {
		super( cont, annotation );
		this.content = content;
	}

	/**
	 * Creates a new {@link ArbitrarySequence} from a {@link String}
	 * representation using the default delimiter.
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
	 *             not be encoded with <code>double</code>s
	 * 
	 * @see ArbitrarySequence#ArbitrarySequence(AlphabetContainer,
	 *      SequenceAnnotation[], SymbolExtractor)
	 */
	public ArbitrarySequence( AlphabetContainer alphabetContainer, String sequence ) throws WrongAlphabetException,
																					WrongSequenceTypeException {
		this( alphabetContainer, null, new SymbolExtractor( sequence, alphabetContainer.getDelim() ) );
	}

	/**
	 * Creates a new {@link ArbitrarySequence} from a {@link String}
	 * representation using the delimiter <code>delim</code>. Annotations for
	 * this sequence are considered by <code>annotation</code>.
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
	 *             not be encoded with <code>double</code>s
	 * @throws IllegalArgumentException
	 *             if the delimiter is empty and the alphabet container is not
	 *             discrete
	 * 
	 * @see ArbitrarySequence#ArbitrarySequence(AlphabetContainer,
	 *      SequenceAnnotation[], SymbolExtractor)
	 */
	public ArbitrarySequence( AlphabetContainer alphabetContainer, SequenceAnnotation[] annotation, String sequence, String delim )
																																	throws WrongAlphabetException,
																																	WrongSequenceTypeException,
																																	IllegalArgumentException {
		this( alphabetContainer, annotation, new SymbolExtractor( sequence, checkDelim( alphabetContainer, delim ) ) );
	}

	private static String checkDelim( AlphabetContainer alphabetContainer, String delim ) throws IllegalArgumentException {
		if( !alphabetContainer.isDiscrete() && delim.length() == 0 ) {
			throw new IllegalArgumentException( "The emtpy delimiter is forbidden for non discrete AlphabetContainers." );
		}
		return delim;
	}

	/**
	 * Creates a new {@link ArbitrarySequence} from a {@link SymbolExtractor}.
	 * Annotations for this sequence are considered by <code>annotation</code>.
	 * 
	 * @param alphabetContainer
	 *            the alphabet container for the sequence
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
	 *             not be encoded with <code>double</code>s
	 * 
	 * @see Sequence#Sequence(AlphabetContainer, SequenceAnnotation[])
	 */
	public ArbitrarySequence( AlphabetContainer alphabetContainer, SequenceAnnotation[] annotation, SymbolExtractor extractor )
																																throws WrongAlphabetException,
																																WrongSequenceTypeException {
		super( alphabetContainer, annotation );
		content = new double[extractor.countElements()];
		for( int k = 0; k < content.length; k++ ) {
			content[k] = alphabetContainer.getCode( k, extractor.nextElement() );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#continuousVal(int)
	 */
	@Override
	public double continuousVal( int pos ) {
		return content[pos];
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#discreteVal(int)
	 */
	@Override
	public int discreteVal( int pos ) {
		return toDiscrete( pos, content[pos] );
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
	protected ArbitrarySequence flatCloneWithoutAnnotation() {
		try {
			return new ArbitrarySequence( alphabetCon, null, content );
		} catch ( Exception doesnothappen ) {
			throw new RuntimeException( doesnothappen );
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#isMultiDimensional()
	 */
	public boolean isMultiDimensional()  {
		return false;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getEmptyContainer()
	 */
	public double[] getEmptyContainer() {
		return new double[1];
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#fillContainer(java.lang.Object, int)
	 */
	public void fillContainer( double[] container, int pos ) {
		container[0] = continuousVal( pos );
	}
	
	public int compareTo( double[] t1, double[] t2 ) {
		if( t1.length == 1 && t2.length == 1 ) {
			return (int)Math.signum( t1[0] - t2[0] );
		} else {
			return t1.length - t2.length;
		}
	}
	
	protected Object getEmptyRepresentation() {
		return new StringBuffer();
	}
	
	protected void addToRepresentation( Object representation, int pos, String delim ) {
		((StringBuffer)representation).append( alphabetCon.getSymbol( pos, continuousVal( pos ) ) + delim );
	}
	protected String getStringRepresentation( Object representation ) {
		return representation.toString();
	}
	
	protected int hashCodeForPos( int pos ) {
		long bits = Double.doubleToLongBits( continuousVal( pos ) );
		return (int)( bits ^ ( bits >>> 32 ) );
	}
}
