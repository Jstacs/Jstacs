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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedList;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotationParser;
import de.jstacs.io.AbstractStringExtractor;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.SymbolExtractor;

/**
 * This class is for any continuous or hybrid sequence.
 * In contrast to {@link ArbitrarySequence}s, the numeric values are represented by <code>float</code>s instead of <code>double</code>s.
 * This may be useful if either the represented data are known to be only of <code>float</code> precision or if <code>double</code> precision is dispensible
 * for the sake of memory consumption.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class ArbitraryFloatSequence extends Sequence<float[]> {

	private float[] content;

	/**
	 * Creates a new {@link ArbitraryFloatSequence} from an array of
	 * <code>float</code>-encoded alphabet symbols. This constructor is
	 * designed for the method
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
	 *             not be encoded with <code>float</code>s
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.StatisticalModel#emitDataSet(int, int...)
	 * @see Sequence#Sequence(AlphabetContainer, SequenceAnnotation[])
	 */
	public ArbitraryFloatSequence( AlphabetContainer alphabetContainer, float[] content ) throws WrongAlphabetException,
																						WrongSequenceTypeException {
		super( alphabetContainer, null );
		this.content = new float[content.length];
		for( int i = 0; i < content.length; i++ ) {
			if( alphabetContainer.isEncodedSymbol( i, content[i] ) ) {
				this.content[i] = content[i];
			} else {
				throw new WrongAlphabetException( "The data of the selected file does not match the entered alphabet: " + content[i] );
			}
		}
	}

	private ArbitraryFloatSequence( AlphabetContainer cont, SequenceAnnotation[] annotation, float[] content ) {
		super( cont, annotation );
		this.content = content;
	}

	/**
	 * Creates a new {@link ArbitraryFloatSequence} from a {@link String}
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
	 *             not be encoded with <code>float</code>s
	 * 
	 * @see ArbitraryFloatSequence#ArbitraryFloatSequence(AlphabetContainer, SequenceAnnotation[], SymbolExtractor)
	 */
	public ArbitraryFloatSequence( AlphabetContainer alphabetContainer, String sequence ) throws WrongAlphabetException,
																					WrongSequenceTypeException {
		this( alphabetContainer, null, new SymbolExtractor( sequence, alphabetContainer.getDelim() ) );
	}

	/**
	 * Creates a new {@link ArbitraryFloatSequence} from a {@link String}
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
	 *             not be encoded with <code>float</code>s
	 * @throws IllegalArgumentException
	 *             if the delimiter is empty and the alphabet container is not
	 *             discrete
	 * 
	 * @see ArbitraryFloatSequence#ArbitraryFloatSequence(AlphabetContainer, SequenceAnnotation[], SymbolExtractor)
	 */
	public ArbitraryFloatSequence( AlphabetContainer alphabetContainer, SequenceAnnotation[] annotation, String sequence, String delim )
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
	 * Creates a new {@link ArbitraryFloatSequence} from a {@link SymbolExtractor}.
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
	 *             not be encoded with <code>float</code>s
	 * 
	 * @see Sequence#Sequence(AlphabetContainer, SequenceAnnotation[])
	 */
	public ArbitraryFloatSequence( AlphabetContainer alphabetContainer, SequenceAnnotation[] annotation, SymbolExtractor extractor )
																																throws WrongAlphabetException,
																																WrongSequenceTypeException {
		super( alphabetContainer, annotation );
		content = new float[extractor.countElements()];
		for( int k = 0; k < content.length; k++ ) {
			content[k] = (float)alphabetContainer.getCode( k, extractor.nextElement() );
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
	protected ArbitraryFloatSequence flatCloneWithoutAnnotation() {
		try {
			return new ArbitraryFloatSequence( alphabetCon, null, content );
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
	public float[] getEmptyContainer() {
		return new float[1];
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#fillContainer(java.lang.Object, int)
	 */
	public void fillContainer( float[] container, int pos ) {
		container[0] = (float)continuousVal( pos );
	}
	
	public int compareTo( float[] t1, float[] t2 ) {
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
	
	/**
	 * This method allows to create a {@link DataSet} containing {@link ArbitraryFloatSequence}s using
	 * a file name. Annotations are parsed by the supplied {@link SequenceAnnotationParser}. The file is 
	 * assumed to be in FastA format.
	 * 
	 * @param con the {@link AlphabetContainer} for the {@link DataSet} and {@link ArbitraryFloatSequence}s
	 * @param filename the file name
	 * @param parser a parser for the annotations of the {@link ArbitraryFloatSequence}s
	 * @return a {@link DataSet} containing {@link ArbitraryFloatSequence}s
	 * 
	 * @throws FileNotFoundException if the file <code>filename</code> could not be found
	 * @throws WrongAlphabetException if the alphabet does not fit the data
	 * @throws WrongSequenceTypeException if the data can not be represented as floats
	 * @throws EmptyDataSetException if not sequences exist in <code>filename</code>
	 * @throws IOException if the file could not be read
	 */
	public static DataSet getDataSet( AlphabetContainer con, String filename, SequenceAnnotationParser parser ) throws FileNotFoundException, WrongAlphabetException, WrongSequenceTypeException, EmptyDataSetException, IOException{
		return getDataSet( con, new SparseStringExtractor( filename, '>', parser ) );
	}
	
	/**
	 * This method allows to create a {@link DataSet} containing {@link ArbitraryFloatSequence}s using
	 * a file name.
	 * 
	 * @param con the {@link AlphabetContainer} for the {@link DataSet} and {@link ArbitraryFloatSequence}s
	 * @param filename the file name
	 * @return a {@link DataSet} containing {@link ArbitraryFloatSequence}s
	 * @throws FileNotFoundException if the file <code>filename</code> could not be found
	 * @throws WrongAlphabetException if the alphabet does not fit the data
	 * @throws WrongSequenceTypeException if the data can not be represented as floats
	 * @throws EmptyDataSetException if not sequences exist in <code>filename</code>
	 * @throws IOException if the file could not be read
	 */
	public static DataSet getDataSet( AlphabetContainer con, String filename ) throws FileNotFoundException, WrongAlphabetException, WrongSequenceTypeException, EmptyDataSetException, IOException{
		return getDataSet( con, new SparseStringExtractor( filename ) );
	}
	
	/**
	 * This method allows to create a {@link DataSet} containing {@link ArbitraryFloatSequence}s.
	 * 
	 * @param con the {@link AlphabetContainer} for the {@link DataSet} and {@link ArbitraryFloatSequence}s
	 * @param se the {@link AbstractStringExtractor}s that handle the {@link DataSet} as {@link String} 
	 * 
	 * @return a {@link DataSet} containing {@link ArbitraryFloatSequence}s
     * @throws WrongAlphabetException if the alphabet does not fit the data
	 * @throws WrongSequenceTypeException if the data can not be represented as floats
	 *  
	 * @throws EmptyDataSetException if a {@link DataSet} with 0 (zero) {@link ArbitraryFloatSequence}s should be created
	 */
	public static DataSet getDataSet( AlphabetContainer con, AbstractStringExtractor... se ) throws WrongAlphabetException, WrongSequenceTypeException, EmptyDataSetException{
		LinkedList<ArbitraryFloatSequence> list = new LinkedList<ArbitraryFloatSequence>();
		SymbolExtractor symEx = new SymbolExtractor( con.getDelim() );
		String s, annot = null;
		for( int i = 0; i < se.length; i++ ) { 
			while( se[i].hasMoreElements() ) {
				SequenceAnnotation[] anns = se[i].getCurrentSequenceAnnotations();
				s = se[i].nextElement();
				symEx.setStringToBeParsed( s );
				list.add(  new ArbitraryFloatSequence( con, anns, symEx ) );
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
