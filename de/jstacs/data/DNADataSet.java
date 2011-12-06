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

import java.io.FileNotFoundException;
import java.io.IOException;

import de.jstacs.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.sequences.annotation.SequenceAnnotationParser;
import de.jstacs.io.SparseStringExtractor;

/**
 * This class exist for convenience to allow the user an easy creation of {@link DataSet}s of DNA {@link Sequence}s.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class DNADataSet extends DataSet {
	
	static {
		try {
			DNA_ALPHABETCONTAINER = new AlphabetContainer( new DNAAlphabet() );
		} catch (Exception doesNotHappen) {
			System.err.println( "Unexpected Error:" );
			doesNotHappen.printStackTrace();
			System.exit(1);
		}		
	}
	
	private static AlphabetContainer DNA_ALPHABETCONTAINER;//TODO final Jan
	
	/**
	 * Creates a new sample of DNA sequence from a FASTA file with file name <code>fName</code>.
	 * 
	 * @param fName the file name
	 * 
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link java.io.File} could not be found
	 * @throws WrongAlphabetException
	 *             if the DNA {@link AlphabetContainer} is not suitable
	 * @throws EmptyDataSetException
	 *             if the {@link DataSet} would be empty
	 * @throws WrongLengthException
	 *             never happens (forwarded from
	 *             {@link DataSet#Sample(AlphabetContainer, de.jstacs.io.AbstractStringExtractor, String, int)}
	 *             )
	 * 
	 * @see de.jstacs.io.AbstractStringExtractor#FASTA
	 * @see DNADataSet#DNASample(String, char)
	 */
	public DNADataSet( String fName ) throws FileNotFoundException, WrongAlphabetException, EmptyDataSetException, WrongLengthException, IOException {
		this( fName, SparseStringExtractor.FASTA );//TODO Jan
	}
	
	/**
	 * Creates a new sample of DNA sequence from a file with file name <code>fName</code>.
	 * 
	 * @param fName the file name
	 * @param ignore the first character of lines that should be treated as comments
	 * 
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link java.io.File} could not be found
	 * @throws WrongAlphabetException
	 *             if the DNA {@link AlphabetContainer} is not suitable
	 * @throws EmptyDataSetException
	 *             if the {@link DataSet} would be empty
	 * @throws WrongLengthException
	 *             never happens (forwarded from
	 *             {@link DataSet#Sample(AlphabetContainer, de.jstacs.io.AbstractStringExtractor, String, int)}
	 *             )
	 * 
	 * @see SparseStringExtractor
	 * @see de.jstacs.io.AbstractStringExtractor#FASTA
	 * @see de.jstacs.io.AbstractStringExtractor#USUALLY
	 * @see DNADataSet#DNASample(String, char, SequenceAnnotationParser)
	 */
	public DNADataSet( String fName, char ignore ) throws FileNotFoundException, WrongAlphabetException, EmptyDataSetException, WrongLengthException, IOException {
		this( fName, ignore, null );
	}
	
	/**
	 * Creates a new sample of DNA sequence from a file with file name <code>fName</code> using the given <code>parser</code>.
	 * 
	 * @param fName the file name
	 * @param ignore the first character of lines that should be treated as comments
	 * @param parser the parser for the {@link de.jstacs.data.sequences.annotation.SequenceAnnotation}
	 * 
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link java.io.File} could not be found
	 * @throws WrongAlphabetException
	 *             if the DNA {@link AlphabetContainer} is not suitable
	 * @throws EmptyDataSetException
	 *             if the {@link DataSet} would be empty
	 * @throws WrongLengthException
	 *             never happens (forwarded from
	 *             {@link DataSet#Sample(AlphabetContainer, de.jstacs.io.AbstractStringExtractor, String, int)}
	 *             )
	 * 
	 * @see SparseStringExtractor
	 * @see de.jstacs.io.AbstractStringExtractor#FASTA
	 * @see de.jstacs.io.AbstractStringExtractor#USUALLY
	 * @see DataSet#Sample(AlphabetContainer, de.jstacs.io.AbstractStringExtractor)
	 */
	public DNADataSet( String fName, char ignore, SequenceAnnotationParser parser ) throws FileNotFoundException, WrongAlphabetException, EmptyDataSetException, WrongLengthException, IOException {
		super( DNA_ALPHABETCONTAINER, new SparseStringExtractor( fName, ignore, parser ) );
	}
}