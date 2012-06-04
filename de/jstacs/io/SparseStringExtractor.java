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

package de.jstacs.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;

import de.jstacs.data.sequences.annotation.NullSequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotationParser;

/**
 * This {@link StringExtractor} reads the {@link String}s from a {@link java.io.File} as
 * the user asks for the {@link String}. Instances of this class are sparse
 * since they do not store an array of {@link String}s internally.
 * 
 * @author Jens Keilwagen
 */
public class SparseStringExtractor extends AbstractStringExtractor {

	private BufferedReader reader;

	private String current, next, old;

	private StringBuffer help;
	
	/**
	 * A parser for the sequence annotation.
	 */
	protected SequenceAnnotationParser annotationParser;

	/**
	 * A constructor that reads the lines from a file.
	 * 
	 * @param fName
	 *            the name of the {@link java.io.File} to be read from
	 *  
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link java.io.File} could not be found
	 * 
	 * @see SparseStringExtractor#SparseStringExtractor(String, char)
	 */
	public SparseStringExtractor( String fName ) throws IOException, FileNotFoundException {
		this( fName, USUALLY );
	}
	
	/**
	 * A constructor that reads the lines from a file.
	 * 
	 * @param f
	 *            the {@link java.io.File} to be read from
	 *  
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * 
	 * @see SparseStringExtractor#SparseStringExtractor(File, char)
	 */
	public SparseStringExtractor( File f ) throws IOException {
		this( f, USUALLY );
	}
	
	/**
	 * A constructor that reads the lines from a file.
	 * 
	 * @param fName
	 *            the name of the {@link java.io.File} to be read from
	 * @param parser
	 * 			  the {@link SequenceAnnotationParser} for the comment lines 
	 *  
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link java.io.File} could not be found
	 * 
	 * @see SparseStringExtractor#SparseStringExtractor(String, char, SequenceAnnotationParser)
	 */
	public SparseStringExtractor( String fName, SequenceAnnotationParser parser ) throws IOException, FileNotFoundException {
		this( fName, USUALLY, parser );
	}
	
	/**
	 * A constructor that reads the lines from a file and ignores
	 * those starting with the comment character <code>ignore</code>.
	 * 
	 * @param fName
	 *            the name of the {@link java.io.File} to be read from
	 * @param ignore
	 *            the first character of lines that should be treated as
	 *            comments
	 * 
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link java.io.File} could not be found
	 * 
	 * @see SparseStringExtractor#SparseStringExtractor(String, char, String, SequenceAnnotationParser)
	 */
	public SparseStringExtractor( String fName, char ignore ) throws IOException, FileNotFoundException {
		this( fName, ignore, null );
	}
	
	/**
	 * A constructor that reads the lines from a file and ignores
	 * those starting with the comment character <code>ignore</code>.
	 * 
	 * @param f
	 *            the {@link java.io.File} to be read from
	 * @param ignore
	 *            the first character of lines that should be treated as
	 *            comments
	 * 
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * 
	 * @see SparseStringExtractor#SparseStringExtractor(File, char, String, SequenceAnnotationParser)
	 */
	public SparseStringExtractor( File f, char ignore ) throws IOException {
		this( f, ignore, f.getName(), null );
	}

	/**
	 * A constructor that reads the lines from a file and ignores
	 * those starting with the comment character <code>ignore</code>.
	 * 
	 * @param fName
	 *            the name of the {@link java.io.File} to be read from
	 * @param ignore
	 *            the first character of lines that should be treated as
	 *            comments
	 * @param parser
	 * 			  the {@link SequenceAnnotationParser} for the comment lines
	 * 
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link java.io.File} could not be found
	 * 
	 * @see SparseStringExtractor#SparseStringExtractor(String, char, String, SequenceAnnotationParser)
	 */
	public SparseStringExtractor( String fName, char ignore, SequenceAnnotationParser parser ) throws IOException, FileNotFoundException {
		this( fName, ignore, fName, parser );
	}

	/**
	 * A constructor that reads the lines from a file and sets the
	 * annotation of the source to <code>annotation</code>.
	 * 
	 * @param fName
	 *            the name of the {@link java.io.File} to be read from
	 * @param annotation
	 *            the annotation for the source
	 * @param parser
	 * 			  the {@link SequenceAnnotationParser} for the comment lines
	 * 
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link java.io.File} could not be found
	 * 
	 * @see SparseStringExtractor#SparseStringExtractor(String, char, String, SequenceAnnotationParser)
	 */
	public SparseStringExtractor( String fName, String annotation, SequenceAnnotationParser parser ) throws IOException, FileNotFoundException {
		this( fName, USUALLY, annotation, parser );
	}

	/**
	 * A constructor that reads the lines from a file, ignores those
	 * starting with the comment character <code>ignore</code> and sets the
	 * annotation of the source to <code>annotation</code>.
	 * 
	 * @param fName
	 *            the name of the {@link java.io.File} to be read from
	 * @param ignore
	 *            the first character of lines that should be treated as
	 *            comments
	 * @param annotation
	 *            the annotation for the source
	 * @param parser
	 * 			  the {@link SequenceAnnotationParser} for the comment lines
	 * 
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link java.io.File} could not be found
	 * 
	 * @see AbstractStringExtractor#AbstractStringExtractor(char)
	 */
	public SparseStringExtractor( String fName, char ignore, String annotation, SequenceAnnotationParser parser ) throws IOException, FileNotFoundException {
		this( new File( fName ) , ignore, annotation, parser );
	}
	
	/**
	 * A constructor that reads the lines from a file, ignores those
	 * starting with the comment character <code>ignore</code> and sets the
	 * annotation of the source to <code>annotation</code>.
	 * 
	 * @param f
	 *            the {@link java.io.File} to be read from
	 * @param ignore
	 *            the first character of lines that should be treated as
	 *            comments
	 * @param annotation
	 *            the annotation for the source
	 * @param parser
	 * 			  the {@link SequenceAnnotationParser} for the comment lines
	 * 
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 */
	public SparseStringExtractor( File f, char ignore, String annotation, SequenceAnnotationParser parser ) throws IOException {
		this( new FileReader( f ), ignore, annotation, parser );
	}
	
	/**
	 * A constructor that reads the lines from a {@link Reader}, ignores those
	 * starting with the comment character <code>ignore</code> and sets the
	 * annotation of the source to <code>annotation</code>.
	 * 
	 * @param reader
	 *            the {@link Reader} to be read from
	 * @param ignore
	 *            the first character of lines that should be treated as
	 *            comments
	 * @param annotation
	 *            the annotation for the source
	 * @param parser
	 * 			  the {@link SequenceAnnotationParser} for the comment lines
	 * 
	 * @throws IOException
	 *             if the {@link java.io.File} could not be read
	 */
	public SparseStringExtractor( Reader reader, char ignore, String annotation, SequenceAnnotationParser parser ) throws IOException {
		super( ignore );
		this.annotationParser = parser == null ? NullSequenceAnnotationParser.DEFAULT_INSTANCE : parser;
		this.annotationParser.clearAnnotation();
		old = null;
		this.reader = new BufferedReader(reader); 
		help = new StringBuffer();
		next = getNext();
		this.annotation = annotation;
	}

	/* (non-Javadoc)
	 * @see java.util.Enumeration#hasMoreElements()
	 */
	public boolean hasMoreElements() {
		return next != null;
	}

	/* (non-Javadoc)
	 * @see java.util.Enumeration#nextElement()
	 */
	public String nextElement() {
		current = next;
		next = getNext();
		return current;
	}

	private String getNext() {
		String res = null, line;
		int anz = 0;
		try {
			annotationParser.clearAnnotation();
			if( old != null ) {
				annotationParser.addToAnnotation( old );
				anz++;
			}
			old = null;
			while( ( line = reader.readLine() ) != null
					&& ( ( anz < 1 || ignore != FASTA ) && ignorePattern.matcher( line ).matches() ) ) {
				annotationParser.addToAnnotation( line );
				anz++;
				//System.out.println( "skip: " + line );
			}
			if( line == null || ignore != FASTA ) {
				if( anz > 0 && line == null ){
					res = "";
				} else {
					res = line;
				}
			} else {
				while( line != null && !ignorePattern.matcher( line ).matches() ) {
					help.append( line );
					line = reader.readLine();
				}
				if( line != null && ignorePattern.matcher( line ).matches() ) {
					old = line;
				}
				res = help.toString();
				help.delete( 0, help.length() );
			}
		} catch ( IOException e ) {
			throw new RuntimeException( "Forwarding IOException: " + e.getMessage() );
		}
		return res;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#finalize()
	 */
	@Override
	protected void finalize() throws Throwable {
		reader.close();
		super.finalize();
	}
	
	public SequenceAnnotation[] getCurrentSequenceAnnotations() {
		return annotationParser.getCurrentAnnotation();
	}
}