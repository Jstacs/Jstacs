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
import java.util.StringTokenizer;

/**
 * This class implements the reader that extracts {@link String}s from either a
 * {@link File} or a {@link String}. Internally the {@link String}s are
 * extracted and stored in an array when creating a new instance.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class StringExtractor extends AbstractStringExtractor {

	private String[] strs;

	private int last, idx, anz;
	
	private boolean first;

	private StringBuffer current;

	private StringExtractor( int initSize, char ignore ) throws IllegalArgumentException {
		super( ignore );
		if( initSize < 1 ) {
			throw new IllegalArgumentException( "The initSize is too small." );
		}
		last = idx = 0;
		strs = new String[initSize];
	}

	/**
	 * A constructor that reads the lines from <code>file</code>.
	 * 
	 * @param file
	 *            the {@link File} to be read from
	 * @param initSize
	 *            the initial number of lines that can be handled
	 * 
	 * @throws IOException
	 *             if the {@link File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link File} could not be found
	 * 
	 * @see StringExtractor#StringExtractor(File, int, char)
	 */
	public StringExtractor( File file, int initSize ) throws IOException, FileNotFoundException {
		this( file, initSize, USUALLY );
	}

	/**
	 * A constructor that reads the lines from <code>file</code> and ignores
	 * those starting with the comment character <code>ignore</code>.
	 * 
	 * @param file
	 *            the {@link File} to be read from
	 * @param initSize
	 *            the initial number of lines that can be handled
	 * @param ignore
	 *            the first character of lines that should be treated as
	 *            comments
	 * 
	 * @throws IOException
	 *             if the {@link File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link File} could not be found
	 * 
	 * @see StringExtractor#StringExtractor(File, int, char, String)
	 */
	public StringExtractor( File file, int initSize, char ignore ) throws IOException, FileNotFoundException {
		this( file, initSize, ignore, file.getName() );
	}

	/**
	 * A constructor that reads the lines from <code>file</code> and sets the
	 * annotation of the source to <code>annotation</code>.
	 * 
	 * @param file
	 *            the {@link File} to be read from
	 * @param initSize
	 *            the initial number of lines that can be handled
	 * @param annotation
	 *            the annotation for the source
	 * 
	 * @throws IOException
	 *             if the {@link File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link File} could not be found
	 * 
	 * @see StringExtractor#StringExtractor(File, int, char, String)
	 */
	public StringExtractor( File file, int initSize, String annotation ) throws IOException, FileNotFoundException {
		this( file, initSize, USUALLY, annotation );
	}

	/**
	 * A constructor that reads the lines from <code>file</code>, ignores those
	 * starting with the comment character <code>ignore</code> and sets the
	 * annotation of the source to <code>annotation</code>.
	 * 
	 * @param file
	 *            the {@link File} to be read from
	 * @param initSize
	 *            the initial number of lines that can be handled
	 * @param ignore
	 *            the first character of lines that should be treated as
	 *            comments
	 * @param annotation
	 *            the annotation for the source
	 * 
	 * @throws IOException
	 *             if the {@link File} could not be read
	 * @throws FileNotFoundException
	 *             if the {@link File} could not be found
	 */
	public StringExtractor( File file, int initSize, char ignore, String annotation ) throws IOException, FileNotFoundException {
		this( initSize, ignore );
		BufferedReader reader = new BufferedReader( new FileReader( file ) );
		String str = null;
		anz = 0;
		first = true;
		while( ( str = reader.readLine() ) != null ) {
			insert( str );
		}
		if( anz > 0 ) {
			strs[last++] = ignore == FASTA ? (current == null ?  "" : current.toString()) : "";
		}
		if( current != null && current.length() > 0 ) {
			strs[last++] = current.toString();
			current.delete( 0, current.length() );
		}
		reader.close();
		this.annotation = annotation;
	}

	/**
	 * A constructor that reads the lines from a {@link String}
	 * <code>content</code> and sets the annotation of the source to
	 * <code>annotation</code>.
	 * 
	 * @param content
	 *            the complete {@link String} with all lines
	 * @param initSize
	 *            the initial number of lines that can be handled
	 * @param annotation
	 *            some annotation for the content
	 * 
	 * @see StringExtractor#StringExtractor(String, int, char, String)
	 */
	public StringExtractor( String content, int initSize, String annotation ) {
		this( content, initSize, USUALLY, annotation );
	}

	/**
	 * A constructor that reads the lines from a {@link String}
	 * <code>content</code>, ignores those starting with the comment character
	 * <code>ignore</code> and sets the annotation of the source to
	 * <code>annotation</code>.
	 * 
	 * @param content
	 *            the complete {@link String} with all lines
	 * @param initSize
	 *            the initial number of lines that can be handled
	 * @param ignore
	 *            the first character of lines that should be treated as
	 *            comments
	 * @param annotation
	 *            some annotation for the content
	 */
	public StringExtractor( String content, int initSize, char ignore, String annotation ) {
		this( initSize, ignore );
		StringTokenizer tok = new StringTokenizer( content, "\n" );
		anz = 0;
		first = true;
		while( tok.hasMoreTokens() ) {
			insert( tok.nextToken() );
		}
		if( anz > 0 ) {
			strs[last++] = ignore == FASTA ? (current == null ?  "" : current.toString()) : "";
		} else if( current != null && current.length() > 0 ) {
			strs[last++] = current.toString();
			current.delete( 0, current.length() );
		}
		this.annotation = annotation;
	}

	private void expand() {
		String[] temp = new String[2 * strs.length];
		System.arraycopy( strs, 0, temp, 0, strs.length );
		strs = temp;
		temp = null;
	}

	/* (non-Javadoc)
	 * @see java.util.Enumeration#nextElement()
	 */
	public String nextElement() {
		if( idx >= last ) {
			throw new IndexOutOfBoundsException();
		}
		return strs[idx++];
	}

	/* (non-Javadoc)
	 * @see java.util.Enumeration#hasMoreElements()
	 */
	public boolean hasMoreElements() {
		return idx < last;
	}

	/**
	 * Returns the number of {@link String}s that have been read.
	 * 
	 * @return the number of {@link String}s that have been read
	 */
	public int getNumberOfElements() {
		return last;
	}

	/**
	 * Returns {@link String} number <code>idx</code> that has been extracted.
	 * 
	 * @param idx
	 *            the number of the {@link String}
	 * 
	 * @return {@link String} number <code>idx</code>
	 */
	public String getElement( int idx ) {
		if( idx < last ) {
			return strs[idx];
		} else {
			throw new IndexOutOfBoundsException();
		}
	}

	private void insert( String str ) {
		if( ignorePattern.matcher( str ).matches() ) {
			if( ignore == FASTA ) {
				if( !first ) {
					if( current != null ) {
						strs[last++] = current.toString();
						current.delete( 0, current.length() );
					} else {
						strs[last++] = "";
					}
				}
				anz = 1;
			} else {
				//ignore;
				anz++;
			}
		} else {
			anz = 0;
			if( ignore == FASTA ) {
				if( current == null ) {
					current = new StringBuffer();
				}
				current.append( str );
			} else {
				strs[last++] = str;
			}
		}
		if( last == strs.length ) {
			expand();
		}
		first = false;
	}
}
