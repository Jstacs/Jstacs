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
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Enumeration;

import de.jstacs.WrongAlphabetException;
import de.jstacs.data.AlphabetContainer;

/**
 * This class enables you to extract elements (symbols) from a given
 * {@link String} similar to a {@link java.util.StringTokenizer}.
 * 
 * <br>
 * <br>
 * 
 * The class has some special functionalities that are not given by a
 * {@link java.util.StringTokenizer}:
 * <ol>
 * <li>It enables you to recycle used objects.
 * <li>It handles the delimiter &quot;&quot; as <code>null</code>. (In this case
 * the characters of the {@link String} are returned.)
 * </ol>
 * 
 * @author Jens Keilwagen
 */
public class SymbolExtractor implements Enumeration<String> {

	private String delim, string;

	private boolean simple;

	private int index;

	// private StringTokenizer st;

	/**
	 * Creates a new {@link SymbolExtractor} using <code>delim</code> as
	 * delimiter. Before invoking other methods one has to use
	 * {@link #setStringToBeParsed(String)}.
	 * 
	 * @param delim
	 *            the delimiter
	 * 
	 * @see SymbolExtractor#setStringToBeParsed(String)
	 * @see SymbolExtractor#SymbolExtractor(String, String)
	 */
	public SymbolExtractor( String delim ) {
		this( null, delim );
	}

	/**
	 * Creates a new {@link SymbolExtractor} using <code>delim</code> as
	 * delimiter and <code>string</code> as the {@link String} to be parsed.
	 * 
	 * @param string
	 *            the {@link String} to be parsed
	 * @param delim
	 *            the delimiter
	 */
	public SymbolExtractor( String string, String delim ) {
		this.delim = delim;
		simple = delim == null | delim.length() == 0;
		setStringToBeParsed( string );
	}

	/**
	 * Sets a new {@link String} to be parsed.
	 * 
	 * @param string
	 *            the {@link String} to be parsed
	 */
	public void setStringToBeParsed( String string ) {
		this.string = string;
		if( simple ) {
			// st = null;
			index = 0;
		} else {
			// st = new StringTokenizer( string, delim );
			index = -1;
		}
	}

	/**
	 * Counts the number of elements (symbols) that can be received initially.
	 * 
	 * @return the number of elements (symbols) that can be received initially
	 */
	public int countElements() {
		if( simple ) {
			return string.length();
		} else {
			int current = -1, prev, no = 0;
			while( current < string.length() ) {
				do {
					prev = current;
					current = string.indexOf( delim, prev + 1 );
				} while( current - prev == 1 );
				if( current == -1 ) {
					return ( string.length() - prev > 1 ) ? ++no : no;
				} else {
					no++;
				}
			}
			return -1; // does not happen
			// return st.countTokens();
		}
	}

	/* (non-Javadoc)
	 * @see java.util.Enumeration#hasMoreElements()
	 */
	public boolean hasMoreElements() {
		if( simple ) {
			return index < string.length();
		} else {
			if( index < string.length() ) {
				int current = index, prev;
				do {
					prev = current;
					current = string.indexOf( delim, prev + 1 );
				} while( current - prev == 1 );
				index = prev;
				if( current == -1 ) {
					return string.length() - prev > 1;
				} else {
					return true;
				}
			} else {
				return false;
			}
			// return st.hasMoreElements();
		}
	}

	/* (non-Javadoc)
	 * @see java.util.Enumeration#nextElement()
	 */
	public String nextElement() {
		if( simple ) {
			return String.valueOf( string.charAt( index++ ) );
		} else {
			int current = index, prev;
			do {
				prev = current;
				current = string.indexOf( delim, prev + 1 );
			} while( current - prev == 1 );
			if( current == -1 ) {
				index = string.length();
			} else {
				index = current;
			}
			return string.substring( prev + 1, index );
			// return st.nextToken();
		}
	}

	/**
	 * This method allows the user to filter the content of a
	 * {@link java.io.File} using a given {@link AlphabetContainer} and a
	 * minimal sequence length. This method is useful to filter, for instance,
	 * data sets with ambiguous nucleotides and/or short sequences. Lines that
	 * can not be parsed with respect to the {@link AlphabetContainer} will be
	 * masked by the character <code>ignore</code>.
	 * 
	 * @param inFile
	 *            the path/name of the input {@link java.io.File}
	 * @param ignore
	 *            the char for comment lines
	 * @param con
	 *            the {@link AlphabetContainer}
	 * @param minLength
	 *            the minimal length of a sequence in the output
	 *            {@link java.io.File}
	 * @param outFile
	 *            the path/name of output {@link java.io.File}
	 * 
	 * @return the number of discarded lines of the input {@link java.io.File}
	 * 
	 * @throws IOException
	 *             if something with the handling of the {@link java.io.File}s
	 *             went wrong
	 */
	public static int filter( String inFile, char ignore, AlphabetContainer con, int minLength, String outFile ) throws IOException {
		BufferedReader r = new BufferedReader( new FileReader( inFile ) );
		BufferedWriter w = new BufferedWriter( new FileWriter( outFile ) );
		String line, ign = "" + ignore;
		SymbolExtractor se = new SymbolExtractor( con.getDelim() );
		int i, lines = 0;
		while( ( line = r.readLine() ) != null ) {
			if( line.startsWith( ign ) ) {
				w.write( line );
			} else {
				se.setStringToBeParsed( line );
				i = 0;
				try {
					while( se.hasMoreElements() ) {
						con.getCode( i, se.nextElement() );
						i++;
					}
					i -= minLength;
				} catch ( WrongAlphabetException wae ) {
					i = -1;
				}
				if( i < 0 ) {
					w.write( ign + line );
					lines++;
				} else {
					w.write( line );
				}
			}

			w.newLine();
		}
		w.close();
		r.close();
		return lines;
	}
}