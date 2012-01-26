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
package de.jstacs.utils;

import java.io.IOException;
import java.io.OutputStream;

/**
 * This class allows to filter individual lines that are then passed to an internal {@link OutputStream}.
 * 
 * @author Jens Keilwagen
 * 
 * @see SafeOutputStream
 */
public class RegExFilterOutputStream extends OutputStream {

	private OutputStream internal;
	private StringBuffer sb;
	private String regex;
	
	/**
	 * Creates a new {@link RegExFilterOutputStream}.
	 *  
	 * @param o the internal {@link OutputStream}
	 * @param regex the regular expression used to filter each line of the output
	 * 
	 * @see String#matches(String)
	 */
	public RegExFilterOutputStream( OutputStream o, String regex ) {
		internal = o;
		sb = new StringBuffer();
		this.regex = regex;
	}

	@Override
	public void write(int b) throws IOException {
		char c = (char) b;
		sb.append( c );
		if( c == '\n' ) {
			String s = sb.toString();
			if( s.matches( regex ) ) {
				internal.write(s.getBytes());
			}
			sb.delete( 0, sb.length() );
		}
	}
}
