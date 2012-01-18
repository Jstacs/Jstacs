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
