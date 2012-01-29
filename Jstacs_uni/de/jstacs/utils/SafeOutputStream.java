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
 * This class is for any output. For example:
 * <ul>
 * <li>It enables you to hide the output (use
 * <code>new SafeOutputStream( null )</code>).
 * <li>It enables you to write the output to the standard output stream (use
 * <code>new SafeOutputStream( System.out )</code>).
 * <li>It enables you to write the output to a file (use
 * <code>new SafeOutputStream( new FileOutputStream( ... ) )</code>).
 * <li>...
 * </ul>
 * 
 * The user has to close the stream by hand. It is not closed by the
 * <code>finalize</code> method.
 * 
 * @author Jens Keilwagen
 */
public class SafeOutputStream extends OutputStream {

	/**
	 * This method returns an instance of {@link SafeOutputStream} for a given {@link OutputStream}.
	 * 
	 * @param out the {@link OutputStream}
	 * 
	 * @return an instance of {@link SafeOutputStream} for a given {@link OutputStream}.
	 */
	public static SafeOutputStream getSafeOutputStream( OutputStream out ) {
		if( out != null && (out instanceof SafeOutputStream) ) {
			return (SafeOutputStream) out;
		} else {
			return new SafeOutputStream( out );
		}	
	}
	
	/**
	 * This stream can be used as default stream.
	 */
	public static final OutputStream DEFAULT_STREAM = System.out;

	private OutputStream ostream;

	/**
	 * Creates a new {@link SafeOutputStream}.
	 * 
	 * @param ostream
	 *            the internal used {@link OutputStream}
	 */
	private SafeOutputStream( OutputStream ostream ) {
		this.ostream = ostream;
	}

	/**
	 * Closes the {@link SafeOutputStream} by closing the {@link OutputStream}
	 * this stream was constructed of. If this {@link SafeOutputStream} was
	 * constructed using <code>System.out</code> this {@link OutputStream} will
	 * be closed. No outputs via <code>System.out</code> will be possible any
	 * more.
	 * 
	 * @throws IOException
	 *             if something went wrong in closing the output stream
	 */
	@Override
	public void close() throws IOException {
		if( ostream != null ) {
			ostream.flush();
			ostream.close();
		}
	}

	/**
	 * Indicates whether the instance is doing something or not.
	 * 
	 * @return <code>true</code> if the internal {@link OutputStream} is
	 *         <code>null</code>
	 */
	public boolean doesNothing() {
		return ostream == null;
	}

	/* (non-Javadoc)
	 * @see java.io.OutputStream#flush()
	 */
	@Override
	public void flush() throws IOException {
		if( ostream != null ) {
			ostream.flush();
		}
	}

	/**
	 * Returns the internal used {@link OutputStream}.
	 * 
	 * @return the internal used {@link OutputStream}
	 */
	public OutputStream getOutputStream() {
		return ostream;
	}

	/**
	 * Writes a line break.
	 * 
	 * @throws IOException
	 *             if something went wrong
	 * 
	 * @see SafeOutputStream#write(int)
	 */
	public void writeln() throws IOException {
		writeln( "" );
	}

	/**
	 * Writes an {@link Object} and a line break.
	 * 
	 * @param s
	 *            the {@link Object} to be written
	 * 
	 * @throws IOException
	 *             if something went wrong
	 * 
	 * @see SafeOutputStream#write(int)
	 * @see Object#toString()
	 */
	public void writeln( Object s ) throws IOException {
		write( s + "\n" );
	}

	/**
	 * Writes an {@link Object} using {@link Object#toString()}.
	 * 
	 * @param s
	 *            the {@link Object} to be written
	 * 
	 * @throws IOException
	 *             if something went wrong
	 * 
	 * @see SafeOutputStream#write(int)
	 */
	public void write( Object s ) throws IOException {
		if( ostream != null ) {
			ostream.write( s.toString().getBytes() );
		}
	}

	/* (non-Javadoc)
	 * @see java.io.OutputStream#write(int)
	 */
	@Override
	public void write( int b ) throws IOException {
		if( ostream != null ) {
			ostream.write( b );
		}
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#finalize()
	 * 
	 * This method does not close the internal stream.
	 * 
	 */
	@Override
	protected void finalize() throws IOException {
		flush();
	}
}
