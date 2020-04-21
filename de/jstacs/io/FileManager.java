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
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * This class is for handling {@link File}s. The most important methods of this
 * class are for writing a {@link StringBuffer} to or reading it from a
 * {@link File}. This is useful for all objects that implement
 * {@link de.jstacs.Storable} and should be saved in or loaded from a
 * {@link File}. Additionally, two write and read methods based on {@link OutputStream}
 * and {@link InputStream} have been implemented to allow for handling compressed data via, e.g.,
 * {@link java.util.zip.GZIPOutputStream} and {@link java.util.zip.GZIPInputStream}.
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.Storable
 * @see File
 */
public class FileManager {

	private FileManager() {}

	/**
	 * This buffer is used to copy {@link File}s.
	 */
	private static byte[] buffer = new byte[100000];
	
	/**
	 * This method copies all {@link File}s and directories, if selected, from a
	 * <code>source</code> {@link File}, i.e. directory, to a
	 * <code>target</code> {@link File}, i.e. directory, that are accepted by the f
	 * {@link FileFilter} <code>filter</code>.
	 * 
	 * @param source
	 *            the source directory denoted as {@link File}
	 * @param target
	 *            the target directory denoted as {@link File}
	 * @param filter
	 *            a {@link FileFilter} for the {@link File}s that enables the
	 *            user to copy only specific {@link File}s
	 * @param newer
	 * 			  a switch allowing to copy only files from the source directory that are newer than those in the target directory
	 * 
	 * @return the number of copied {@link File}s
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>source</code> and <code>target</code> are not
	 *             directories
	 * @throws IOException
	 *             if something went wrong while copying the {@link File}s and
	 *             directories
	 * 
	 * @see File
	 * @see FileFilter
	 */
	public static int copy( File source, File target, FileFilter filter, boolean newer ) throws IllegalArgumentException,
			IOException {
		if( !source.isDirectory() || (target.exists() && !target.isDirectory()) ) {
			throw new IllegalArgumentException( "The source and the target have to be directories. (" + source.getAbsolutePath() + ", " + target.getAbsolutePath() );
		}
		File[] files = filter == null ? source.listFiles() : source.listFiles( filter );
		int anz = 0;
		for( File f : files ) {
			if( f.isDirectory() ) {
				anz += copy( f, new File( target.getAbsolutePath() + "/" + f.getName() ), filter, newer );
			} else {
				String current = target.getAbsolutePath() + "/" + f.getName();
				File h = new File( current );
				if( !h.exists() || !newer || f.lastModified() > h.lastModified() ) {
					if( !target.exists() ) {
						target.mkdirs();
					}
					copy( f.getAbsolutePath(), current );
					//System.out.println(current);
					anz++;
				}
			}
		}
		return anz;
	}

	/**
	 * This method copies a {@link File} in a faster manner.
	 * 
	 * @param from
	 *            the {@link File} name of the original file
	 * @param to
	 *            the {@link File} name of the copied file
	 * 
	 * @throws IOException
	 *             if something went wrong
	 * 
	 * @see FileManager#copy(String, String, byte[])
	 */
	public static void copy( String from, String to ) throws IOException {
		copy( from, to, buffer );
	}

	/**
	 * This method copies a {@link File} in a faster manner using a specified
	 * buffer.
	 * 
	 * @param from
	 *            the {@link File} name of the original file
	 * @param to
	 *            the {@link File} name of the copied file
	 * @param buffer
	 *            an array for reading the content of the original {@link File},
	 *            the size of the array determines how many <code>byte</code>s
	 *            are read at once
	 * 
	 * @throws IOException
	 *             if something went wrong
	 */
	public static synchronized void copy( String from, String to, byte[] buffer ) throws IOException {
		FileInputStream in = null;
		FileOutputStream out = null;
		try {
			in = new FileInputStream( from );
			out = new FileOutputStream( to );
			int amountRead;
			while( ( amountRead = in.read( buffer ) ) > -1 ) {
				out.write( buffer, 0, amountRead );
			}
		} finally {
			if( in != null ) {
				in.close();
			}
			if( out != null ) {
				out.close();
			}
		}
		new File( to ).setLastModified( new File( from ).lastModified() );
	}

	/**
	 * This method reads a {@link StringBuffer} from a file with the given filen name.
	 * 
	 * @param fName
	 *            the name of the file to be read
	 * 
	 * @return a {@link StringBuffer} with the content of the {@link File}
	 * 
	 * @throws IOException
	 *             if something went wrong with the {@link File}
	 * 
	 * @see #writeFile(File, CharSequence)
	 * @see #read(Reader)
	 * @see FileReader
	 */
	public static StringBuffer readFile( String fName ) throws IOException {
		return read( new FileReader( fName ) );
	}

	
	/**
	 * This method reads a {@link StringBuffer} from a given {@link File}.
	 * 
	 * @param file
	 *            the {@link File} to be read
	 * 
	 * @return a {@link StringBuffer} with the content of the {@link File}
	 * 
	 * @throws IOException
	 *             if something went wrong with the {@link File}
	 * 
	 * @see #writeFile(File, CharSequence)
	 * @see #read(Reader)
	 * @see FileReader
	 */
	public static StringBuffer readFile( File file ) throws IOException {
		return read( new FileReader( file ) );
	}

	/**
	 * This method reads a {@link StringBuffer} from a given {@link InputStream}.
	 * 
	 * @param inputStream
	 *            the {@link InputStream} to be read
	 * 
	 * @return a {@link StringBuffer} with the content of the {@link InputStream}
	 * 
	 * @throws IOException
	 *             if something went wrong with the {@link InputStream}
	 * 
	 * @see #writeOutputStream(OutputStream, CharSequence)
	 * @see #read(Reader)
	 * @see InputStreamReader
	 */
	public static StringBuffer readInputStream( InputStream inputStream ) throws IOException {
		return read( new InputStreamReader( inputStream ) );
	}
	
	/**
	 * This method reads a {@link StringBuffer} from a given {@link Reader}.
	 * 
	 * @param reader
	 *            the {@link Reader} to be read
	 * 
	 * @return a {@link StringBuffer} with the content of the {@link Reader}
	 * 
	 * @throws IOException
	 *             if something went wrong with the {@link Reader}
	 * 
	 * @see #write(Writer, CharSequence)
	 * @see BufferedReader
	 */
	public static StringBuffer read( Reader reader ) throws IOException {
		BufferedReader r = new BufferedReader( reader, 100000 );
		StringBuffer res = new StringBuffer( 1000000 );
		String help;
		String c="";
		while( ( help = r.readLine() ) != null ) {
			res.append( c + help );
			c="\n";
		}
		r.close();
		return res;
	}

	/**
	 * This method saves a {@link CharSequence} to a given {@link File}.
	 * 
	 * @param outputFile
	 *            the {@link File} into which the output should be written
	 * @param buffer
	 *            the buffer to be written in the {@link File}
	 * 
	 * @throws IOException
	 *             if something went wrong with the {@link File}
	 * 
	 * @see #readFile(File)
	 * @see #write(Writer, CharSequence)
	 */
	public static void writeFile( File outputFile, CharSequence buffer ) throws IOException {
		write( new FileWriter( outputFile ), buffer );
	}
	
	/**
	 * This method saves a {@link CharSequence} to a {@link File} with user-specified file name. If the file already exists and the content differs from that stored in <code>buffer</code>, a BAK-file will be created. 
	 * 
	 * @param fName
	 *            the name of the file to be written
	 * @param buffer
	 *            the buffer to be written in the {@link File}
	 *            
	 * @throws IOException
	 *             if something went wrong with the {@link File}
	 * 
	 * @see #writeFile(String, CharSequence)
	 */
	public static void overwriteFile( String fName, CharSequence buffer ) throws IOException {
		StringBuffer old = null;
		File outputFile = new File( fName );
		if( outputFile.exists() ) {
			old = FileManager.readFile(outputFile);
		}
		if( old == null || !buffer.toString().equals(old.toString()) ) {
			if( old != null ) {
				Path oldPath = Paths.get(outputFile.getAbsolutePath()+".bak");
				Files.delete(oldPath);
				Files.move(Paths.get(fName), oldPath);
			}
			FileManager.writeFile(outputFile, buffer);
		}
	}
	
	
	/**
	 * This method saves a {@link CharSequence} to a {@link File} with user-specified file name.
	 * 
	 * @param fName
	 *            the name of the file to be written
	 * @param buffer
	 *            the buffer to be written in the {@link File}
	 * 
	 * @throws IOException
	 *             if something went wrong with the {@link File}
	 * 
	 * @see #readFile(File)
	 * @see #writeFile(File, CharSequence)
	 */
	public static void writeFile( String fName, CharSequence buffer ) throws IOException {
		writeFile( new File( fName ), buffer );
	}
		
	/**
	 * This method saves a {@link CharSequence} to a given {@link OutputStream}.
	 * 
	 * @param outStream
	 *            the {@link OutputStream} into which the output should be written
	 * @param buffer
	 *            the buffer to be written in the {@link OutputStream}
	 * 
	 * @throws IOException
	 *             if something went wrong with the {@link OutputStream}
	 * 
	 * @see #readInputStream(InputStream)
	 * @see #write(Writer, CharSequence)
	 * @see OutputStreamWriter
	 */
	public static void writeOutputStream( OutputStream outStream, CharSequence buffer ) throws IOException {
		write( new OutputStreamWriter( outStream ), buffer );
	}
	
	/**
	 * This method saves a {@link CharSequence} to a given {@link Writer}.
	 * 
	 * @param writer
	 *            the {@link Writer} into which the output should be written
	 * @param buffer
	 *            the buffer to be written in the {@link Writer}
	 * 
	 * @throws IOException
	 *             if something went wrong with the {@link Writer}
	 * 
	 * @see #read(Reader)
	 * @see BufferedWriter
	 */
	public static void write( Writer writer, CharSequence buffer ) throws IOException {
		BufferedWriter w = new BufferedWriter( writer );
		w.write( buffer.toString() );
		w.close();
	}
}
