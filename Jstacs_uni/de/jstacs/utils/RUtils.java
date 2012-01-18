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

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;

import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RFileInputStream;
import org.rosuda.REngine.Rserve.RFileOutputStream;

/**
 * This is a collection of various silly things you can do when accessing
 * <a href="http://www.rforge.net/Rserve/">Rserve</a>.<br>
 * <br>
 * 
 * This class is based on an implementation of Andreas Stephanik and was further
 * developed by Jens Keilwagen.
 * 
 * @author Andreas Stephanik, Jens Keilwagen
 * 
 * @see REnvironment
 */
public class RUtils {

	/**
	 * Installs an R script on the Rserve server.
	 * 
	 * Do not forget to remove the R script by
	 * {@link RConnection#removeFile(String targetname)} at the end of your
	 * session.
	 * 
	 * @param sourcePath
	 *            path of your R script on your local box, this script will be
	 *            copied to the Rserve server (e.g. pdw-24 /tmp/Rserve/conn*)
	 * @param targetName
	 *            just the desired name of the R script on the Rserve server
	 *            (mostly the same as on local box, but without path
	 *            delimiters!), should not contain any path delimiters, since
	 *            Rserve may restrict the access to local working directory
	 * @param rconnection
	 *            an already open connection to the Rserve server
	 * 
	 * @throws Exception
	 *             if something went wrong while installing the R script on the
	 *             Rserve server
	 */
	public static void installRScript( String sourcePath, String targetName, RConnection rconnection ) throws Exception {
		copyFileToServer( sourcePath, targetName, rconnection );
		rconnection.voidEval( "source(\"" + targetName + "\")" );
	}

	/**
	 * Copies a file to the R side. Normally it is send to /tmp/Rserve/conn* .
	 * 
	 * @param source
	 *            the source file
	 * @param targetName
	 *            should not contain any path delimiters, since Rserve may
	 *            restrict the access to local working directory
	 * @param rconnection
	 *            the connection to R
	 * 
	 * @throws Exception
	 *             if something went wrong while copying the file
	 * 
	 */
	public static void copyFileToServer( File source, String targetName, RConnection rconnection ) throws Exception {
		if( !source.exists() ) {
			throw new IOException( "\"" + source.getAbsolutePath() + "\" does not exit." );
		}
		if( !source.isFile() ) {
			throw new IOException( "The source has o be a file." );
		}
		if( !source.canRead() ) {
			throw new IOException( "Can't read the source file." );
		}

		/* create a file on the R side */
		RFileOutputStream os = rconnection.createFile( targetName );

		/* Transfer the r script */
		PrintStream ps = new PrintStream( os );
		FileInputStream fileInputStream = new FileInputStream( source );
		BufferedReader bufferedReader = new BufferedReader( new InputStreamReader( fileInputStream ) );
		String line = bufferedReader.readLine();
		if( line != null ) {
			while( line != null ) {
				ps.println( line );
				line = bufferedReader.readLine();
			}
		}
		bufferedReader.close();
		fileInputStream.close();
		ps.close();
	}

	/**
	 * Copies a file to the R side. Normally it is send to /tmp/Rserve/conn* .
	 * 
	 * @param sourcePath
	 *            the source path
	 * @param targetName
	 *            should not contain any path delimiters, since Rserve may
	 *            restrict the access to local working directory
	 * @param rconnection
	 *            the connection to R
	 * 
	 * @throws Exception
	 *             if something went wrong while copying the file
	 * 
	 * @see RUtils#copyFileToServer(File, String, RConnection)
	 */
	public static void copyFileToServer( String sourcePath, String targetName, RConnection rconnection ) throws Exception {
		copyFileToServer( new File( sourcePath ), targetName, rconnection );
	}

	/**
	 * This method copies a file from the server to the client.
	 * 
	 * @param sourcePath
	 *            the server path name
	 * @param targetPath
	 *            the client path name
	 * @param c
	 *            the connection to R
	 * 
	 * @return the number of copied bytes
	 * 
	 * @throws IOException
	 *             if the file could not be copied
	 */
	public static int copyFileFromServer( String sourcePath, String targetPath, RConnection c ) throws IOException {
		BufferedOutputStream out = new BufferedOutputStream( new FileOutputStream( targetPath ) );
		int code = copyFileFromServer( sourcePath, out, c );
		out.close();
		return code;
	}
	
	/**
	 * Pipes a {@link RFileInputStream} of the given <code>sourcePath</code> into the given {@link OutputStream} <code>out</code>
	 *  
	 * @param sourcePath the server path name
	 * @param out the stream which is used to write the content
	 * @param c the connection to R
	 * 
	 * @return the number of copied bytes
	 * 
	 * @throws IOException if the file could not be copied
	 */
	public static int copyFileFromServer( String sourcePath, OutputStream out, RConnection c ) throws IOException {
		RFileInputStream is = c.openFile( sourcePath );
		int bufSize = 65536, n, code = 0;
		byte[] buf = new byte[bufSize];
		while( ( n = is.read( buf ) ) > -1 ) {
			out.write( buf, 0, n );
			code+=n;
		}
		return code;
	}

	/**
	 * This method opens an {@link RConnection}.
	 * 
	 * @param rServeHostName
	 *            the name of the server with RServe
	 * @param loginName
	 *            the login (if needed)
	 * @param passwd
	 *            the password (if needed)
	 * 
	 * @return a connection to R
	 * 
	 * @throws Exception
	 *             if no connection could be established
	 */
	public static RConnection openRConnection( String rServeHostName, String loginName, String passwd ) throws Exception {
		RConnection c = new RConnection( rServeHostName );
		// if server requires authentication, send one
		if( c.needLogin() ) {
			c.login( loginName, passwd );
		}
		return c;
	}
}
