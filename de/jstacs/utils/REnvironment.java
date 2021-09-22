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

import java.awt.Canvas;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.MediaTracker;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.image.BufferedImage;
import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.Iterator;

import javax.imageio.ImageIO;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.rosuda.REngine.REXP;
import org.rosuda.REngine.RList;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;

import de.jstacs.io.RegExFilenameFilter;
import de.jstacs.io.RegExFilenameFilter.Directory;

/**
 * This is an environment class for <a href="http://www.r-project.org/">R</a> and <a href="http://www.rforge.net/Rserve/">Rserve</a> that helps you to handle all things in a
 * more gentle way.
 * 
 * @author Jens Keilwagen
 * 
 * @see RUtils
 */
public class REnvironment {

	private RConnection c;

	private HashSet<String> filesOnTheServer;

	private int initSize;

	private boolean windows, pngOkay;

	/**
	 * Creates a new {@link REnvironment}.
	 * 
	 * @param server
	 *            the name of the server
	 * @param defaultPort
	 *            if <code>true</code> the default port is used, otherwise the
	 *            given port
	 * @param port
	 *            the given port, that will be used if
	 *            <code>defaultPort=false</code>
	 * @param user
	 *            the login
	 * @param passwd
	 *            the password for the login
	 * @param initSize
	 *            the initial number of slots for files on the server (has to be
	 *            at least 1)
	 * @param pngOkay
	 *            if png can be used for creating BufferedImages  			
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	private REnvironment( String server, boolean defaultPort, int port, String user, String passwd, int initSize, boolean pngOkay ) throws Exception {
		if( defaultPort ) {
			c = new RConnection( server );
		} else {
			c = new RConnection( server, port );
		}
		if( c.needLogin() ) { // if server requires authentication, send one
			c.login( user, passwd );
		}
		if( initSize < 1 ) {
			throw new IOException( "initSize has to be at least 1" );
		}
		this.initSize = initSize;
		resetFilesOnTheServer();
		windows = eval( ".Platform$OS.type" ).asString().equalsIgnoreCase( "windows" );
		this.pngOkay = pngOkay;
	}
	
	/**
	 * This constructor creates a connection to an R instance on the &quot;localhost&quot; with empty login and password.
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see REnvironment#REnvironment(String, String, String)
	 */
	public REnvironment() throws Exception {
		this( "localhost", "", "" );
	}

	/**
	 * Creates a new {@link REnvironment} using the default port.
	 * 
	 * @param server
	 *            the name of the server
	 * @param user
	 *            the login
	 * @param passwd
	 *            the password for the login
	 * @param initSize
	 *            the initial number of slots for files on the server (has to be
	 *            at least 1)
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public REnvironment( String server, String user, String passwd, int initSize ) throws Exception {
		this( server, true, -1, user, passwd, initSize, true );
	}

	/**
	 * Creates a new {@link REnvironment} with initial size 10 using the default
	 * port.
	 * 
	 * @param server
	 *            the name of the server
	 * @param user
	 *            the login
	 * @param passwd
	 *            the password for the login
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see REnvironment#REnvironment(String, String, String, int)
	 */
	public REnvironment( String server, String user, String passwd ) throws Exception {
		this( server, true, -1, user, passwd, 10, true );
	}

	/**
	 * Creates a new {@link REnvironment} using the given port.
	 * 
	 * @param server
	 *            the name of the server
	 * @param port
	 *            the port
	 * @param user
	 *            the login
	 * @param passwd
	 *            the password for the login
	 * @param initSize
	 *            the initial number of slots for files on the server (has to be
	 *            at least 1)
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public REnvironment( String server, int port, String user, String passwd, int initSize ) throws Exception {
		this( server, false, port, user, passwd, initSize, true );
	}

	/**
	 * Creates a new {@link REnvironment} with initial size 10 using the given
	 * port.
	 * 
	 * @param server
	 *            the name of the server
	 * @param port
	 *            the port
	 * @param user
	 *            the login
	 * @param passwd
	 *            the password for the login
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see REnvironment#REnvironment(String, String, String, int)
	 */
	public REnvironment( String server, int port, String user, String passwd ) throws Exception {
		this( server, false, port, user, passwd, 10, true );
	}

	/**
	 * Copies a file from the server to the local machine.
	 * 
	 * @param serverFileName
	 *            name of the file on the server
	 * @param clientFileName
	 *            name of the file on the client (local machine)
	 * @param overwriteExistingFile
	 *            whether to overwrite the file, if already existing on the
	 *            client
	 * 
	 * @return <code>true</code> if the file is copied successfully
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public boolean copyFileFromServer( String serverFileName, String clientFileName, boolean overwriteExistingFile ) throws Exception {
		if( !( new File( clientFileName ) ).exists() || overwriteExistingFile ) {
			RUtils.copyFileFromServer( serverFileName, clientFileName, c );
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Copies a file from the local machine to the server.
	 * 
	 * @param serverFileName
	 *            name of the file on the server
	 * @param clientFileName
	 *            name of the file on the client (local machine)
	 * @param overwriteExistingFile
	 *            whether to overwrite the file, if already existing on the
	 *            server
	 * 
	 * @return <code>true</code> if the file is copied successfully
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public boolean copyFileToServer( String clientFileName, String serverFileName, boolean overwriteExistingFile ) throws Exception {
		if( filesOnTheServer.contains( serverFileName ) ) {
			if( overwriteExistingFile ) {
				RUtils.copyFileToServer( clientFileName, serverFileName, c );
				return true;
			} else {
				return false;
			}
		} else {
			RUtils.copyFileToServer( clientFileName, serverFileName, c );
			filesOnTheServer.add( serverFileName );
			return true;
		}
	}

	/**
	 * Closes the {@link REnvironment} and removes all files from the server.
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public void close() throws Exception {
		deleteAllFilesAtTheServer();
		c.close();
	}

	/**
	 * Creates a matrix of <code>int</code>s.
	 * 
	 * @param matrixName
	 *            the name of the matrix
	 * @param matrix
	 *            the matrix of <code>int</code>s
	 * 
	 * @return an R expression
	 * 
	 * @throws RserveException
	 *             if something with Rserve went wrong
	 * @throws IllegalArgumentException
	 *             if <code>matrix[i].length != matrix[j].length</code>
	 */
	public REXP createMatrix( String matrixName, int[][] matrix ) throws RserveException, IllegalArgumentException {
		if( matrix == null ) {
			return c.eval( matrixName + " = NULL;" );
		} else {		
			StringBuffer cmd = new StringBuffer( matrix.length * matrix[0].length * 20 );
			cmd.append( matrixName + " = matrix( c(" );
			int rows = matrix.length, columns = matrix[0].length, counter1 = 0, counter2;
			while( counter1 < rows ) {
				counter2 = 0;
				if( columns != matrix[0].length ) {
					throw new IllegalArgumentException( "The matrix is not rectangular" );
				}
				while( counter2 < columns ) {
					if( counter1 == 0 && counter2 == 0 ) {
						cmd.append( matrix[counter1][counter2] );
					} else {
						cmd.append( ", " + matrix[counter1][counter2] );
					}
					counter2++;
				}
				counter1++;
			}
			return c.eval( cmd + "), nrow = " + rows + ", ncol = " + columns + ", byrow = TRUE );" );
		}
	}

	/**
	 * Creates a matrix of <code>double</code>s.
	 * 
	 * @param matrixName
	 *            the name of the matrix
	 * @param matrix
	 *            the matrix of <code>double</code>s
	 * 
	 * @return an R expression
	 * 
	 * @throws RserveException
	 *             if something with Rserve went wrong
	 * @throws IllegalArgumentException
	 *             if <code>matrix[i].length != matrix[j].length</code>
	 */
	public REXP createMatrix( String matrixName, double[][] matrix ) throws RserveException, IllegalArgumentException {
		if( matrix == null ) {
			return c.eval( matrixName + " = NULL;" );
		} else {		
			StringBuffer cmd = new StringBuffer( matrix.length * matrix[0].length * 20 );
			cmd.append( matrixName + " = matrix( c(" );
			int rows = matrix.length, columns = matrix[0].length, counter1 = 0, counter2;
			while( counter1 < rows ) {
				counter2 = 0;
				if( columns != matrix[0].length ) {
					throw new IllegalArgumentException( "The matrix is not rectangular" );
				}
				while( counter2 < columns ) {
					if( counter1 == 0 && counter2 == 0 ) {
						cmd.append( getDoubleVal( matrix[counter1][counter2] ) );
					} else {
						cmd.append( ", " + getDoubleVal( matrix[counter1][counter2] ) );
					}
					counter2++;
				}
				counter1++;
			}
			cmd.append( "), nrow = " + rows + ", ncol = " + columns + ", byrow = TRUE );" );
			return c.eval( cmd.toString() );
		}
	}

	private String getDoubleVal( double val ) {
		if( Double.isInfinite( val ) ) {
			if( Double.NEGATIVE_INFINITY == val ) {
				return "-Inf";
			} else {
				return "Inf";
			}
		} else {
			return "" + val;
		}
	}

	/**
	 * Creates a vector of {@link String}s.
	 * 
	 * @param vectorName
	 *            the name of the vector
	 * @param vector
	 *            the vector of {@link String}s
	 * 
	 * @return an R expression
	 * 
	 * @throws RserveException
	 *             if something with Rserve went wrong
	 */
	public REXP createVector( String vectorName, String[] vector ) throws RserveException {
		if( vector == null ) {
			return c.eval( vectorName + " = NULL;" );
		} else {
			StringBuffer cmd = new StringBuffer( vector.length * 100 );
			cmd.append( vectorName + " = c(" );
			if( vector != null && vector.length > 0 ) {
				cmd.append( "\"" + vector[0] + "\"" );
				for( int i = 1; i < vector.length; i++ ) {
					cmd.append( ", \"" + vector[i] + "\"" );
				}
			}
			return c.eval( cmd + ");" );
		}
	}

	/**
	 * Creates a vector of <code>int</code>s.
	 * 
	 * @param vectorName
	 *            the name of the vector
	 * @param vector
	 *            the vector of <code>int</code>s
	 * 
	 * @return an R expression
	 * 
	 * @throws RserveException
	 *             if something with Rserve went wrong
	 */
	public REXP createVector( String vectorName, int[] vector ) throws RserveException {
		if( vector == null ) {
			return c.eval( vectorName + " = NULL;" );
		} else {
			StringBuffer cmd = new StringBuffer( vector.length * 20 );
			cmd.append( vectorName + " = c(" );
			if( vector != null && vector.length > 0 ) {
				cmd.append( vector[0] );
				for( int i = 1; i < vector.length; i++ ) {
					cmd.append( ", " + vector[i] );
				}
			}
			return c.eval( cmd + ");" );
		}
	}

	/**
	 * Creates a vector of <code>long</code>s.
	 * 
	 * @param vectorName
	 *            the name of the vector
	 * @param vector
	 *            the vector of <code>long</code>s
	 * 
	 * @return an R expression
	 * 
	 * @throws RserveException
	 *             if something with Rserve went wrong
	 */
	public REXP createVector( String vectorName, long[] vector ) throws RserveException {
		if( vector == null ) {
			return c.eval( vectorName + " = NULL;" );
		} else {
			StringBuffer cmd = new StringBuffer( vector.length * 20 );
			cmd.append( vectorName + " = c(" );
			if( vector != null && vector.length > 0 ) {
				cmd.append( vector[0] );
				for( int i = 1; i < vector.length; i++ ) {
					cmd.append( ", " + vector[i] );
				}
			}
			return c.eval( cmd + ");" );
		}
	}

	/**
	 * Creates a vector of <code>double</code>s.
	 * 
	 * @param vectorName
	 *            the name of the vector
	 * @param vector
	 *            the vector of <code>double</code>s
	 * 
	 * @return an R expression
	 * 
	 * @throws RserveException
	 *             if something with Rserve went wrong
	 */
	public REXP createVector( String vectorName, double[] vector ) throws RserveException {
		if( vector == null ) {
			return c.eval( vectorName + " = NULL;" );
		} else {
			StringBuffer cmd = new StringBuffer( vector.length * 20 );
			cmd.append( vectorName + " = c(" );
			if( vector != null && vector.length > 0 ) {
				cmd.append( getDoubleVal( vector[0] ) );
				for( int i = 1; i < vector.length; i++ ) {
					cmd.append( ", " + getDoubleVal( vector[i] ) );
				}
			}
			return c.eval( cmd + ");" );
		}
	}

	/**
	 * Deletes all files that have been copied to the server or created on the
	 * server.
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public void deleteAllFilesAtTheServer() throws Exception {
		Iterator<String> it = filesOnTheServer.iterator();
		Exception e = null;
		int removed = 0, resisted = 0;
		while( it.hasNext() ) {
			try {
				c.removeFile( it.next() );
				removed++;
			}catch( Exception ex ) {
				e = ex;
				resisted++;
			}
		}
		resetFilesOnTheServer();
		if( e != null ) {
			System.out.println( resisted + " files resisted (" + removed + " files could be removed)" );
			throw e;
		}
	}

	/**
	 * Evaluates the {@link String} as R commands.
	 * 
	 * @param cmd
	 *            the {@link String} to be evaluated
	 * 
	 * @return the result from R
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public REXP eval( CharSequence cmd ) throws Exception {
		return c.eval( cmd.toString() );
	}

	/**
	 * Returns information about the version of R that is used.
	 * 
	 * @return information about the version of R that is used
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public String getVersionInformation() throws Exception {
		RList l = eval( "version" ).asList();
		String[] s = l.keys();
		String erg = "";
		for( int counter1 = 0; counter1 < s.length; counter1++ ) {
			erg += s[counter1] + ": \"" + l.at( s[counter1] ).asString() + "\"\n";
		}
		return erg;
	}

	/**
	 * Installs a script on the server.
	 * 
	 * @param clientFileName
	 *            the name of the scriptfile on the client
	 * @param serverFileName
	 *            the name of the scriptfile on the server
	 * @param overwriteExistingFile
	 *            if <code>true</code> the method is enabled to overwrite an
	 *            existing file
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public void installScript( String clientFileName, String serverFileName, boolean overwriteExistingFile ) throws Exception {
		copyFileToServer( clientFileName, serverFileName, overwriteExistingFile );
		voidEval( "source(\"" + serverFileName + "\")" );
	}

	/**
	 * Creates a buffered image from a given plot command.
	 * 
	 * <br>
	 * <br>
	 * 
	 * If you use a java version below 1.6 your image is internally encoded as
	 * jpeg.<br>
	 * If you use a java version at least 1.6 your image is internally encoded
	 * as png.
	 * 
	 * @param pltcmd
	 *            the plot command
	 * 
	 * @return the buffered image
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see REnvironment#plot(CharSequence, double, double)
	 * @see REnvironment#showImage(String, BufferedImage)
	 * @see ImageIO#write(java.awt.image.RenderedImage, String, File)
	 */
	public BufferedImage plot( CharSequence pltcmd ) throws Exception {
		return plot( pltcmd, -1, -1 );
	}

	/**
	 * Creates a buffered image with given dimension from a given plot command.
	 * 
	 * <br>
	 * <br>
	 * 
	 * If you use a java version below 1.6 your image is internally encoded as
	 * jpeg.<br>
	 * If you use a java version at least 1.6 your image is internally encoded
	 * as png.
	 * 
	 * @param pltcmd
	 *            the plot command
	 * @param width
	 *            the width of the image (in pixel)
	 * @param height
	 *            the height of the image (in pixel)
	 * 
	 * @return the buffered image
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see REnvironment#showImage(String, BufferedImage)
	 * @see ImageIO#write(java.awt.image.RenderedImage, java.lang.String,
	 *      java.io.File)
	 */
	public BufferedImage plot( CharSequence pltcmd, double width, double height ) throws Exception {
		String serverFileName = "help-" + System.currentTimeMillis();
		if( windows ) {
			if( pngOkay ) {

				if( height > 0 && width > 0 ) {
					plot( eval( "try( png( filename = \"" + serverFileName + "\", width = " + width + ", height = " + height + " ) )" ),
							pltcmd,
							"png" );
				} else {
					plot( eval( "try( png( filename = \"" + serverFileName + "\" ) )" ), pltcmd, "png" );
				}
			} else {
				if( height > 0 && width > 0 ) {
					plot( eval( "try( jpeg( filename = \"" + serverFileName
								+ "\", width = "
								+ width
								+ ", height = "
								+ height
								+ ", quality = 100) )" ), pltcmd, "jpeg" );
				} else {
					plot( eval( "try( jpeg( filename = \"" + serverFileName + "\" ), quality = 100 )" ), pltcmd, "jpeg" );
				}
			}
		} else {
			String file, type;
			if( pngOkay ) {
				file = "png";
				type = "png16m";
			} else {
				file = "jpeg";
				type = "jpeg";
			}

			if( height > 0 && width > 0 ) {
				plot( eval( "try( bitmap( file = \"" + serverFileName
							+ "\", width = "
							+ width
							/ 72
							+ ", height = "
							+ height
							/ 72
							+ ", type=\""
							+ type
							+ "\", taa=4, gaa=4, pointsize=12 ) )" ), pltcmd, file );
			} else {
				plot( eval( "try( bitmap( file = \"" + serverFileName + "\", type=\"" + type + "\", taa=4, gaa=4, pointsize=12 ) )" ), pltcmd, file );
			}
		}
		ByteArrayOutputStream bais = new ByteArrayOutputStream();
		BufferedOutputStream out = new BufferedOutputStream( bais );
		RUtils.copyFileFromServer(serverFileName, out, c);
		out.flush();
		ByteArrayInputStream in = new ByteArrayInputStream( bais.toByteArray() );
		BufferedImage i = ImageIO.read( in );
		out.close();
		in.close();
		filesOnTheServer.add( serverFileName );
		return i;
	}
	
	/**
	 * This method creates an image with given dimension from a given plot command and pipes it to the output stream.
	 * 
	 * @param pltcmd the plotting command
	 * @param width the width of the image
	 * @param height the height of the image
	 * @param formatName the name of the image format used for writing the image (see {@link ImageIO})
	 * @param out the {@link OutputStream} used for writing the image
	 * 
	 * @throws IOException see {@link ImageIO#write(java.awt.image.RenderedImage, String, OutputStream)}
	 * @throws Exception see {@link #plot(CharSequence, double, double)}
	 * 
	 * @see ImageIO#write(java.awt.image.RenderedImage, String, OutputStream)
	 */
	public void plot( CharSequence pltcmd, double width, double height, String formatName, OutputStream out ) throws IOException, Exception {
		ImageIO.write( plot( pltcmd, width, height ), formatName, out );
	}

	/**
	 * Creates a pdf file from a given plot command.
	 * 
	 * @param pltcmd
	 *            the plot command
	 * @param fileName
	 *            the name of the pdf file
	 * @param overwriteExistingFile
	 *            if <code>true</code> the method is enabled to overwrite an
	 *            existing file
	 * 
	 * @return <code>true</code> if the plot is done
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public boolean plotToPDF( CharSequence pltcmd, String fileName, boolean overwriteExistingFile ) throws Exception {
		return plotToPDF( pltcmd, -1, -1, fileName, overwriteExistingFile );
	}

	/**
	 * Creates a pdf file with given dimension from a given plot command.
	 * 
	 * @param pltcmd
	 *            the plot command
	 * @param width
	 *            the width of the image
	 * @param height
	 *            the height of the image
	 * @param fileName
	 *            the name of the pdf file
	 * @param overwriteExistingFile
	 *            if <code>true</code> the method is enabled to overwrite an
	 *            existing file
	 * 
	 * @return <code>true</code> if the plot is done
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public boolean plotToPDF( CharSequence pltcmd, double width, double height, String fileName, boolean overwriteExistingFile ) throws Exception {
		String graphicTyp = "pdf", serverFileName = "tmp-" + System.currentTimeMillis() + "." + graphicTyp;
		if( height > 0 && width > 0 ) {
			return plotToFile( eval( "try( pdf( \"" + serverFileName + "\", width = " + width + ", height = " + height + " ) )" ),
					pltcmd,
					fileName,
					serverFileName,
					graphicTyp,
					overwriteExistingFile );
		} else {
			return plotToFile( eval( "try( pdf( \"" + serverFileName + "\" ) )" ),
					pltcmd,
					fileName,
					serverFileName,
					graphicTyp,
					overwriteExistingFile );
		}
	}

	/**
	 * Creates a tex file from a given plot command.
	 * 
	 * @param pltcmd
	 *            the plot command
	 * @param fileName
	 *            the name of the tex file
	 * @param overwriteExistingFile
	 *            if <code>true</code> the method is enabled to overwrite an
	 *            existing file
	 * 
	 * @return <code>true</code> if the plot is done
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public boolean plotToTexFile( CharSequence pltcmd, String fileName, boolean overwriteExistingFile ) throws Exception {
		return plotToTexFile( pltcmd, -1, -1, fileName, overwriteExistingFile );
	}

	/**
	 * Creates a tex file with given dimension from a given plot command.
	 * 
	 * @param pltcmd
	 *            the plot command
	 * @param width
	 *            the width of the image
	 * @param height
	 *            the height of the image
	 * @param fileName
	 *            the name of the tex file
	 * @param overwriteExistingFile
	 *            if <code>true</code> the method is enabled to overwrite an
	 *            existing file
	 * 
	 * @return <code>true</code> if the plot is done
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public boolean plotToTexFile( CharSequence pltcmd, double width, double height, String fileName, boolean overwriteExistingFile ) throws Exception {
		String graphicTyp = "tex", serverFileName = "tmp-" + System.currentTimeMillis() + "." + graphicTyp;
		if( height > 0 && width > 0 ) {
			return plotToFile( eval( "try( pictex( \"" + serverFileName + "\", width = " + width + ", height = " + height + " ) )" ),
					pltcmd,
					fileName,
					serverFileName,
					graphicTyp,
					overwriteExistingFile );
		} else {
			return plotToFile( eval( "try( pictex( \"" + serverFileName + "\" ) )" ),
					pltcmd,
					fileName,
					serverFileName,
					graphicTyp,
					overwriteExistingFile );
		}
	}

	/**
	 * Evaluates the {@link String} as R commands.
	 * 
	 * @param cmd
	 *            the {@link String} to be evaluated
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public void voidEval( CharSequence cmd ) throws Exception {
		c.voidEval( cmd.toString() );
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#finalize()
	 */
	@Override
	protected void finalize() throws Throwable {
		close();
		super.finalize();
	}

	private boolean plotToFile( REXP xp, CharSequence pltcmd, String fileName, String serverFileName, String graphicTyp,
			boolean overwriteExistingFile ) throws Exception {
		filesOnTheServer.add( serverFileName );
		plot( xp, pltcmd, graphicTyp );
		boolean b = copyFileFromServer( serverFileName, fileName, overwriteExistingFile );
		return b;
	}

	private void plot( REXP xp, CharSequence pltcmd, String graphicTyp ) throws Exception {
		if( !xp.isNull() ) {
			// if there's a string then we have a problem, R sent an error
			String msg = "Can't open " + graphicTyp + " graphics device:\n\t" + xp.asString();
			xp = eval( "if (exists(\"last.warning\") && length(last.warning)>0) names(last.warning)[1] else 0" );
			if( xp.asString() != null ) {
				msg += "\t=> " + xp.asString();
			}
			throw new IOException( msg );
		}
		// ok, so the device should be fine - let's plot
		voidEval( pltcmd );
		voidEval( "dev.off()" );
	}

	private void resetFilesOnTheServer() {
		if( filesOnTheServer == null ) {
			filesOnTheServer = new HashSet<String>( initSize );
		} else {
			filesOnTheServer.clear();
		}
	}

	/**
	 * Enables you to show an image.
	 * 
	 * @param title
	 *            the title of the frame
	 * @param img
	 *            the image
	 * 
	 * @return a frame showing the image
	 * 
	 * @throws InterruptedException
	 *             if the process is interrupted
	 * 
	 * @see REnvironment#showImage(String, BufferedImage, int)
	 */
	public static JFrame showImage( String title, BufferedImage img ) throws InterruptedException {
		return showImage( title, img, JFrame.DISPOSE_ON_CLOSE );
	}

	/**
	 * Enables you to show an image.
	 * 
	 * @param title
	 *            the title of the frame
	 * @param img
	 *            the image
	 * @param defaultCloseOperation
	 *            the variable used to control the window-closing operation
	 * 
	 * @return a frame showing the image
	 * 
	 * @throws InterruptedException
	 *             if the process is interrupted
	 * 
	 * @see javax.swing.WindowConstants
	 */
	public static JFrame showImage( String title, BufferedImage img, int defaultCloseOperation ) throws InterruptedException {
		JFrame f = new JFrame( title );
		f.add( new ImageDisplay( img, true ) );
		f.setDefaultCloseOperation( defaultCloseOperation );
		f.pack();
		f.setVisible( true );
		return f;
	}

	/**
	 * The class for displaying an image.
	 */
	private static class ImageDisplay extends Canvas {

		private static final long serialVersionUID = 1L;

		private boolean rescale;

		private BufferedImage img;

		/**
		 * Constructor for an {@link ImageDisplay}.
		 * 
		 * @param img
		 *            the image to be displayed
		 * @param rescale
		 *            a switch whether to allow rescaling
		 * 
		 * @throws InterruptedException
		 *             if the process is interrupted
		 */
		public ImageDisplay( BufferedImage img, boolean rescale ) throws InterruptedException {
			this.img = img;
			MediaTracker mediaTracker = new MediaTracker( this );
			mediaTracker.addImage( img, 0 );
			mediaTracker.waitForID( 0 );
			setSize( img.getWidth( this ), img.getHeight( this ) );
			setRescale( rescale );
			addMouseListener( new MyMouseListener() );
		}

		/**
		 * Constructor for an {@link ImageDisplay}.
		 * 
		 * @param img
		 *            the image to be displayed
		 * @param rescale
		 *            a switch whether to allow rescaling
		 * @param width
		 *            the initial frame width
		 * @param height
		 *            the initial frame height
		 * 
		 * @throws InterruptedException
		 *             if the process is interrupted
		 */
		public ImageDisplay( BufferedImage img, boolean rescale, int width, int height ) throws InterruptedException {
			this.img = img;
			MediaTracker mediaTracker = new MediaTracker( this );
			mediaTracker.addImage( img, 0 );
			mediaTracker.waitForID( 0 );
			setSize( width, height );
			setRescale( rescale );
		}

		private class MyMouseListener implements MouseListener {

			/* (non-Javadoc)
			 * @see java.awt.event.MouseListener#mouseClicked(java.awt.event.MouseEvent)
			 */
			public void mouseClicked( MouseEvent e ) {
				JFileChooser chooser = new JFileChooser();
				chooser.setAcceptAllFileFilterUsed( false );
				chooser.setFileFilter( new RegExFilenameFilter( "png-image", Directory.ALLOWED, false, ".*\\.png" ) );
				chooser.setMultiSelectionEnabled( false );
				if( chooser.showSaveDialog( new JFrame() ) == JFileChooser.APPROVE_OPTION ) {
					try {
						File f = chooser.getSelectedFile();
						String name = f.getAbsolutePath();
						if( !name.endsWith( ".png" ) ) {
							f = new File( name + ".png" );
						}
						ImageIO.write( img, "png", f );
					} catch ( Exception ex ) {
						JOptionPane.showMessageDialog( new JFrame(), "Could not save the image:" + ex.getClass().getName() );
					}
				}
			}

			/* (non-Javadoc)
			 * @see java.awt.event.MouseListener#mousePressed(java.awt.event.MouseEvent)
			 */
			public void mousePressed( MouseEvent e ) {}

			/* (non-Javadoc)
			 * @see java.awt.event.MouseListener#mouseReleased(java.awt.event.MouseEvent)
			 */
			public void mouseReleased( MouseEvent e ) {}

			/* (non-Javadoc)
			 * @see java.awt.event.MouseListener#mouseEntered(java.awt.event.MouseEvent)
			 */
			public void mouseEntered( MouseEvent e ) {}

			/* (non-Javadoc)
			 * @see java.awt.event.MouseListener#mouseExited(java.awt.event.MouseEvent)
			 */
			public void mouseExited( MouseEvent e ) {}
		}

		private void setRescale( boolean rescale ) {
			this.rescale = rescale;
		}

		/* (non-Javadoc)
		 * @see java.awt.Canvas#paint(java.awt.Graphics)
		 */
		@Override
		public void paint( Graphics g ) {
			if( rescale ) {
				g.drawImage( img.getScaledInstance( getWidth(), getHeight(), Image.SCALE_SMOOTH ), 0, 0, null );
			} else {
				g.drawImage( img, 0, 0, null );
			}
		}
	}
}