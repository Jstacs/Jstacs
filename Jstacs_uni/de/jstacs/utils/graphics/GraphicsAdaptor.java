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

package de.jstacs.utils.graphics;


import java.awt.Graphics2D;
import java.io.File;
import java.io.IOException;

/**
 * Generic class for different adaptors for plotting graphics to a file
 * using different graphics formats.
 * 
 * @author Jan Grau
 *
 */
public abstract class GraphicsAdaptor {

	
	/**
	 * Returns a {@link Graphics2D} object for this {@link GraphicsAdaptor} of the given width and height.
	 * This object may be used for plotting using the standard {@link Graphics2D} methods.
	 * @param width the width
	 * @param height the height
	 * @return the {@link Graphics2D} object
	 * @throws IOException if the graphics device could not be created (only applies to some adaptors)
	 */
	public abstract Graphics2D getGraphics(int width, int height) throws IOException;
	
	/**
	 * Generates the outputs and saves the results to a file with the supplied name.
	 * @param filename the file name
	 * @throws IOException if the output could not be written
	 */
	public void generateOutput(String filename) throws IOException{
		generateOutput( new File(filename) );
	}
	
	/**
	 * Generates the outputs and saves the results to the supplied file.
	 * @param file the file
	 * @throws IOException if the output could not be written
	 */
	public abstract void generateOutput(File file) throws IOException;
	
	/**
	 * Returns the file extension for the graphics file format of this {@link GraphicsAdaptor}.
	 * @return the file extension
	 */
	public abstract String getGraphicsExtension();
	
	//TODO createNewGraphics
	
}
