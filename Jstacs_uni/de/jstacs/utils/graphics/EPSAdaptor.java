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
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.xmlgraphics.java2d.ps.EPSDocumentGraphics2D;

/**
 * {@link GraphicsAdaptor} for the EPS format.
 * @author Jan Grau
 *
 */
public class EPSAdaptor extends GraphicsAdaptor {

	/**
	 * The EPS document
	 */
	protected EPSDocumentGraphics2D graphics;
	/**
	 * The stream for saving the results
	 */
	protected ByteArrayOutputStream stream;
	
	/**
	 * Creates a new adaptor for plotting to an EPS device.
	 */
	public EPSAdaptor(){
		graphics = new EPSDocumentGraphics2D( true );
		graphics.setGraphicContext(new org.apache.xmlgraphics.java2d.GraphicContext());
		
	}
	
	@Override
	public Graphics2D getGraphics( int width, int height ) throws IOException {
		stream = new ByteArrayOutputStream();		
		graphics.setupDocument( stream, width, height );
		
		return graphics;
	}

	@Override
	public void generateOutput( File file ) throws IOException {
		
		stream.flush();
		stream.close();
		
		FileOutputStream out = new FileOutputStream( file );
		out.write( stream.toByteArray() );
		graphics.finish();
		out.close();
	}

	@Override
	public String getGraphicsExtension() {
		return "eps";
	}
	
	

}
