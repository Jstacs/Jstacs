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
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

import org.apache.batik.dom.svg.SVGDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.transcoder.Transcoder;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

/**
 * {@link GraphicsAdaptor} for the SVG format.
 * @author Jan Grau
 *
 */
public class SVGAdaptor extends GraphicsAdaptor {

	/**
	 * The internal graphics object
	 */
	protected SVGGraphics2D graphics;
	/**
	 * The SVG document representation, may be used in sub-classes for different {@link Transcoder}s
	 */
	protected Document document;
	
	/**
	 * Creates a new adaptor for plotting to an SVG device.
	 */
	public SVGAdaptor(){
		DOMImplementation domImpl = SVGDOMImplementation.getDOMImplementation();
		
		String svgNS = SVGDOMImplementation.SVG_NAMESPACE_URI;
	    document = domImpl.createDocument(svgNS, "svg", null);
	    graphics = new SVGGraphics2D(document);
	}
	
	@Override
	public Graphics2D getGraphics( int width, int height ) {
		
		return graphics;
		
	}

	@Override
	public void generateOutput( File file ) throws IOException {
		
		boolean useCSS = true; // we want to use CSS style attributes
	    Writer out = new FileWriter( file );
	    graphics.stream(out, useCSS);
	    out.close();
	}

	@Override
	public String getGraphicsExtension() {
		return "svg";
	}
	
	

}
