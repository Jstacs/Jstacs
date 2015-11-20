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
