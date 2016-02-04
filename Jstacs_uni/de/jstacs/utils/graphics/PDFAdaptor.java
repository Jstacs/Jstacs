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
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.batik.transcoder.TranscoderException;
import org.apache.batik.transcoder.TranscoderInput;
import org.apache.batik.transcoder.TranscoderOutput;
import org.apache.fop.svg.PDFTranscoder;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * {@link GraphicsAdaptor} for the PDF format.
 * @author Jan Grau
 *
 */
public class PDFAdaptor extends SVGAdaptor {

	private int width;
	private int height;
	
	/**
	 * Creates a new adaptor for plotting to a PDF device.
	 */
	public PDFAdaptor(){
		super();
	}
	
	
	@Override
	public Graphics2D getGraphics( int width, int height ) {
		this.width = width;
		this.height = height;
		return super.getGraphics( width, height );
		
	}

	@Override
	public void generateOutput( File file ) throws IOException {
		
	//	ByteArrayOutputStream pipedOut = new ByteArrayOutputStream();//TODO true pipe
		
	//	graphics.stream( new OutputStreamWriter( pipedOut ), true );
		
	//	pipedOut.close();
		
		Document doc=graphics.getDOMFactory();
		Element rootE=doc.getDocumentElement();
		graphics.getRoot(rootE);
		TranscoderInput input=new TranscoderInput(doc);
		
		PDFTranscoder trans = new PDFTranscoder();
		//Map<TranscodingHints.Key, Object> hints = new HashMap<TranscodingHints.Key, Object>();
		trans.addTranscodingHint( PDFTranscoder.KEY_WIDTH, (float)width );
		trans.addTranscodingHint( PDFTranscoder.KEY_HEIGHT, (float)height );
		trans.addTranscodingHint( PDFTranscoder.KEY_STROKE_TEXT, false );

    	//trans.setTranscodingHints( hints );
		
	  //  ByteArrayInputStream instream = new ByteArrayInputStream( pipedOut.toByteArray() );
	    
	   // TranscoderInput input = new TranscoderInput( document );
	    
	    
	    
	    FileOutputStream stream = new FileOutputStream( file );
	    TranscoderOutput output = new TranscoderOutput( stream );
	    
	    try {
			trans.transcode( input, output );
		} catch ( TranscoderException e ) {
			throw new IOException( e );
		}
	    
	    stream.flush();
	    stream.close();
		
	}

	@Override
	public String getGraphicsExtension() {
		return "pdf";
	}
	
	
	
}
