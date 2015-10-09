package de.jstacs.utils.graphics;


import java.awt.Graphics2D;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Map;

import org.apache.batik.transcoder.TranscoderException;
import org.apache.batik.transcoder.TranscoderInput;
import org.apache.batik.transcoder.TranscoderOutput;
import org.apache.batik.transcoder.TranscodingHints;
import org.apache.fop.svg.PDFTranscoder;
import org.w3c.dom.Document;
import org.w3c.dom.Element;


public class PDFAdaptor extends SVGAdaptor {

	private int width;
	private int height;
	
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
