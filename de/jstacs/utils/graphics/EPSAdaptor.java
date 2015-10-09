package de.jstacs.utils.graphics;


import java.awt.Graphics2D;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.xmlgraphics.java2d.ps.EPSDocumentGraphics2D;


public class EPSAdaptor extends GraphicsAdaptor {

	protected EPSDocumentGraphics2D graphics;
	protected ByteArrayOutputStream stream;
	
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
