package de.jstacs.utils.graphics;


import java.awt.Graphics2D;
import java.io.File;
import java.io.IOException;


public abstract class GraphicsAdaptor {

	
	public abstract Graphics2D getGraphics(int width, int height) throws IOException;
	
	public void generateOutput(String filename) throws IOException{
		generateOutput( new File(filename) );
	}
	
	public abstract void generateOutput(File file) throws IOException;
	
	public abstract String getGraphicsExtension();
	
}
