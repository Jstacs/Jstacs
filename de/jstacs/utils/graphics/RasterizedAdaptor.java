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
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

/**
 * {@link GraphicsAdaptor} for rasterized formats, namely PNG or JPEG.
 * @author Jan Grau
 *
 */
public class RasterizedAdaptor extends GraphicsAdaptor {

	private BufferedImage img;
	/**
	 * The graphics object used for plotting
	 */
	protected Graphics2D graphics;
	private String type;
	
	/**
	 * Creates a new {@link RasterizedAdaptor} for different formats.
	 * Currently, PNG and JPEG are supported
	 * @param type "png" or "jpg"
	 */
	public RasterizedAdaptor(String type){
		this.type = type;
	}
	
	@Override
	public Graphics2D getGraphics( int width, int height ) {
		img = new BufferedImage( width, height, BufferedImage.TYPE_INT_RGB);
		graphics = (Graphics2D)img.getGraphics();
		
		graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		graphics.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		
		return graphics;
	}

	@Override
	public void generateOutput( File file ) throws IOException {
		
		ImageIO.write( img, type, file );

	}

	@Override
	public String getGraphicsExtension() {
		return type;
	}
	
	/**
	 * Returns the internal image as a {@link BufferedImage}
	 * @return the image
	 */
	public BufferedImage getImage(){
		return img;
	}
	

}
