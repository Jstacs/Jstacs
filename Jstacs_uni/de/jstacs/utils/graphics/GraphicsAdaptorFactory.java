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


/**
 * Factory class for {@link GraphicsAdaptor}s
 * @author Jan Grau
 *
 */
public class GraphicsAdaptorFactory {

	/**
	 * The allowed output formats
	 * @author dev
	 *
	 */
	public enum OutputFormat{
		/**
		 * SVG format
		 */
		SVG,
		/**
		 * PDF format
		 */
		PDF,
		/**
		 * EPS format
		 */
		EPS,
		/**
		 * PNG format
		 */
		PNG,
		/**
		 * JPEG format
		 */
		JPEG
	}
	
	/**
	 * Returns an appropriat {@link GraphicsAdaptor} for the given format
	 * @param key the format key
	 * @return a new {@link GraphicsAdaptor}
	 */
	public static GraphicsAdaptor getAdaptor(OutputFormat key){
		switch(key){
			case SVG : return new SVGAdaptor();
			case PDF : return new PDFAdaptor();
			case EPS : return new EPSAdaptor();
			case PNG : return new RasterizedAdaptor( "png" );
			case JPEG: return new RasterizedAdaptor( "jpg" );
			default: throw new IllegalArgumentException();
		}
	}
	
	
}
