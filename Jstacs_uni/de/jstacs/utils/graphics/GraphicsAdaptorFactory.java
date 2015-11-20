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
