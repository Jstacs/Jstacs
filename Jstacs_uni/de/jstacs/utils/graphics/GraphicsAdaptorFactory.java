package de.jstacs.utils.graphics;



public class GraphicsAdaptorFactory {

	public enum OutputFormat{
		SVG,
		PDF,
		EPS,
		PNG,
		JPEG
	}
	
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
