package de.jstacs.results.savers;

import java.io.File;

import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.utils.graphics.GraphicsAdaptor;
import de.jstacs.utils.graphics.PDFAdaptor;

/**
 * {@link ResultSaver} for {@link PlotGeneratorResult}s.
 * The plots are saved to disk using the {@link GraphicsAdaptor#generateOutput(File)} method.
 * As images are typically not well represented as strings, the {@link ResultSaver#writeOutput(de.jstacs.results.Result, StringBuffer)} method 
 * is marked as {@link Deprecated}.
 * 
 * @author Jan Grau
 *
 */
public class PlotGeneratorResultSaver implements ResultSaver<PlotGeneratorResult> {

	/**
	 * Registers this {@link ResultSaver} in the {@link ResultSaverLibrary}
	 */
	public static void register(){
		ResultSaverLibrary.register( PlotGeneratorResult.class, new PlotGeneratorResultSaver() );
	}

	private PlotGeneratorResultSaver() {
	}

	@Override
	public boolean isAtomic() {
		return true;
	}

	@Override
	public String[] getFileExtensions( PlotGeneratorResult result ) {
		return new String[]{"pdf"};
	}

	@Override
	public boolean writeOutput( PlotGeneratorResult result, File path ) {

		try{
			GraphicsAdaptor ga = new PDFAdaptor();

			result.getValue().generatePlot( ga );

			ga.generateOutput( path );
			
			return true;
		}catch(Exception e){
			e.printStackTrace( );
			return false;
		}

	}

	@Override
	@Deprecated
	public boolean writeOutput( PlotGeneratorResult result, StringBuffer buf ) {
		throw new RuntimeException( "Impossible for images." );
	}

	
	
}
