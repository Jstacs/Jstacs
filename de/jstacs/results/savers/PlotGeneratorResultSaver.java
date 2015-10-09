package de.jstacs.results.savers;

import java.io.File;

import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.utils.graphics.GraphicsAdaptor;
import de.jstacs.utils.graphics.PDFAdaptor;


public class PlotGeneratorResultSaver implements ResultSaver<PlotGeneratorResult> {

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
	public boolean writeOutput( PlotGeneratorResult result, StringBuffer buf ) {
		throw new RuntimeException( "Impossible for images." );
	}

	
	
}
