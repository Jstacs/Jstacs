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

package de.jstacs.results.savers;

import java.io.File;

import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.utils.graphics.GraphicsAdaptor;
import de.jstacs.utils.graphics.PDFAdaptor;
import de.jstacs.utils.graphics.RasterizedAdaptor;

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

	
	public enum Format{
		PDF,
		PNG
	}
	
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
		if (result.getFormat() == PlotGeneratorResultSaver.Format.PDF) {
            return new String[] { "pdf" };
        }else{
        	return new String[] { "png" };
        }
	}

	@Override
	public boolean writeOutput( PlotGeneratorResult result, File path ) {

		GraphicsAdaptor ga = null;
        if (result.getFormat() == PlotGeneratorResultSaver.Format.PDF) {
            ga = (GraphicsAdaptor)new PDFAdaptor();
        }
        else {
            ga = (GraphicsAdaptor)new RasterizedAdaptor("png");
        }
        
		try{

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
