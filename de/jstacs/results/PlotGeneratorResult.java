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

package de.jstacs.results;

import java.awt.Graphics2D;

import de.jstacs.DataType;
import de.jstacs.Storable;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.savers.PlotGeneratorResultSaver;
import de.jstacs.utils.graphics.GraphicsAdaptor;

/**
 * Class for a {@link Result} that may be used to generate plots for different output formats using
 * {@link GraphicsAdaptor} sub-classes. Graphics are generated on demand.
 * 
 * @author Jan Grau
 *
 */
public class PlotGeneratorResult extends Result {

	/**
	 * Interface for a class that may generate a plot using the specified {@link GraphicsAdaptor}.
	 * @author Jan Grau
	 *
	 */
	public static interface PlotGenerator extends Storable{
		
		/**
		 * Generates the plot using the {@link Graphics2D} device of the supplied {@link GraphicsAdaptor}.
		 * @param ga the graphics adaptor
		 * @throws Exception if the plot could not be generated
		 */
		public void generatePlot(GraphicsAdaptor ga) throws Exception;
		
	}
	
	
	private PlotGenerator gen;
	private boolean isStatic;
	private PlotGeneratorResultSaver.Format outputFormat;
	
	/**
	 * Creates a new {@link PlotGeneratorResult} with the given name, comment, {@link PlotGenerator}.
	 * @param name the name of the result
	 * @param comment a comment on the result
	 * @param gen the object that may generate the plot
	 * @param isStatic if <code>true</code>, the plot is considered static and may be cached.
	 */
	public PlotGeneratorResult( String name, String comment, PlotGenerator gen, boolean isStatic ) {
		this(name, comment, gen, isStatic, PlotGeneratorResultSaver.Format.PDF);
	}
	
	
	/**
	 * Creates a new {@link PlotGeneratorResult} with the given name, comment, {@link PlotGenerator}.
	 * @param name the name of the result
	 * @param comment a comment on the result
	 * @param gen the object that may generate the plot
	 * @param isStatic if <code>true</code>, the plot is considered static and may be cached.
	 * @param outputFormat the output format
	 */
	public PlotGeneratorResult( String name, String comment, PlotGenerator gen, boolean isStatic, PlotGeneratorResultSaver.Format outputFormat ) {
		super( name, comment, DataType.IMAGE );
		this.gen = gen;
		this.isStatic = isStatic;
		this.outputFormat = outputFormat;
	}

	/**
	 * Creates a new {@link PlotGeneratorResult} from its XML representation.
	 * @param rep the XML representation
	 * @throws NonParsableException if XML could not be parsed
	 */
	public PlotGeneratorResult( StringBuffer rep ) throws NonParsableException {
		super( rep );
	}

	@Override
	public String getXMLTag() {
		return "PlotGeneratorResult";
	}

	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		XMLParser.appendObjectWithTags( buf, gen, "generator" );
		XMLParser.appendObjectWithTags( buf, isStatic, "isStatic" );
		XMLParser.appendObjectWithTags(buf, this.outputFormat.toString(), "outputFormat");
	}

	@Override
	protected void extractFurtherInfos( StringBuffer buf ) throws NonParsableException {
		gen = (PlotGenerator)XMLParser.extractObjectForTags( buf, "generator" );
		isStatic = (Boolean)XMLParser.extractObjectForTags( buf, "isStatic" );
		try {
            String outform = (String)XMLParser.extractObjectForTags(buf, "outputFormat");
            this.outputFormat = PlotGeneratorResultSaver.Format.valueOf(outform);
        }
        catch (NonParsableException ex) {
            this.outputFormat = PlotGeneratorResultSaver.Format.PDF;
        }
	}

	@Override
	public PlotGenerator getValue() {
		return gen;
	}

	/**
	 * Returns <code>true</code> if the plot is considered static and may be cached.
	 * @return if the plot is considered static and may be cached.
	 */
	public boolean isStatic() {
		return isStatic;
	}
	
	public PlotGeneratorResultSaver.Format getFormat() {
        return this.outputFormat;
    }

}
