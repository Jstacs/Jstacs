package de.jstacs.results;

import de.jstacs.DataType;
import de.jstacs.Storable;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.graphics.GraphicsAdaptor;


public class PlotGeneratorResult extends Result {

	public static interface PlotGenerator extends Storable{
		
		public void generatePlot(GraphicsAdaptor ga) throws Exception;
		
	}
	
	
	private PlotGenerator gen;
	private boolean isStatic;
	
	public PlotGeneratorResult( String name, String comment, PlotGenerator gen, boolean isStatic ) {
		super( name, comment, DataType.IMAGE );
		this.gen = gen;
		this.isStatic = isStatic;
	}

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
	}

	@Override
	protected void extractFurtherInfos( StringBuffer buf ) throws NonParsableException {
		gen = (PlotGenerator)XMLParser.extractObjectForTags( buf, "generator" );
		isStatic = (Boolean)XMLParser.extractObjectForTags( buf, "isStatic" );
	}

	@Override
	public PlotGenerator getValue() {
		return gen;
	}

	
	public boolean isStatic() {
		return isStatic;
	}

}
