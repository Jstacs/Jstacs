package projects.xanthogenomes;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Rectangle2D;

import de.jstacs.clustering.hierachical.ClusterTree;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.utils.graphics.GraphicsAdaptor;
import projects.xanthogenomes.TALEFamilyBuilder.TALEFamily;


public class FamilyGroupPlotter extends TALEFamilyTreePlotter {

	public static class FamilyGroupPlotGenerator implements PlotGenerator{

		private ClusterTree<TALEFamily> groupTree;
		
		public FamilyGroupPlotGenerator(ClusterTree<TALEFamily> groupTree) {
			this.groupTree = groupTree;
		}
		
		public FamilyGroupPlotGenerator(StringBuffer xml) throws NonParsableException {
			this.groupTree = (ClusterTree<TALEFamily>)XMLParser.extractObjectForTags( xml, "groupTree" );
		}
		
		@Override
		public StringBuffer toXML() {
			StringBuffer xml = new StringBuffer();
			XMLParser.appendObjectWithTags( xml, groupTree, "groupTree" );
			return xml;
		}

		@Override
		public void generatePlot( GraphicsAdaptor adaptor ) throws Exception {
			FamilyGroupPlotter plotter = new FamilyGroupPlotter( 30 );
			
			Graphics2D dummy = adaptor.getGraphics( 10, 10 );
			
			int[] dim = plotter.getDimension( dummy, groupTree );
			
			Graphics2D graphics = adaptor.getGraphics( dim[0], dim[1] );
			
			graphics.setColor( Color.white );
			
			graphics.fillRect( 0, 0, dim[0], dim[1] );
			
			graphics.setColor( Color.black );
			
			graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			graphics.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
			
			plotter.plot( graphics, groupTree );
		}
		
	}
	
	public FamilyGroupPlotter(int lineHeight){
		super(lineHeight);
	}
	
	private double[] getIDFontDimensions(Graphics2D graphics, TALEFamily[] members){
		double fontHeight = 0;
		double fontWidth = 0;
		
		for(int i=0;i<members.length;i++){
			String id = members[i].getFamilyId();
			Rectangle2D rect = graphics.getFontMetrics().getStringBounds( id, graphics );
			if(rect.getHeight() > fontHeight){
				fontHeight = rect.getHeight();
			}
			if(rect.getWidth() > fontWidth){
				fontWidth = rect.getWidth();
			}
		}
		return new double[]{fontWidth,fontHeight};
	}
	
	private void plotIDs(Graphics2D graphics, TALEFamily[] members, int xoff, int yoff){
		
		for(int i=0;i<members.length;i++){
			graphics.drawString( members[i].getFamilyId(), xoff, yoff );
			yoff += 2*lineHeight;
			graphics.setColor( Color.BLACK );
		}
		
	}
	
	public int[] getDimension(Graphics2D graphics, ClusterTree<TALEFamily> tree){
		graphics = (Graphics2D)graphics.create();
		
		graphics.setFont( new Font(Font.MONOSPACED, Font.PLAIN, graphics.getFont().getSize()) );
		
		int height = getHeight( tree.getNumberOfElements() );
		int treeWidth = height;
		double[] idDim = getIDFontDimensions( graphics, tree.getClusterElements() );
		
		double rat = lineHeight/idDim[1];
		
		int idWidth = (int)Math.ceil( idDim[0]*rat );
		
		return new int[]{treeWidth + idWidth + 2*lineHeight, height + (tree.getNumberOfElements() > 1 ? 2*lineHeight : 0)};
	}
	
	public PlotGenerator getPlotGenerator(ClusterTree<TALEFamily> tree){
		return new FamilyGroupPlotGenerator(tree);
	}
	
	public void plot(Graphics2D graphics, ClusterTree<TALEFamily> tree){
		
		graphics = (Graphics2D)graphics.create();
		
		graphics.setFont( new Font(Font.MONOSPACED, Font.PLAIN, graphics.getFont().getSize()) );
		
		TALEFamily[] members = tree.getClusterElements();
		
		int numMem = members.length;
		
		int height = getHeight( numMem );
		
		int treeWidth = height;
		
		
		
		int xoff = treeWidth;
		int yoff = (int) Math.ceil(lineHeight*1.5);
		
		double[] idDim = getIDFontDimensions( graphics, members );
		
		double rat = lineHeight/idDim[1];
		
		int idWidth = (int)Math.ceil( idDim[0]*rat );
		
		graphics.setFont( new Font(graphics.getFont().getFontName(), graphics.getFont().getStyle(), (int)Math.floor(graphics.getFont().getSize()*rat) ) );
		
		plotTree( graphics, treeWidth, tree, 0, 0 );
		
		plotIDs(graphics, members, xoff, yoff);
		
		
	}
	
}
