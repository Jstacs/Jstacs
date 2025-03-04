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

package projects.motifComp;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Rectangle2D;
import java.text.DecimalFormat;
import java.util.HashMap;

import de.jstacs.clustering.distances.DeBruijnMotifComparison;
import de.jstacs.clustering.hierachical.ClusterTree;
import de.jstacs.clustering.hierachical.PWMSupplier;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.sequenceScores.statisticalModels.StatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.PFMWrapperTrainSM;
import de.jstacs.utils.NiceScale;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter;
import de.jstacs.utils.graphics.GraphicsAdaptor;


public class MotifTreePlotter {

	public static class MotifTreePlotGenerator implements PlotGenerator{

		private ClusterTree<StatisticalModel> tree;
		private int lineHeight;
		private String repID;
		private Pair<double[][],int[][]> pair;
		private boolean labelInnerNodes;
		private Double cutHeight;
		
		public MotifTreePlotGenerator(ClusterTree<StatisticalModel> tree, int lineHeight, int n) throws Exception {
			this(tree,lineHeight,n,null, false, null);
		}
		
		public MotifTreePlotGenerator(ClusterTree<StatisticalModel> tree, int lineHeight, int n, String repID, boolean labelInnerNodes, Double cutHeight) throws Exception {
			this.tree = tree;
			this.lineHeight = lineHeight;
			this.repID = repID;
			this.pair = DeBruijnMotifComparison.getClusterRepresentative(tree, n);
			this.labelInnerNodes = labelInnerNodes;
			this.cutHeight = cutHeight;
		}
		
		public MotifTreePlotGenerator(StringBuffer xml) throws NonParsableException{
			tree = (ClusterTree<StatisticalModel>) XMLParser.extractObjectForTags(xml, "tree");
			lineHeight = (Integer) XMLParser.extractObjectForTags(xml, "lineHeight");
			cutHeight = (Double) XMLParser.extractObjectForTags(xml, "cutHeight");
			labelInnerNodes = (boolean) XMLParser.extractObjectForTags(xml, "labelInnerNodes");
			double[][] rep = (double[][]) XMLParser.extractObjectForTags(xml, "rep");
			int[][] offs = (int[][]) XMLParser.extractObjectForTags(xml, "offs");
			this.pair = new Pair<double[][], int[][]>(rep, offs);
			repID = (String) XMLParser.extractObjectForTags(xml, "repID");
		}
		
		@Override
		public StringBuffer toXML() {
			StringBuffer buf = new StringBuffer();
			XMLParser.appendObjectWithTags(buf, tree, "tree");
			XMLParser.appendObjectWithTags(buf, lineHeight, "lineHeight");
			XMLParser.appendObjectWithTags(buf, cutHeight, "cutHeight");
			XMLParser.appendObjectWithTags(buf, labelInnerNodes, "labelInnerNodes");
			XMLParser.appendObjectWithTags(buf, pair.getFirstElement(),"rep");
			XMLParser.appendObjectWithTags(buf, pair.getSecondElement(), "offs");
			XMLParser.appendObjectWithTags(buf, repID, "repID");
			return buf;
		}

		public double[][] getRepresentative(){
			return pair.getFirstElement();
		}
		
		
		public HashMap<String,int[]> getOffsets() throws CloneNotSupportedException{
			HashMap<String, int[]> res = new HashMap<String, int[]>();
			StatisticalModel[] members = tree.getClusterElements();
			int[][] offs = ArrayHandler.clone(pair.getSecondElement());
			
			int minShift = 0;
			for(int i=0;i<offs.length;i++){
				if(offs[i][0] < minShift){
					minShift = offs[i][0];
				}
			}
			for(int i=0;i<offs.length;i++){
				offs[i][0] -= minShift;
			}
			
			int maxLen = 0;
			for(int i=0;i<offs.length;i++) {
				int len = offs[i][0] + members[i].getLength();
				if(len > maxLen) {
					maxLen = len;
				}
			}
			for(int i=0;i<offs.length;i++) {
				res.put( ((PFMWrapperTrainSM)members[i]).getName(), new int[] {offs[i][0],maxLen-(((PFMWrapperTrainSM)members[i]).getLength()+offs[i][0]),offs[i][1]});
			}
			return res;
		}

		
		@Override
		public void generatePlot(GraphicsAdaptor ga) throws Exception {
			MotifTreePlotter plotter = new MotifTreePlotter(lineHeight);
			
			double[][] rep = pair.getFirstElement();
			int[][] offs = pair.getSecondElement();
			
			rep = cutDown(rep);
			
			int[] dim = plotter.getDimension(tree, rep, offs);
			
		/*	if(dim[0] > 400){
				double fac = 400.0/dim[0];
				dim[0] *= fac;
				dim[1] *= fac;
			}*/
			
			plotter.plot(ga.getGraphics(dim[0], dim[1]), tree, rep, offs, repID, labelInnerNodes, cutHeight);
		}
		
	}
	
	public static double[][] cutDown(double[][] rep) {
		int i=0;
		while( SeqLogoPlotter.getICScale(rep[i]) < 0.125 ){
			i++;
		}
		int j=rep.length-1;
		while( SeqLogoPlotter.getICScale(rep[j]) < 0.125 ){
			j--;
		}
		if(j<=i){
			return rep;
		}else{
			double[][] temp = new double[j-i+1][];
			for(int k=i;k<=j;k++){
				temp[k-i] = rep[k];
			}
			return temp;
		}
	}
	
	protected static DecimalFormat format = new DecimalFormat();
	protected static DecimalFormat labFormat = new DecimalFormat("#.###");
	protected int lineHeight;
	
	public MotifTreePlotter(int lineHeight){
		this.lineHeight = lineHeight;
	}
	
	
	protected int getHeight(int numMem){
		return numMem*lineHeight + (numMem-1)*lineHeight/2 + lineHeight;
	}
	
	
	protected int plotTree(Graphics2D graphics, int treeDim, ClusterTree<StatisticalModel> tree, int xoff, int yoff, boolean labelInnerNodes, Double cutHeight){
		xoff += lineHeight/2;
		graphics = (Graphics2D)graphics.create();
		
		graphics.setStroke(new BasicStroke(lineHeight/40f));
		
		//ClusterTree<TALE> tree = family.getTree();
		
		double minDist = tree.getMinimumDistance();
		double maxDist = tree.getMaximumDistance();
		
		treeDim -= lineHeight/2;
		double rat = (treeDim-lineHeight*2)/(maxDist - minDist)/4;
		
	//	System.out.println("max: "+maxDist+" min: "+minDist+" rat: "+rat+" dim: "+treeDim);
		
		if(cutHeight != null){
			int x = lineHeight+(int)Math.round((maxDist - cutHeight)*rat) + xoff;
			
			Color back = graphics.getColor();
			graphics.setColor(Color.RED);
			graphics.drawLine(x, treeDim+lineHeight, x, 0);
			graphics.setColor(back);
		}
		
		
		int loc = plotEdges(graphics,tree,xoff,0,lineHeight/4, 0,treeDim, maxDist, rat, labelInnerNodes);
		
		//graphics.drawLine( 0, loc, lineHeight, loc );
		
		//System.out.println("family size: "+family.getFamilySize());
		
		
		int textHeight = (int)Math.round( graphics.getFontMetrics().getStringBounds( "0.9", graphics ).getHeight() );
		double fac = lineHeight/3.0/textHeight;
		graphics.setFont(new Font(graphics.getFont().getName(),graphics.getFont().getStyle(), (int)Math.floor(graphics.getFont().getSize()*fac)));	
		
		if(tree.getNumberOfElements() > 1){
			graphics.drawLine( xoff, treeDim+lineHeight, xoff+treeDim/4+lineHeight, treeDim+lineHeight );
		
			if(tree.getNumberOfElements()==2 || minDist == maxDist){
				
				
				
				int x = lineHeight+(int)Math.round((maxDist - tree.getDistance())*rat) + xoff;
				graphics.drawLine( x, treeDim+lineHeight, x, treeDim+lineHeight+lineHeight/6 );
				String label = format.format( tree.getDistance() );
				int textCenter = (int)Math.round( graphics.getFontMetrics().getStringBounds( label, graphics ).getCenterX() );
				graphics.drawString( label, x-textCenter, treeDim+lineHeight+lineHeight/4 + textHeight*2 );
			}else{
				NiceScale scale = new NiceScale( minDist, maxDist );
				scale.setMaxTicks( Math.min( tree.getNumberOfElements(), 5 ) );
				double min = scale.getNiceMin();
				double max = scale.getNiceMax();
				double tick = scale.getTickSpacing();
				//System.out.println(minDist+" "+maxDist+" -> "+min+" "+max+" "+tick);
				
				
				int lastEnd = Integer.MAX_VALUE;
				int firstStart = lineHeight+(int)Math.round((maxDist - max)*rat) + xoff;
				int textCenter = (int)Math.round( graphics.getFontMetrics().getStringBounds( format.format( max ), graphics ).getCenterX() );
				firstStart += (int)(textCenter*1.25);
				
				for(double i=min;i<=max+tick/2.0;i+=tick){
					int x = lineHeight+(int)Math.round((maxDist - i)*rat) + xoff;
					if(x > xoff){
						graphics.drawLine( x, treeDim+lineHeight, x, treeDim+lineHeight+lineHeight/6 );
						String label = format.format( i );
						//int textWidth = (int)Math.round( graphics.getFontMetrics().getStringBounds( label, graphics ).getWidth() );
						textCenter = (int)Math.round( graphics.getFontMetrics().getStringBounds( label, graphics ).getCenterX() );
						
						if((x+textCenter < lastEnd && x-textCenter > firstStart) || 
								i >= max - tick/2.0){
							graphics.drawString( label, x-textCenter, treeDim+lineHeight+lineHeight/3 + textHeight*2 );
							lastEnd = x-(int)(textCenter*1.25);
						}
					}
				}
				
			}
		}
		
		return loc;
		
	}
	
	private int plotEdges(Graphics2D graphics, ClusterTree tree, int xoff, int reloff, int yoff, int numAtop, int treeWidth, double maxDist, double rat, boolean labelInnerNodes){
		//System.out.println(tree.getDistance());
		if(tree.getNumberOfElements() == 1){
			graphics.drawLine( xoff+reloff, yoff+((numAtop)*lineHeight*3)/2 + lineHeight, xoff+(int)(lineHeight*1.5)+(int)Math.round(maxDist*rat), yoff+((numAtop)*lineHeight*3)/2 + lineHeight );
			
			return yoff+((numAtop)*lineHeight*3)/2 + lineHeight;
		}else{
		
			ClusterTree[] subs = tree.getSubTrees();
			
			int minPrev = Integer.MAX_VALUE;
			int maxPrev = 0;
			int newXoff = lineHeight+(int)Math.round((maxDist - tree.getDistance())*rat);
			for(int i=0;i<subs.length;i++){
				int prev = plotEdges( graphics, subs[i], xoff, newXoff, yoff, numAtop, treeWidth, maxDist, rat, labelInnerNodes );
				if(prev < minPrev){
					minPrev = prev;
				}
				if(prev > maxPrev){
					maxPrev = prev;
				}
				numAtop += subs[i].getNumberOfElements();
			}
			graphics.drawLine( xoff+reloff, minPrev+(maxPrev-minPrev)/2, xoff+newXoff, minPrev+(maxPrev-minPrev)/2 );

			graphics.drawLine( xoff+newXoff, minPrev, xoff+newXoff, maxPrev );
			
			if(labelInnerNodes){
				
				Font f = graphics.getFont();
				int size = (int)Math.ceil(lineHeight/10.0);
				graphics.setFont(new Font(graphics.getFont().getName(),graphics.getFont().getStyle(), size));
				String lab = " "+labFormat.format(tree.getDistance());
				graphics.drawString(lab, xoff+newXoff,minPrev+(maxPrev-minPrev)/2 );
				
				graphics.setFont(f);
			}
			
			return minPrev+(maxPrev-minPrev)/2;
			
		}
	}
	
	private void plotLeaves(Graphics2D graphics, StatisticalModel[] members, int[][] offs, int xoff, int yoff, String repID){
		graphics = (Graphics2D) graphics.create();
		yoff +=2*lineHeight;
		
		for(int i=0;i<members.length;i++){
			if(members[i] instanceof PWMSupplier){
				String name = "";
				boolean isRep = false;
				if(members[i] instanceof PFMWrapperTrainSM){
					name = ((PFMWrapperTrainSM)members[i]).getName();
					if(name != null && name.equals(repID)){
						isRep = true;
					}
				}
				
				double[][] pwm = ((PWMSupplier)members[i]).getPWM();
				if(offs[i][1] < 0){
					pwm = PFMComparator.getReverseComplement(DNAAlphabet.SINGLETON, pwm);
					name += " (RC)";
				}
				
				double w = SeqLogoPlotter.getColumnWidth( lineHeight, pwm.length );
				
				int width = SeqLogoPlotter.getWidth(lineHeight, pwm);
				
				SeqLogoPlotter.plotLogo(graphics, (int)(xoff+offs[i][0]*w), yoff, width, lineHeight, pwm, null, "Position", "bits");
				
				if(isRep){
					Stroke stroke = graphics.getStroke();
					graphics.setStroke(new BasicStroke(lineHeight/40f));
					graphics.setColor(Color.RED);
					//graphics.drawLine((int)(xoff+(offs[i][0]-0.5)*w), (int)(yoff), (int)(xoff+(offs[i][0]-0.5)*w), (int)(yoff-lineHeight));
					graphics.drawRect((int)(xoff+(offs[i][0]-0.5)*w), (int)(yoff-lineHeight*1.4), width+(int)w, (int)(lineHeight*1.4));
					graphics.setColor(Color.BLACK);
					graphics.setStroke(stroke);
				}
				
				if(name != null){
					int textWidth = (int)Math.round( graphics.getFontMetrics().getStringBounds( name, graphics ).getWidth() );
					
					double fac = 0.8*width/(double)textWidth;
					
					int size = (int)Math.floor(graphics.getFont().getSize()*fac);
					if(size > lineHeight/5){
						size = lineHeight/5;
					}
					graphics.setFont(new Font(graphics.getFont().getName(),graphics.getFont().getStyle(), size));
										
					int textCenter = (int)Math.round( graphics.getFontMetrics().getStringBounds( name, graphics ).getCenterX() );
					
					
					
					graphics.drawString(name, (int)(xoff+offs[i][0]*w-textCenter+width/2), yoff-(int)Math.floor(lineHeight*1.1));
					
					
					
				}
				
				
			}
			
			yoff += (3*lineHeight)/2;
		}
		
	}
	
	
	
	public int[] getDimension(ClusterTree<StatisticalModel> tree, double[][] rep, int[][] offs){
				
		int height = getHeight( tree.getNumberOfElements() );
		int treeWidth = height;
		
		int minOff = offs[0][0];
		for(int i=1;i<offs.length;i++){
			if(offs[i][0]<minOff){
				minOff = offs[i][0];
			}
		}
		
		int maxWidth = 0;
		StatisticalModel[] models = tree.getClusterElements();
		for(int i=0;i<models.length;i++){
			int temp = models[i].getLength()+(offs[i][0]-minOff);
			if(temp > maxWidth){
				maxWidth = temp;
			}
		}
		
		
		int repWidth = SeqLogoPlotter.getWidth(lineHeight, rep.length);
		maxWidth = SeqLogoPlotter.getWidth(lineHeight, maxWidth);
		
		return new int[]{treeWidth/4 + repWidth + maxWidth +lineHeight*2, height+2*lineHeight };
	}
	
	public void plot(Graphics2D graphics, ClusterTree<StatisticalModel> tree, double[][] rep, int[][] offs, String repID, boolean labelInnerNodes, Double cutHeight){
		graphics = (Graphics2D) graphics.create();
		int[] dim = getDimension(tree, rep, offs);
		graphics.setColor(Color.WHITE);
		graphics.fillRect(0, 0, dim[0], dim[1]);
		graphics.setColor(Color.BLACK);
		
		int minShift = 0;
		for(int i=0;i<offs.length;i++){
			if(offs[i][0] < minShift){
				minShift = offs[i][0];
			}
		}
		for(int i=0;i<offs.length;i++){
			offs[i][0] -= minShift;
		}
		
		graphics = (Graphics2D)graphics.create();		
		
		StatisticalModel[] members = tree.getClusterElements();
		
		int numMem = tree.getNumberOfElements();
		
		int height = getHeight( numMem );
		
		int treeHeight = height;
		
		
		
		int repWidth = SeqLogoPlotter.getWidth(lineHeight, rep.length);
		
		int loc = plotTree(graphics, treeHeight, tree, repWidth, 0, labelInnerNodes, cutHeight);
		
		if(rep != null){
			
			int textWidth = (int)Math.round( graphics.getFontMetrics().getStringBounds( "consensus", graphics ).getWidth() );
			//int w = SeqLogoPlotter.getColumnWidth( lineHeight );
			int width = SeqLogoPlotter.getWidth(lineHeight, rep);
			double fac = 0.8*width/(double)textWidth;
			
			int size = (int)Math.floor(graphics.getFont().getSize()*fac);
			if(size > lineHeight/4){
				size = lineHeight/4;
			}
			graphics.setFont(new Font(graphics.getFont().getName(),graphics.getFont().getStyle(), size));
								
			int textCenter = (int)Math.round( graphics.getFontMetrics().getStringBounds( "consensus", graphics ).getCenterX() );
			
			graphics.drawString("consensus", lineHeight/4-textCenter+width/2, loc+lineHeight/2-(int)Math.floor(lineHeight*1.1));
			
			SeqLogoPlotter.plotLogo(graphics, lineHeight/4, loc+lineHeight/2, repWidth, lineHeight, rep, null, "Position", "bits");
		}
		
		plotLeaves(graphics, members, offs, repWidth+treeHeight/4+(int)(lineHeight*1.5), 0, repID);

		
	}
	
	
}
