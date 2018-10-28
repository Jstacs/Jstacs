package projects.xanthogenomes;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.text.DecimalFormat;

import de.jstacs.clustering.hierachical.ClusterTree;
import de.jstacs.utils.NiceScale;
import de.jstacs.utils.Pair;
import projects.xanthogenomes.TALE.Type;
import projects.xanthogenomes.TALEFamilyBuilder.TALEFamily;


public class TALEFamilyTreePlotter {

	protected static DecimalFormat format = new DecimalFormat();
	protected int lineHeight;
	
	public TALEFamilyTreePlotter(int lineHeight){
		this.lineHeight = lineHeight;
	}
	
	
	protected int getHeight(int numMem){
		return numMem*lineHeight + (numMem-1)*lineHeight + lineHeight;
	}
	
	private double[] getIDFontDimensions(Graphics2D graphics, TALE[] members){
		double fontHeight = 0;
		double fontWidth = 0;
		
		for(int i=0;i<members.length;i++){
			String id = members[i].getId();
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
	
	private double[] getAlignmentFontDimensions(Graphics2D graphics, TALEFamily family){
		Pair<TALE[],String[]> pair = family.getInducedMultipleAlignment();
		String[] als = pair.getSecondElement();
		
		Rectangle2D rect = graphics.getFontMetrics().getStringBounds( als[0], graphics );
		return new double[]{rect.getWidth(),rect.getHeight()};
	}
	
	protected void plotTree(Graphics2D graphics, int treeDim, ClusterTree tree, int xoff, int yoff){
		
		graphics = (Graphics2D)graphics.create();
		graphics.setStroke( new BasicStroke( lineHeight/15f ) );
		graphics.setFont( new Font(Font.SANS_SERIF,Font.PLAIN,(graphics.getFont().getSize()*2)/3) );
		
		//ClusterTree<TALE> tree = family.getTree();
		
		double minDist = tree.getMinimumDistance();
		double maxDist = tree.getMaximumDistance();
		
		treeDim -= lineHeight/2;
		double rat = (treeDim-lineHeight*2)/(maxDist - minDist);
		
	//	System.out.println("max: "+maxDist+" min: "+minDist+" rat: "+rat+" dim: "+treeDim);
		
		int loc = plotEdges(graphics,tree,xoff+lineHeight,lineHeight/4, 0,treeDim, maxDist, rat);
		
		graphics.drawLine( 0, loc, lineHeight, loc );
		
		//System.out.println("family size: "+family.getFamilySize());
		
		if(tree.getNumberOfElements() > 1){
			graphics.drawLine( lineHeight/2, treeDim+lineHeight, treeDim-lineHeight/2, treeDim+lineHeight );
		
			if(tree.getNumberOfElements()==2 || minDist == maxDist){
				int x = lineHeight+(int)Math.round((maxDist - tree.getDistance())*rat);
				graphics.drawLine( x, treeDim+lineHeight, x, treeDim+lineHeight+lineHeight/4 );
				String label = format.format( tree.getDistance() );
				int textCenter = (int)Math.round( graphics.getFontMetrics().getStringBounds( label, graphics ).getCenterX() );
				graphics.drawString( label, x-textCenter, treeDim+lineHeight*2 );
			}else{
				NiceScale scale = new NiceScale( minDist, maxDist );
				scale.setMaxTicks( Math.min( tree.getNumberOfElements(), 5 ) );
				double min = scale.getNiceMin();
				double max = scale.getNiceMax();
				double tick = scale.getTickSpacing();
				//System.out.println(minDist+" "+maxDist+" -> "+min+" "+max+" "+tick);
				
				for(double i=min;i<=max+tick/2.0;i+=tick){
					if(i<1E-6 && i>-1E-6){
						i=0;
					}
					//System.out.println(i+" "+(i+tick)+" <= "+max+", "+(max+tick/2.0));
					int x = lineHeight+(int)Math.round((maxDist - i)*rat);
					if(x > lineHeight/2){
						graphics.drawLine( x, treeDim+lineHeight, x, treeDim+lineHeight+lineHeight/4 );
						String label = format.format( i );
						int textCenter = (int)Math.round( graphics.getFontMetrics().getStringBounds( label, graphics ).getCenterX() );
						graphics.drawString( label, x-textCenter, treeDim+lineHeight*2 );
					}
				}
				
			}
		}
		
		
	}
	
	private int plotEdges(Graphics2D graphics, ClusterTree tree, int xoff, int yoff, int numAtop, int treeWidth, double maxDist, double rat){
		
		if(tree.getNumberOfElements() == 1){
			graphics.drawLine( xoff, yoff+(numAtop)*lineHeight*2 + lineHeight, treeWidth, yoff+(numAtop)*lineHeight*2 + lineHeight );
			
			return yoff+(numAtop)*lineHeight*2 + lineHeight;
		}else{
		
			ClusterTree[] subs = tree.getSubTrees();
			
			int minPrev = Integer.MAX_VALUE;
			int maxPrev = 0;
			int newXoff = lineHeight+(int)Math.round((maxDist - tree.getDistance())*rat);
			for(int i=0;i<subs.length;i++){
				int prev = plotEdges( graphics, subs[i], newXoff, yoff, numAtop, treeWidth, maxDist, rat );
				if(prev < minPrev){
					minPrev = prev;
				}
				if(prev > maxPrev){
					maxPrev = prev;
				}
				numAtop += subs[i].getNumberOfElements();
			}
			
			graphics.drawLine( xoff, minPrev+(maxPrev-minPrev)/2, newXoff, minPrev+(maxPrev-minPrev)/2 );
			graphics.drawLine( newXoff, minPrev, newXoff, maxPrev );
			
			return minPrev+(maxPrev-minPrev)/2;
			
		}
	}
	
	private void plotIDs(Graphics2D graphics, TALE[] members, int xoff, int yoff){
		
		for(int i=0;i<members.length;i++){
			if(members[i].isNew()){
				graphics.setColor( Color.BLUE );
			}
			graphics.drawString( members[i].getId(), xoff, yoff );
			yoff += 2*lineHeight;
			graphics.setColor( Color.BLACK );
		}
		
	}
	
	private void plotTALEAlignment(Graphics2D graphics, TALEFamily family, int xoff, int yoff){
		Pair<TALE[],String[]> pair = family.getInducedMultipleAlignment();
		
		String[] als = pair.getSecondElement();
		TALE[] members = pair.getFirstElement();
		
		String space = " ";
		int spacew = graphics.getFontMetrics().stringWidth(space);
		for(int i=0;i<als.length;i++){
			if(members[i].isNew()){
				graphics.setColor( Color.BLUE );
			}
			//graphics.drawString( als[i], xoff, yoff );
			String[] al = als[i].trim().split(" ");
			int currOff = xoff+spacew;
			for(int j=0,r=0;j<al.length;j++){
				Color bef = graphics.getColor();
				if(! "--".equals(al[j]) && r < members[i].getNumberOfRepeats()-1){
					//System.out.println(j+" "+al[j]);
					if(members[i].getRepeat(r).getType() == Type.LONG){
						graphics.setColor(Color.RED);
					}else if(members[i].getRepeat(r).getType() == Type.SHORT){
						graphics.setColor(Color.GREEN);
					}else if(members[i].getRepeat(r).getType() == Type.THIRTYFIVE){
						graphics.setColor(Color.GRAY);
					}
					r++;
				}
				graphics.drawString(al[j], currOff, yoff);
				int w = graphics.getFontMetrics().stringWidth(al[j]);
				currOff += w+spacew;
				graphics.setColor(bef);
			}
			graphics.setColor( Color.BLACK );
			yoff += lineHeight;
			if(i<als.length-1){
				/*String s1 = als[i];
				String s2 = als[i+1];
				StringBuffer comp = new StringBuffer();
				for(int j=0;j<s1.length();j++){
					if(s1.charAt( j ) == ' ' || s2.charAt( j ) == ' ' || s1.charAt( j ) == '-' || s2.charAt( j ) == '-'){
						comp.append( " " );
					}else if(s1.charAt( j ) == s2.charAt( j )){
						comp.append( "|" );
					}else{
						comp.append( ":" );
					}
				}
				graphics.drawString( comp.toString(), xoff, yoff );*/
				String[] s1 = als[i].trim().split(" ");
				String[] s2 = als[i+1].trim().split(" ");
				currOff = xoff+spacew;
				for(int j=0,r1=0,r2=0;j<s1.length;j++){
					for(int k=0;k<2;k++){
						Color bef = graphics.getColor();
						if(s1[j].charAt(k) == '-' || s2[j].charAt(k) == '-'){
							currOff += spacew;
						}else if( s1[j].charAt(k) == s2[j].charAt(k)){
							String cod1 = members[i].getCodon(r1,k);
							String cod2 = members[i+1].getCodon(r2,k);
							if(cod1 != null && cod2 != null && !cod1.equals(cod2)){
								graphics.setColor(Color.RED);
							}
							graphics.drawString("|", currOff, yoff);
							currOff += graphics.getFontMetrics().stringWidth("|");
						}else{
							//graphics.setColor(Color.RED);
							graphics.drawString(":", currOff, yoff);
							currOff += graphics.getFontMetrics().stringWidth(":");
						}
						graphics.setColor(bef);
					}
					if(! "--".equals(s1[j])){
						r1++;
					}
					if(! "--".equals(s2[j])){
						r2++;
					}
					currOff += spacew;
				}
				yoff += lineHeight;
			}
		}
	}
	
	public int[] getDimension(Graphics2D graphics, TALEFamily family){
		graphics = (Graphics2D)graphics.create();
		
		graphics.setFont( new Font(Font.MONOSPACED, Font.PLAIN, graphics.getFont().getSize()) );
		
		int height = getHeight( family.getFamilySize() );
		int treeWidth = height;
		double[] idDim = getIDFontDimensions( graphics, family.getFamilyMembers() );
		
		double rat = lineHeight/idDim[1];
		
		int idWidth = (int)Math.ceil( idDim[0]*rat );
		
		double[] alDim = getAlignmentFontDimensions( graphics, family );
		
		return new int[]{treeWidth + idWidth + (int)Math.ceil( alDim[0]*rat ) + 2*lineHeight, height + (family.getFamilySize() > 1 ? 2*lineHeight : 0)};
	}
	
	public void plot(Graphics2D graphics, TALEFamily family){
		
		graphics = (Graphics2D)graphics.create();
		
		graphics.setFont( new Font(Font.MONOSPACED, Font.PLAIN, graphics.getFont().getSize()) );
		
		ClusterTree<TALE> tree = family.getTree();
		
		TALE[] members = tree.getClusterElements();
		
		int numMem = family.getFamilySize();
		
		int height = getHeight( numMem );
		
		int treeWidth = height;
		
		
		
		int xoff = treeWidth;
		int yoff = (int) Math.ceil(lineHeight*1.5);
		
		double[] idDim = getIDFontDimensions( graphics, members );
		
		double rat = lineHeight/idDim[1];
		
		int idWidth = (int)Math.ceil( idDim[0]*rat );
		
		graphics.setFont( new Font(graphics.getFont().getFontName(), graphics.getFont().getStyle(), (int)Math.floor(graphics.getFont().getSize()*rat) ) );
		
		plotTree( graphics, treeWidth, family.getTree(), 0, 0 );
		
		plotIDs(graphics, members, xoff, yoff);
		
		xoff += idWidth;
		
		plotTALEAlignment( graphics, family, xoff, yoff );
		
		
	}
	
	
}
