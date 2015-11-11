package de.jstacs.clustering.hierachical;

import java.awt.BasicStroke;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GraphicsConfiguration;
import java.awt.geom.Rectangle2D;

import projects.motifComp.PWMSupplier;
import de.jstacs.utils.SeqLogoPlotter;


public class TreeDrawer {
	
	public interface LeafDrawer{
		
		public void drawLeaf(Graphics2D graphics, int xoff, int yoff, int height, int margin);
		
		public int getWidth(Graphics2D graphics, int height);
		
	}
	
	private int elementHeight;
	private int margin;
	private int totalHeight;
	private int externalHeight;
	private boolean hang;
	private Coordinates coords;
	private static int gScale = 100;
	private Graphics2D graphics;
	
	public TreeDrawer(int pwmHeight, int margin, int totalHeight, ClusterTree tree, Graphics2D graphics){
		pwmHeight *= gScale;
		margin *= gScale;
		totalHeight *= gScale;
		this.elementHeight = pwmHeight;
		this.externalHeight = totalHeight;
		this.graphics = graphics;
		this.totalHeight = totalHeight-margin/2-elementHeight-elementHeight/2;
		this.margin = margin;
		//System.out.println("maxDist: "+tree.getMaximumDistance());
		//System.out.println("minDist: "+tree.getMinimumDistance());
		coords = new Coordinates( tree, (double)(this.totalHeight-margin/2)/(tree.getMaximumDistance()-Math.min( 0, tree.getMinimumDistance() )), Math.min( 0, tree.getMinimumDistance() ) );
		coords.setX( 0 );
	}
	
	
	public Dimension getDimension(){
		return new Dimension( coords.width/gScale, externalHeight/gScale );
	}
	
	private int[] getLeafDimensions(ClusterTree leaf){
		
		Object element = leaf.getClusterElements()[0];
		
		if(element instanceof PWMSupplier){
			double[][] pwm = ((PWMSupplier)element).getPWM();

			int width = SeqLogoPlotter.getWidth( elementHeight, pwm );
			return new int[]{width+margin,elementHeight};
		}else if(element instanceof LeafDrawer){
			int w = ((LeafDrawer)element).getWidth( graphics, elementHeight );
			return new int[]{w+margin,elementHeight};
		}else{
			return new int[]{elementHeight*element.toString().length(),elementHeight};
		}
	}
	
	private void plotLeaf(ClusterTree leaf, int xoff, int yoff){
		Object element = leaf.getClusterElements()[0];

		if(element instanceof PWMSupplier){
			double[][] pwm = ((PWMSupplier)element).getPWM();

			int width = SeqLogoPlotter.getWidth( elementHeight, pwm );
			graphics.drawRect( xoff, yoff, width, elementHeight );
			SeqLogoPlotter.plotLogo( graphics, xoff+margin/4, yoff+elementHeight-margin/4-margin, width-margin/2, elementHeight-margin/2-margin, pwm, null, "Position", "bits" );

			Font font = new Font(graphics.getFont().getName(),Font.BOLD,(elementHeight-margin/2-margin)/10);
			graphics.setFont( font );
			
			String name = ((PWMSupplier)element).getName();
			Rectangle2D rect = graphics.getFontMetrics().getStringBounds( name, graphics );
			graphics.drawString( name, (float)(xoff+margin/4 + width/2 - rect.getCenterX()), yoff+elementHeight-margin/4 );
		}else if(element instanceof LeafDrawer){
			((LeafDrawer)element).drawLeaf( graphics, xoff, yoff, elementHeight, margin );
		}else{
			String name = leaf.toString();
			Rectangle2D rect = graphics.getFontMetrics().getStringBounds( name, graphics );
			graphics.drawString( name, (float)(xoff+margin/4 + rect.getWidth()/2 - rect.getCenterX()), yoff+elementHeight-margin/4 );
		}
	}
	
	public void plotTree(boolean hang){
		this.hang = hang;
		graphics = (Graphics2D)graphics.create();
		graphics.scale( 1.0/gScale, 1.0/gScale );
		graphics.setStroke( new BasicStroke( elementHeight/30f ) );
		plotTree( null, coords );
	}
	
	private void plotTree(Coordinates parent, Coordinates coords){
		
		if(coords.children == null){
			if(parent != null && hang){
				plotLeaf( coords.tree, coords.x, totalHeight-parent.height+elementHeight/2 );
			}else{
				plotLeaf( coords.tree, coords.x, totalHeight-coords.height+elementHeight/2 );
			}
		}else{
			for(int i=0;i<coords.children.length;i++){
				if(i > 0){
					graphics.drawLine( coords.children[i-1].x+coords.children[i-1].width/2, totalHeight-coords.height, coords.children[i].x+coords.children[i].width/2, totalHeight-coords.height );
				}
				if(coords.children[i].children == null && hang){
					graphics.drawLine( coords.children[i].x+coords.children[i].width/2, totalHeight-coords.height, coords.children[i].x+coords.children[i].width/2, totalHeight-coords.height+elementHeight/2 );
				}else{
					graphics.drawLine( coords.children[i].x+coords.children[i].width/2, totalHeight-coords.height, coords.children[i].x+coords.children[i].width/2, totalHeight-coords.children[i].height+(coords.children[i].children == null ? elementHeight/2 : 0) );
				}
				plotTree( coords, coords.children[i] );
			}
		}
		
	}
	
	private class Coordinates{
		
		private Coordinates[] children;
		private ClusterTree tree;
		private int x;
		private int width;
		private int height;
		private int[] xoffs;
		
		public Coordinates(ClusterTree tree, double scale, double minDist ){
			//System.out.println("scale: "+scale+", minDist: "+minDist);
			if(tree.getNumberOfElements() == 1){
				children = null;
				this.tree = tree;
				int[] dim = getLeafDimensions( tree );
				this.width = dim[0];
				this.height = (int) minDist;
			}else{
				ClusterTree[] subs = tree.getSubTrees();
				children = new Coordinates[subs.length];
				xoffs = new int[subs.length];
				this.width = margin/2;
				for(int i=0;i<subs.length;i++){
					children[i] = new Coordinates( subs[i], scale, minDist );
					xoffs[i] = this.width;//TODO margin
					int width = children[i].width;
					this.width += width;
				}
				this.height = (int)((tree.getDistance()-minDist)*scale);
			}
			
		}
		
		public void setX(int x){
			if(children == null){
				this.x = x;
				//System.out.println("leaf: "+x+", "+height);
			}else{
				for(int i=0;i<children.length;i++){
					this.x = x;
					//System.out.println("this: "+x+", "+height);
					//System.out.println("child "+i+" start");
					children[i].setX( x+xoffs[i] );
					//System.out.println("child "+i+" end");
				}
			}
		}
		
	}
	
}
