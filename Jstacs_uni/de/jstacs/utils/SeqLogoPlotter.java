package de.jstacs.utils;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import javax.imageio.ImageIO;

import de.jstacs.data.alphabets.DNAAlphabet;

/**
 * Class with static methods for plotting sequence logos of DNA motifs, i.e., position weight matrices defined over a {@link DNAAlphabet}.
 * In general, sequence logos can be plotted to any {@link Graphics2D} object, e.g., for on-screen printing using the <code>plotLogo</code> methods.
 * 
 * For convenience, the method {@link SeqLogoPlotter#plotLogoToPNG(String, int, double[][])} can be used to directly plot a sequence logo
 * to a PNG file with a given height and automatically chosen aspect ratio.
 * 
 * @author Jan Grau
 *
 */
public class SeqLogoPlotter {

	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>.
	 * 
	 * The positions of the sequence logo are numbered contiguously from 1 to <code>ps.length</code>. The label of the
	 * x-axis is set to &quot;Position&quot;, and the label of the y-axis is set to &quot;bits&quot;.
	 * 
	 * The sequence logo is written to the PNG file given in <code>path</code>.
	 * 
	 * @param path the path to the PNG file written
	 * @param height the height of the PNG image (in pixels)
	 * @param ps the position weight matrix
	 * @throws IOException if the file could not be written
	 */
	public static void plotLogoToPNG(String path, int height, double[][] ps) throws IOException{
		Pair<BufferedImage, Graphics2D> pair = getBufferedImageAndGraphics( height, ps );
		Graphics2D g = pair.getSecondElement();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		plotLogo( g, height, ps );
		ImageIO.write( pair.getFirstElement(), "png", new File(path) );
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>.
	 * 
	 * The positions of the sequence logo are numbered contiguously from 1 to <code>ps.length</code>. The label of the
	 * x-axis is set to &quot;Position&quot;, and the label of the y-axis is set to &quot;bits&quot;.
	 * 
	 * The sequence logo is return as a {@link BufferedImage}.
	 * 
	 * @param height the height of the PNG image (in pixels)
	 * @param ps the position weight matrix
	 * @return the sequence logo
	 */
	public static BufferedImage plotLogoToBufferedImage(int height, double[][] ps) {
		Pair<BufferedImage, Graphics2D> pair = getBufferedImageAndGraphics( height, ps );
		Graphics2D g = pair.getSecondElement();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		plotLogo( g, height, ps );
		return pair.getFirstElement();
	}
	
	/**
	 * Creates a new {@link BufferedImage} with given height and width chosen automatically according to the number of rows
	 * of <code>ps</code>, and returns this {@link BufferedImage} and its {@link Graphics2D} object.
	 * 
	 * @param height the height (in pixels)
	 * @param ps the position weight matrix
	 * @return the created {@link BufferedImage} and its {@link Graphics2D} object
	 */
	protected static Pair<BufferedImage,Graphics2D> getBufferedImageAndGraphics(int height,double[][] ps){
		int w = getWidth( height, ps );
		BufferedImage img = new BufferedImage( w, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = (Graphics2D)img.getGraphics();
		
		return new Pair<BufferedImage, Graphics2D>( img, g );
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>.
	 * 
	 * The positions of the sequence logo are numbered contiguously from 1 to <code>ps.length</code>. The label of the
	 * x-axis is set to &quot;Position&quot;, and the label of the y-axis is set to &quot;bits&quot;.
	 * 
	 * The sequence logo is written to the {@link Graphics2D} object given in <code>g</code>.
	 * 
	 * @param g the {@link Graphics2D} object
	 * @param h the height of the sequence logo
	 * @param ps the position weight matrix
	 */
	public static void plotLogo(Graphics2D g, int h, double[][] ps){
		plotLogo( g, h, ps, null, "Position", "bits" );
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>.
	 * 
	 * The sequence logo is written to the {@link Graphics2D} object given in <code>g</code>.
	 * 
	 * @param g the {@link Graphics2D} object
	 * @param height the height of the sequence logo
	 * @param ps the position weight matrix
	 * @param labels the labels of the positions of the sequence logo, if <code>null</code> the positions are numbered contiguously from 1 to <code>ps.length</code>
	 * @param labX the label of the x-axis
	 * @param labY the label of the y-axis
	 */
	public static void plotLogo(Graphics2D g, int height, double[][] ps, String[] labels, String labX, String labY){
		int w = getWidth( height, ps );
		plotLogo( g, w, height, ps, labels, labX, labY );
	}
	
	/**
	 * Returns the automatically chosen width for a given height and position weight matrix.
	 * @param height the height
	 * @param ps the position weight matrix
	 * @return the width
	 */
	public static int getWidth(int height, double[][] ps){
		return (int)(height/6.0*(ps.length+1.5));
	}
	
	/**
	 * Returns the automatically chosen height for a given width and position weight matrix.
	 * @param width the width
	 * @param ps the position weight matrix
	 * @return the height
	 */
	public static int getHeight(int width, double[][] ps){
		return (int)(width*6.0/(ps.length+1.5));
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * 
	 * The sequence logo is written to the {@link Graphics2D} object given in <code>g</code>.
	 * 
	 * @param g the {@link Graphics2D} object
	 * @param w the width of the sequence logo
	 * @param h the height of the sequence logo
	 * @param ps the position weight matrix
	 * @param labels the labels of the positions of the sequence logo, if <code>null</code> the positions are numbered contiguously from 1 to <code>ps.length</code>
	 * @param labX the label of the x-axis
	 * @param labY the label of the y-axis
	 */
	public static void plotLogo(Graphics2D g, int w, int h, double[][] ps, String[] labels, String labX, String labY){
		plotLogo( g, 0, h, w, h, ps, labels, labX, labY );
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * 
	 * The sequence logo is written to the {@link Graphics2D} object given in <code>g</code>.
	 * 
	 * @param g the {@link Graphics2D} object
	 * @param x the x-coordinate of the bottom left corner of the sequence logo
	 * @param y the y-coordinate of the bottom left corner of the sequence logo (<code>-h</code> results in a sequence logo spanning from <code>0</code> to <code>h</code>)
	 * @param w the width of the sequence logo
	 * @param h the height of the sequence logo
	 * @param ps the position weight matrix
	 * @param labels the labels of the positions of the sequence logo, if <code>null</code> the positions are numbered contiguously from 1 to <code>ps.length</code>
	 * @param labX the label of the x-axis
	 * @param labY the label of the y-axis
	 */
	public static void plotLogo(Graphics2D g, int x, int y, int w, int h, double[][] ps, String[] labels, String labX, String labY){
		g.setColor( Color.WHITE );
		g.fillRect( x, y-h, w, h );
		
		Font font = new Font(g.getFont().getName(),Font.BOLD,h/17);
		g.setFont( font );
		
		if(labels == null){
			labels = new String[ps.length];
			for(int i=0;i<ps.length;i++){
				labels[i] = (i+1)+"";
			}
		}
		

		double wl = h*0.25;
		
		double w2 = (w-wl)/(ps.length);
		double x2 = x + wl*0.9;
		
		double h2 = h*0.7;
		double y2 = y*0.75;
		
		g.setColor( Color.BLACK );
		g.setStroke( new BasicStroke( h/400+1 ) );
		g.drawLine( (int)x2, (int)(y2*1.04)+1, (int)(x2+w2*ps.length), (int)(y2*1.04)+1 );
		g.drawLine( (int)(x2*0.94)-1, (int)y2, (int)(x2*0.94)-1, (int)(y2-h2) );
		String[] labs = {"0", "0.5", "1", "1.5", "2"};
		for(int i=0;i<=4;i++){
			g.drawLine( (int)(x2*0.7), (int)(y2-i*h2/4.0), (int)(x2*0.94)-1, (int)(y2-i*h2/4.0) );
			Rectangle2D rect = g.getFontMetrics().getStringBounds( labs[i], g );
			g.drawString( labs[i], (int)(x2*0.6-rect.getWidth()-2), (int)(y2-i*h2/4.0 - rect.getCenterY()) );
		}
		AffineTransform back = g.getTransform();
		g.rotate( -Math.PI/2 );
		Rectangle2D rect = g.getFontMetrics().getStringBounds( labY, g );

		g.drawString(labY,-(int)(y2-2*h2/4.0 + rect.getCenterX()), (int)(x+rect.getWidth()/2d));
		g.setTransform( back );
		
		rect = g.getFontMetrics().getStringBounds( labX, g );
		g.drawString( labX, (int)(x2+w2*ps.length/2.0-rect.getCenterX()), (int)(y-0.3*rect.getHeight()) );
		for(int i=0;i<ps.length;i++){
			plotLogo( g, x2, y2, w2, h2, ps[i] );
			g.setColor( Color.BLACK );
			rect = g.getFontMetrics().getStringBounds( labels[i], g );
			g.drawString( labels[i], (float)(x2+w2/2d-rect.getCenterX()), (float)(y2+2*rect.getHeight()) );
			x2 += w2;
		}
	}
	
	/**
	 * Plots the TALgetter logo for the binding specificities given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to specificities of the RVDs. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>. In addition, the importance of RVDs is plotted as given in <code>imp</code>.
	 * 
	 * The labels of the RVDs are given in <code>lab</code>. The label of the
	 * x-axis is set to &quot;RVD&quot;, and the label of the y-axes are set to &quot;bits&quot; and &quot;Importance&quot;, respectively.
	 * 
	 * The TALgetter logo is written to the PNG file given in <code>path</code>.
	 * 
	 * @param path the path to the PNG file written
	 * @param height the height of the PNG image (in pixels)
	 * @param ps the binding specificities of RVDs
	 * @param imp the importance of RVDs
	 * @param lab the amino acids of the RVDs in one-letter code
	 * @throws IOException if the file could not be written
	 */
	public static void plotTALgetterLogoToPNG(String path, int height, double[][] ps, double[] imp, String[] lab) throws IOException{
		Pair<BufferedImage, Graphics2D> pair = getBufferedImageAndGraphics( height, ps );
		Graphics2D g = pair.getSecondElement();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		int w = getWidth( height, ps );
		
		plotTALgetterLogo( g, 0, height, w, height, ps, imp, lab, "RVD", "bits", "Importance" );
		ImageIO.write( pair.getFirstElement(), "png", new File(path) );
	}
	
	/**
	 * Plots the TALgetter logo for the binding specificities given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to specificities of the RVDs. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>. In addition, the importance of RVDs is plotted as given in <code>imp</code>.
	 * 
	 * The labels of the RVDs are given in <code>lab</code>. The label of the
	 * x-axis is set to &quot;RVD&quot;, and the label of the y-axes are set to &quot;bits&quot; and &quot;Importance&quot;, respectively.
	 * 
	 * The TALgetter logo is returned as {@link BufferedImage}.
	 * 
	 * @param path the path to the PNG file written
	 * @param height the height of the PNG image (in pixels)
	 * @param ps the binding specificities of RVDs
	 * @param imp the importance of RVDs
	 * @param lab the amino acids of the RVDs in one-letter code
	 * @return the TALgetter logo
	 */
	public static BufferedImage plotTALgetterLogoToBufferedImage(int height, double[][] ps, double[] imp, String[] lab) {
		Pair<BufferedImage, Graphics2D> pair = getBufferedImageAndGraphics( height, ps );
		Graphics2D g = pair.getSecondElement();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		int w = getWidth( height, ps );
		
		plotTALgetterLogo( g, 0, height, w, height, ps, imp, lab, "RVD", "bits", "Importance" );
		return pair.getFirstElement();
	}
	
	/**
	 * Plots the TALgetter logo for the binding specificities given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to specificities of the RVDs. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>. In addition, the importance of RVDs is plotted as given in <code>imp</code>.
	 * 
	 * The labels of the RVDs are given in <code>lab</code>. The label of the
	 * x-axis is set to &quot;RVD&quot;, and the label of the y-axes are set to &quot;bits&quot; and &quot;Importance&quot;, respectively.
	 * 
	 * The TALgetter logo is returned as {@link BufferedImage}.
	 * 
	 * @param g the {@link Graphics2D} object
	 * @param x the x-coordinate of the bottom left corner of the TALgetter logo
	 * @param y the y-coordinate of the bottom left corner of the TALgetter logo (<code>-h</code> results in a logo spanning from <code>0</code> to <code>h</code>)
	 * @param w the width of the TALgetter logo
	 * @param h the height of the TALgetter logo
	 * @param ps the binding specificities of RVDs
	 * @param imp the importance of RVDs
	 * @param labels the amino acids of the RVDs in one-letter code
	 * @param labX the label of the x-axis
	 * @param labY the label of the y-axis
	 * @param labY2 the label of the second (importance) y-axis
	 */
	public static void plotTALgetterLogo(Graphics2D g, int x, int y, int w, int h, double[][] ps, double[] imp, String[] labels, String labX, String labY, String labY2){
		g.setColor( Color.WHITE );
		g.fillRect( x, y-h, w, h );
		

		Font font = new Font(g.getFont().getName(),Font.BOLD,h/17);
		g.setFont( font );
		
		if(labels == null){
			labels = new String[ps.length];
			for(int i=0;i<ps.length;i++){
				labels[i] = (i+1)+"";
			}
		}
		

		double wl = h*0.5;
		
		double w2 = (w-wl)/(ps.length);
		double x2 = x + wl*0.45;
		
		double h2 = h*0.7;
		double y2 = y*0.75;
		
		g.setColor( Color.BLACK );
		g.setStroke( new BasicStroke( h/400+1 ) );
		g.drawLine( (int)x2, (int)(y2*1.04)+1, (int)(x2+w2*ps.length), (int)(y2*1.04)+1 );
		g.drawLine( (int)(x2*0.94)-1, (int)y2, (int)(x2*0.94)-1, (int)(y2-h2) );
		String[] labs = {"0", "0.5", "1", "1.5", "2"};
		for(int i=0;i<=4;i++){
			g.drawLine( (int)(x2*0.7), (int)(y2-i*h2/4.0), (int)(x2*0.94)-1, (int)(y2-i*h2/4.0) );
			Rectangle2D rect = g.getFontMetrics().getStringBounds( labs[i], g );
			g.drawString( labs[i], (int)(x2*0.6-rect.getWidth()-2), (int)(y2-i*h2/4.0 - rect.getCenterY()) );
		}
		AffineTransform back = g.getTransform();
		g.rotate( -Math.PI/2 );
		Rectangle2D rect = g.getFontMetrics().getStringBounds( labY, g );

		g.drawString(labY,-(int)(y2-2*h2/4.0 + rect.getCenterX()), (int)(x+rect.getHeight()));
		g.setTransform( back );
		
		
		rect = g.getFontMetrics().getStringBounds( labX, g );
		g.drawString( labX, (int)(x2+w2*ps.length/2.0-rect.getCenterX()), (int)(y-0.3*rect.getHeight()) );
		for(int i=0;i<ps.length;i++){
			plotLogo( g, x2, y2, w2, h2, ps[i] );
			g.setColor( Color.BLACK );
			rect = g.getFontMetrics().getStringBounds( labels[i], g );
			g.drawString( labels[i], (float)(x2+w2/2d-rect.getCenterX()), (float)(y2+2*rect.getHeight()) );
			x2 += w2;
		}
		
		g.drawLine( (int)(x2+w2*0.1)+1, (int)y2, (int)(x2+w2*0.1)+1, (int)(y2-h2) );
		labs = new String[]{"0", "0.5", "1"};
		Rectangle2D rect2 = g.getFontMetrics().getStringBounds( "0.5", g );
		for(int i=0;i<=2;i++){
			g.drawLine( (int)(x2+w2*0.34)+2, (int)(y2-i*h2/2.0), (int)(x2+w2*0.1)+1, (int)(y2-i*h2/2.0) );
			//rect = g.getFontMetrics().getStringBounds( labs[i], g );
			g.drawString( labs[i], (int)((x2+w2*0.5)+2), (int)(y2-i*h2/2.0 - rect.getCenterY()) );
		}
		
		back = g.getTransform();
		g.rotate( -Math.PI/2 );
		rect = g.getFontMetrics().getStringBounds( labY2, g );

		g.drawString(labY2,-(int)(y2-2*h2/4.0 + rect.getCenterX()), (int)(w-rect.getHeight()/2d));
		g.setTransform( back );
		
		x2 = x + wl*0.45 + 2*w2;
		

		g.setColor( Color.GRAY );
		for(int i=1;i<imp.length;i++){
			g.drawLine( (int)(x2-w2/2d + w2/20d),(int)(y2-h2*imp[i-1]+w2/20d) , (int)(x2+w2/2d - w2/20d), (int)(y2-h2*imp[i]+w2/20d) );
			x2 += w2;
		}
		
		x2 = x + wl*0.45 + w2;
		
		g.setColor( Color.BLUE );
		for(int i=0;i<imp.length;i++){
			
			g.fillRect( (int)(x2+w2/2d - w2/20d), (int)(y2-h2*imp[i]) , (int)(w2/10d), (int)(w2/10d) );
			x2 += w2;
		}
		
	}
	
	private static void plotLogo(Graphics2D g, double x, double y, double w, double h, double[] p){
		
		//y += h;
		
		double ic = Math.log( 4 )/Math.log( 2 );
		for(int i=0;i<p.length;i++){
			if(p[i] > 0){
				ic += p[i]*Math.log( p[i] )/Math.log( 2 );
			}
		}
		ic /= 2.0;
		h *= ic;
		//System.out.println("h: "+h);
		
		double[] mp = p.clone();
		for(int i=0;i<mp.length;i++){
			mp[i] *= -1;
		}
		
		int[] r = ToolBox.rank( mp, false );
		int[] order = new int[r.length];
		for(int i=0;i<r.length;i++){
			order[r[i]] = i;
		}
		
		
		for(int i=0;i<order.length;i++){
			double curr = p[order[i]];
			if(order[i] == 0){
				g.setColor( Color.GREEN );
				g.fill(getA( x, y, w, h*curr ));
			}else if(order[i] == 1){
				g.setColor( Color.BLUE );
				g.fill(getC( x, y, w, h*curr ));
			}else if(order[i] == 2){
				g.setColor( Color.ORANGE );
				g.fill(getG( x, y, w, h*curr ));
			}else{
				g.setColor( Color.RED );
				g.fill(getT( x, y, w, h*curr ));
			}
			//System.out.println("y: "+y);
			y -= h*curr;
		}
		//System.out.println(y);
	}
	
	private static Area getC( double x, double y, double w, double h ){
		
		Shape s = new Ellipse2D.Double( 0, -90, 90, 90 );
		Area a1 = new Area( s );
		
		Shape s2 = new Ellipse2D.Double( 15, -75, 60, 60 );
		Area a2 = new Area(s2);
		
		a1.subtract( a2 );
		
		Shape s3 = new Rectangle2D.Double( 65, -60, 30, 30 );
		Area a3 = new Area( s3 );
		
		a1.subtract( a3 );
		
		AffineTransform t = new AffineTransform();
		t.scale( 1d/88d, 1d/90d );
		t.scale( w, h );
		a1.transform( t );
		t = new AffineTransform();
		t.translate( x, y );
		a1.transform( t );
		return a1;
	}
	
	private static Area getT( double x, double y, double w, double h ){
		
		Shape s = new Rectangle2D.Double( 37.5, -100, 15, 100 );
		Area a1 = new Area( s );
		
		Shape s2 = new Rectangle2D.Double(0,-100,90,15);
		Area a2 = new Area( s2 );
		
		a1.add( a2 );
		
		AffineTransform t = new AffineTransform();
		t.scale( 1d/90d, 1d/100d );
		t.scale( w, h );
		a1.transform( t );
		t = new AffineTransform();
		t.translate( x, y );
		a1.transform( t );
		return a1;
		
	}
	
	private static Area getG( double x, double y, double w, double h ){
		
		Shape s = new Ellipse2D.Double( 0, -90, 90, 90 );
		Area a1 = new Area( s );
		
		Shape s2 = new Ellipse2D.Double( 15, -75, 60, 60 );
		Area a2 = new Area(s2);
		
		a1.subtract( a2 );
		
		Shape s3 = new Rectangle2D.Double( 65, -60, 30, 30 );
		Area a3 = new Area( s3 );
		
		a1.subtract( a3 );
		
		Shape s4 = new Rectangle2D.Double( 55, -40, 35, 15 );
		Area a4 = new Area( s4 );
		
		a1.add( a4 );
		
		Shape s5 = new Rectangle2D.Double( 80, -40, 10, 40 );
		Area a5 = new Area( s5 );
		
		a1.add( a5 );
		
		AffineTransform t = new AffineTransform();
		t.scale( 1d/90d, 1d/90d );
		t.scale( w, h );
		a1.transform( t );
		t = new AffineTransform();
		t.translate( x, y );
		a1.transform( t );
		return a1;
	}
	
	private static Area getA( double x, double y, double w, double h ){
		
		Shape s = new Polygon( new int[]{0,40,50,90,75,45,45,15,0}, new int[]{0,-100,-100,0,0,-80,-80,0,0}, 9 );
		Area a = new Area( s );
		
		Shape s2 = new Polygon(new int[]{20,70,70,20},new int[]{-35,-35,-50,-50},4);
		Area a2 = new Area( s2 );
		
		a.add( a2 );
		
		AffineTransform t = new AffineTransform();
		t.scale( 1.0/90.0, 1.0/100.0 );
		t.scale( w, h );
		a.transform( t );
		t = new AffineTransform();
		t.translate( x, y );
		a.transform( t );
		return a;
	}
	
}
