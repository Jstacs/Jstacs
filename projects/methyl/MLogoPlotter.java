package projects.methyl;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Locale;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.utils.NiceScale;
import de.jstacs.utils.SeqLogoPlotter;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.SeqLogoPlotter.SeqLogoPlotGenerator;
import de.jstacs.utils.graphics.GraphicsAdaptor;

public class MLogoPlotter extends SeqLogoPlotter {

	
	/**
	 * {@link PlotGenerator} for plotting sequence logos.
	 * 
	 * @author Jan Grau
	 *
	 */
	public static class MLogoPlotGenerator implements PlotGenerator{

		private double[][] pwm;
		private double[] CpG;
		private double[] MpH;
		private int height;
		
		/**
		 * Creates a new {@link MLogoPlotGenerator} for the given PWM using the specified height of the plot.
		 * The corresponding width is computed automatically using the {@link SeqLogoPlotter#getWidth(int, double[][])} method.
		 * 
		 * @param pwm the PWM
		 * @param height the height of the plot
		 */
		public MLogoPlotGenerator(double[][] pwm, double[] CpG, double[] MpH, int height){
			this.pwm = pwm;
			this.CpG = CpG;
			this.MpH = MpH;
			this.height = height;
		}
		
		/**
		 * Creates a {@link MLogoPlotGenerator} from its XML representation.
		 * @param xml the XML representation
		 * @throws NonParsableException if XML could not be parsed
		 */
		public MLogoPlotGenerator(StringBuffer xml) throws NonParsableException{
			this.pwm = (double[][]) XMLParser.extractObjectForTags(xml, "pwm");
			this.CpG = (double[]) XMLParser.extractObjectForTags(xml, "CpG");
			this.MpH = (double[]) XMLParser.extractObjectForTags(xml, "MpH");
			this.height = (Integer) XMLParser.extractObjectForTags(xml, "height");
		}
		
		@Override
		public StringBuffer toXML() {
			StringBuffer xml = new StringBuffer();
			XMLParser.appendObjectWithTags(xml, pwm, "pwm");
			XMLParser.appendObjectWithTags(xml, CpG, "CpG");
			XMLParser.appendObjectWithTags(xml, MpH, "MpH");
			XMLParser.appendObjectWithTags(xml, height, "height");
			return xml;
		}

		@Override
		public void generatePlot(GraphicsAdaptor ga) throws Exception {
			
			int width = SeqLogoPlotter.getWidth(height/2, pwm);
			
			//SeqLogoPlotter.plotLogo(ga.getGraphics(width, height), height, pwm);
			MLogoPlotter.plotMLogo(ga.getGraphics(width, height), 0, height, width, height, pwm, CpG, MpH, null, "", "");
			
		}
		
	}
	
	private static String[] viridis_a = {"#000004","#010107","#02020B","#030311","#050417","#07061C","#090721","#0C0926","#0F0B2C","#120D32","#150E37","#180F3E","#1C1044","#1F114A","#221150","#261257","#2A115D","#2F1163","#331068","#38106C","#3C0F71","#400F74","#451077","#491078","#4E117B","#51127C","#56147D","#5A167E","#5D177F","#611980","#661A80","#6A1C81","#6D1D81","#721F81","#762181","#792282","#7D2482","#822581","#862781","#8A2981","#8E2A81","#922B80","#962C80","#9B2E7F","#9F2F7F","#A3307E","#A7317D","#AB337C","#AF357B","#B3367A","#B83779","#BC3978","#C03A76","#C43C75","#C83E73","#CD4071","#D0416F","#D5446D","#D8456C","#DC4869","#DF4B68","#E34E65","#E65163","#E95562","#EC5860","#EE5C5E","#F1605D","#F2655C","#F4695C","#F66D5C","#F7735C","#F9785D","#F97C5D","#FA815F","#FB8661","#FC8A62","#FC9065","#FD9567","#FD9A6A","#FE9E6C","#FEA36F","#FEA873","#FEAC76","#FEB27A","#FEB67D","#FEBB81","#FEC085","#FEC488","#FEC98D","#FECD90","#FED395","#FED799","#FDDC9E","#FDE1A2","#FDE5A7","#FDEBAB","#FCEFB1","#FCF4B6","#FCF8BA","#FCFDBF"};
	
	public static void plotMLogo(Graphics2D g, int x, int yp, int w, int he, double[][] ps, double[] CpG, double[] MpH,String[] labels, String labX, String labY){
		
		int h=he/2;
		int y=yp-h;
		
		g = (Graphics2D)g.create();
		/*g.scale( 1.0/(h*100), 1.0/(h*100) );
		x *= h*100;
		y *= h*100;
		w *= h*100;
		h *= h*100;*/
		g.setColor( Color.WHITE );
		g.fillRect( x, yp-he, w, he );
		
		Font font = new Font(g.getFont().getName(),Font.PLAIN,h/14);//17
		g.setFont( font );
		
		if(labels == null){
			labels = new String[ps.length];
			for(int i=0;i<ps.length;i++){
				labels[i] = (i+1)+"";
			}
		}
		

		double wl = h*0.2;
		
		double w2 = (w-wl*1.5)/(ps.length);
		double x2 = wl*0.9;
		
		double h2 = h*0.75;
		double y2 = y - h*0.1;
		double y3 = y - h*0.8;
		
		g.setColor( Color.BLACK );
		g.setStroke( new BasicStroke( h/100f ) );
		//g.drawLine( x+(int)x2, (int)(y3+0.05*h+h), x+(int)(x2+w2*ps.length), (int)(y3+0.05*h+h) );
	//	g.drawLine( (int)(wl+w2/2d), (int)(y3+0.05*h+h), x+(int)(wl+w2*(ps.length-0.5)), (int)(y3+0.05*h+h) );
		g.drawLine( x+(int)(x2*0.94), (int)y2, x+(int)(x2*0.94), (int)(y2-h2) );
		String[] labs = {"0", "0.5", "1", "1.5", "2"};
		for(int i=0;i<=4;i++){
			g.drawLine( x+(int)(x2*0.7), (int)(y2-i*h2/4.0), x+(int)(x2*0.94), (int)(y2-i*h2/4.0) );
			Rectangle2D rect = g.getFontMetrics().getStringBounds( labs[i], g );
			g.rotate(-Math.PI/2);
			//g.drawString( labs[i], x+(int)(x2*0.6-rect.getWidth()), (int)(y2-i*h2/4.0 - rect.getCenterY()) );
			g.drawString( labs[i], -(int)(y2-i*h2/4.0 + rect.getCenterX()), x+(int)(x2*0.6) );
			g.rotate(Math.PI/2);
		}
		AffineTransform back = g.getTransform();
		g.rotate( -Math.PI/2 );
		Rectangle2D rect = g.getFontMetrics().getStringBounds( labY, g );

		g.drawString(labY,-(int)(y2-2*h2/4.0 + rect.getCenterX()), (int)(x+rect.getWidth()/2d));
		g.setTransform( back );
		x2 = wl;
		rect = g.getFontMetrics().getStringBounds( labX, g );
		g.drawString( labX, x+(int)(x2+w2*ps.length/2.0-rect.getCenterX()), (int)(y-0.3*rect.getHeight()+h) );
		for(int i=0;i<ps.length;i++){
			plotLogo( g, x+x2, y2, w2, h2, ps[i] );
			g.setColor( Color.BLACK );
			
			g.drawLine( (int)(x2+w2/2d), (int)(y3+0.1*h+h), (int)(x2+w2/2d), (int)(y3+0.1*h+0.05*h+h ) );
			
			g.rotate(-Math.PI/2);
			rect = g.getFontMetrics().getStringBounds( labels[i], g );
			//g.drawString( labels[i], x+(float)(x2+w2/2d-rect.getCenterX()), (float)(y3+1.5*rect.getHeight()+h) );
			g.drawString( labels[i], -(float)(y3+0.2*h+rect.getWidth()+h), x+(float)(x2+w2/2d - rect.getCenterY()) );
			g.rotate(Math.PI/2);
			
			x2 += w2;
		}
		
		
		y2 = y3+0.1*h;
		x2 = wl*0.9;
		
		//previous
		/*g.setColor(Color.BLACK);
		
		drawAxis(g, x, x2, y2, h2, h, w, w2, CpG, true);
				
		
		
		g.setColor(Color.MAGENTA);
		
		drawAxis(g, x, x2, y2, h2, h, w, w2, MpH, false);	*/
		
		drawSens(g, x, wl, y2, h2, h, w, w2, CpG, MpH);
		
		g.drawLine( (int)(wl+w2/2d), (int)(y3+0.1*h+h), x+(int)(wl+w2*(ps.length-0.5)), (int)(y3+0.1*h+h) );
		
	}
	
	
	private static void drawSens(Graphics2D g, int x, double wl, double y2, double h2, int h, int w, double w2, double[] CpG, double[] MpH) {
		drawRects(g,x,wl,y2,h2,h,w2,CpG, true);
		drawRects(g,x,wl,y2,h2,h,w2,MpH, false);
	}
	
	private static void drawRects(Graphics2D g, int x, double wl, double y2, double h2, int h, double w2, double[] vals, boolean CpG) {
		
		int yoff = (int)(y2+h-w2);
		if(CpG) {
			yoff = (int)(y2+h-w2 - w2);
		}
		
	//	double wl = h*0.2;
	//	g.drawLine((int)(x+wl), (int)y2, (int)(x+wl+w2*5), (int)y2+h);
		int[] xpoints1 = new int[vals.length];
		for(int i=0;i<vals.length;i++) {
			xpoints1[i] = x+(int)(wl + w2*(i+0.5));
		}
		
		double h3 = w2;
		
		double min = ToolBox.min(vals);
		if(CpG) {
			min = 0;
		}
		double max = ToolBox.max(vals);
		
		Color back = g.getColor();
		for(int i=0;i<vals.length-1;i++) {
			if(CpG) {
				float val = (float)( (max-vals[i])/(max-min)*0.6+0.3 );
				Color c = new Color(val, val, val);
				g.setColor(c);
			}else {
				double val = (vals[i]-min)/(max-min);
				int idx = (int)((val*(viridis_a.length-1)));
				Color c = Color.decode(viridis_a[idx]);
				g.setColor(c);
			}
			g.fillRect(xpoints1[i], (int)(yoff), xpoints1[i+1]-xpoints1[i], (int)h3);
		}
		g.setColor(back);
		
		
		String lab2 = "";
		if(CpG) {
			lab2 = "CpG";
		}else {
			lab2 = "MH";
		}
		
		Rectangle2D rect = g.getFontMetrics().getStringBounds( lab2, g );
		g.drawString(lab2, x+(int)(wl-rect.getCenterX()-w2/2.0), yoff+(int)(h3/2.0-rect.getCenterY()));
		
		int wstep = (int)(w2*(vals.length-3)/2.0/100.0);
		int hstep = (int)(w2/3.0*2.0);
		
		int xoff = 0;
		if(CpG) {
			xoff = xpoints1[0];
		}else {
			xoff = xpoints1[xpoints1.length-1] - (int)(wstep*100);
		}
		
		yoff = (int)( y2+ h + 1.7*h3 );
		
		
		
		int xcurr = xoff;
		
		back = g.getColor();
		for(int i=0;i<100;i++) {
			if(CpG) {
				float val = (float)( (100.0-i)/100.0*0.6+0.3 );
				Color c = new Color(val, val, val);
				g.setColor(c);
			}else {
				Color c = Color.decode(viridis_a[i]);
				g.setColor(c);
			}
			g.fillRect(xcurr, yoff, wstep, hstep);
			xcurr += wstep;
		}
		g.setColor(back);
		
		
		
		NiceScale scale = new NiceScale(min,max);
		scale.setMaxTicks(5);
		
		double niceMin = scale.getNiceMin();
		double niceMax = scale.getNiceMax();
		double ticks = scale.getTickSpacing();
		
		double pScale = (xcurr - xoff)/(max-min);
		
		String tickString = ticks+"";
		if(Math.abs(ticks)>0.0001) {
			if(tickString.matches("[0-9]+\\\\.[0-9]*[1-9]0000+[0-9]$")) {
				tickString = tickString.replaceFirst("0000+[0-9]$", "");
			}
		}
		int signif = tickString.replaceFirst("^.*\\.","").length();
		
		
		DecimalFormat dc = (DecimalFormat) DecimalFormat.getInstance(Locale.US);
		dc.setMaximumFractionDigits(signif);
		dc.setMinimumFractionDigits(signif);
		
		for(double p = niceMin;p<=niceMax;p+=ticks) {
			if(p>=min && p<=max) {
				int pos = xoff + (int)( (p-min) *pScale);
				g.drawLine(pos, yoff+hstep, pos, yoff+hstep+hstep/2);
				
				String lab = dc.format(p);//p+"";
				if(lab.matches("^-0\\.?0*$")) {
					lab = lab.substring(1);
				}
				
				rect = g.getFontMetrics().getStringBounds( lab, g );
				g.drawString(lab, (int)(pos-rect.getCenterX()), (int)(yoff+1.5*hstep+rect.getHeight()));
			}
		}
		
		g.drawLine(xoff, yoff+hstep, xcurr, yoff+hstep);
		
		String lab = "";
		if(CpG) {
			lab = "CpG";
		}else {
			lab = "Methylation sensitivity (MH)";
		}
		rect = g.getFontMetrics().getStringBounds( lab, g );
		g.drawString(lab, (int)((xoff+xcurr)/2.0 - rect.getCenterX()), (int)(yoff+2.5*hstep + rect.getHeight()));
		
		
	}
	
	
	private static void drawAxis(Graphics2D g, int x, double x2, double y2, double h2, int h, int w, double w2, double[] vals, boolean left) {
				
		double min = ToolBox.min(vals);
		double max = ToolBox.max(vals);
		NiceScale scale = new NiceScale(min,max);
		scale.setMaxTicks(4);
		
		double niceMin = scale.getNiceMin();
		double niceMax = scale.getNiceMax();
		double ticks = scale.getTickSpacing();
		if(niceMin < min) min = niceMin;
		if(niceMax > max) max = niceMax;
		
		
		int niceYoff = (int)y2+h;
		double niceYstep = h2/(max-min);
		
		/*double maxWidth = 0;
		for(double i=niceMin;i<niceMax;i+= ticks) {
			String lab = i+"";
			Rectangle2D rect = g.getFontMetrics().getStringBounds( lab, g );
			if(rect.getWidth()>maxWidth) {
				maxWidth = rect.getWidth();
			}
		}*/
		
		Integer lastEnd = null;
		for(double i=niceMin;i<=niceMax;i+= ticks ) {
			double xoff = 0;
			if(!left) {
				xoff = w2*(vals.length+0.5)+0.24*x2;
			}
			g.drawLine(x+(int)(x2*0.7+xoff), (int)(niceYoff - (i-min)*niceYstep), x+(int)(x2*0.94+xoff), (int)(niceYoff - (i-min)*niceYstep));
			
			String lab = i+"";
			if(Math.abs(ticks)>0.0001) {
				if(lab.matches("[0-9]+\\.[0-9]*[1-9]0000+[0-9]$")) {
					lab = lab.replaceFirst("0000+[0-9]$", "");
				}
			}
			
			g.rotate(-Math.PI/2);
			Rectangle2D rect = g.getFontMetrics().getStringBounds( lab, g );
			if(left) {
				//g.drawString(lab, x+(int)(x2*0.6-rect.getWidth()), (int)(niceYoff - (i-min)*niceYstep - rect.getCenterY()));
				int currStart = -(int)(niceYoff - (i-min)*niceYstep + rect.getCenterX()*1.25 );
				if(lastEnd == null || currStart>lastEnd) {
					g.drawString(lab, -(int)(niceYoff - (i-min)*niceYstep + rect.getCenterX()), x+(int)(x2*0.6));
					lastEnd = -(int)(niceYoff - (i-min)*niceYstep - rect.getCenterX()*1.25 );
				}
			}else {
				//g.drawString(lab, x+(int)(x2*1.0+maxWidth-rect.getWidth()+xoff), (int)(niceYoff - (i-min)*niceYstep - rect.getCenterY()));
				int currStart = -(int)(niceYoff - (i-min)*niceYstep + rect.getCenterX()*1.25 );
				if(lastEnd == null || currStart>lastEnd) {
					g.drawString(lab, -(int)(niceYoff - (i-min)*niceYstep + rect.getCenterX()), x+(int)(x2*1.0+rect.getHeight()+xoff));
					lastEnd = -(int)(niceYoff - (i-min)*niceYstep - rect.getCenterX()*1.25 );
				}
			}
			g.rotate(Math.PI/2);
			
		}
		
		if(!left && min <= 0 && max >= 0) {
			Graphics2D gc = (Graphics2D) g.create();
			Stroke dashed = new BasicStroke(3, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{9}, 0);
			gc.setStroke(dashed);
			gc.drawLine(x+(int)(x2*0.94), (int)(niceYoff - (0.0-min)*niceYstep) , x+(int)(x2*0.94 + w2*(vals.length+0.5)), (int)(niceYoff - (0.0-min)*niceYstep));
			gc.dispose();
		}
		
		double wl = h*0.4;
		
		int[] xpoints = new int[vals.length-1];
		int[] ypoints = new int[vals.length-1];
		for(int i=0;i<vals.length-1;i++) {
			xpoints[i] = x+(int)(wl + w2*(i+0.5));
			ypoints[i] = (int)(niceYoff - (vals[i]-min)*niceYstep);
		}
		
		
		
		if(left) {
			Graphics2D g2 = (Graphics2D) g.create();
			g2.drawLine( x+(int)(x2*0.94), (int)y2+h, x+(int)(x2*0.94), (int)(y2-h2+h) );
			
			g2.rotate( -Math.PI/2 );
			Rectangle2D rect = g2.getFontMetrics().getStringBounds( "CpG", g );
			g2.drawString("CpG",-(int)(y2-2*h2/4.0 + rect.getCenterX() + h), (int)(x+rect.getWidth()/2d));
			g2.dispose();
		}else {
			Graphics2D g2 = (Graphics2D) g.create();
			g2.drawLine( x+(int)(x2*0.94 + w2*(vals.length+0.5)), (int)y2+h, x+(int)(x2*0.94 + w2*(vals.length+0.5)), (int)(y2-h2+h) );
			
			g2.rotate( -Math.PI/2 );
			String lab = "Methylation sensitivity";
			Rectangle2D rect = g2.getFontMetrics().getStringBounds( lab, g );
			g2.drawString(lab,-(int)(y2-2*h2/4.0 + rect.getCenterX() + h), (int)(x+w-rect.getHeight()));
			g2.dispose();
		}
		
		Stroke back = g.getStroke();
		BasicStroke temp = new BasicStroke(((BasicStroke)back).getLineWidth()*1.5f,BasicStroke.CAP_BUTT,BasicStroke.JOIN_BEVEL);
		g.setStroke(temp);
		
		g.drawPolyline(xpoints, ypoints, vals.length-1);
		
		g.setStroke(back);
	}
	
}
