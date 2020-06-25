package projects.methyl;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;

import de.jstacs.classifiers.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.utils.graphics.GraphicsAdaptor;

public class CurvePlotter implements PlotGenerator {

	private DoubleTableResult dtr;
	private String xlab;
	private String ylab;
	
	CurvePlotter(DoubleTableResult dtr, String xlab, String ylab){
		this.dtr = dtr;
		this.xlab = xlab;
		this.ylab = ylab;
	}
	
	public CurvePlotter(StringBuffer xml) throws NonParsableException {
		fromXML(xml);
	}

	public void fromXML(StringBuffer xml) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, "CurvePlotter");
		dtr = (DoubleTableResult) XMLParser.extractObjectForTags(xml, "dtr");
		xlab = (String) XMLParser.extractObjectForTags(xml, "xlab");
		ylab = (String) XMLParser.extractObjectForTags(xml, "ylab");
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, dtr, "dtr");
		XMLParser.appendObjectWithTags(xml, xlab, "xlab");
		XMLParser.appendObjectWithTags(xml, ylab, "ylab");
		XMLParser.addTags(xml, "CurvePlotter");
		return xml;
	}
	
	

	@Override
	public void generatePlot(GraphicsAdaptor ga) throws Exception {
		
		plot(ga.getGraphics(1180, 1180));
		
		

	}

	private void plot(Graphics2D graphics) {
		graphics = (Graphics2D) graphics.create();
		graphics.setColor(Color.WHITE);
		graphics.fillRect(0, 0, 1180, 1180);
		graphics.setColor(Color.BLACK);
		
		graphics.setStroke(new BasicStroke(5f));
		Font font = new Font(graphics.getFont().getName(),Font.PLAIN,30);//17
		graphics.setFont( font );
		
		graphics.translate(150, 1030);
		
		graphics.drawLine(0, 15, 1000, 15);
		graphics.drawLine(-15,0,-15,-1000);
		
		String[] labs = new String[] {"0","0.25","0.5","0.75","1"};
		
		for(int i=0;i<labs.length;i++) {
			int y = - i*250;
			graphics.drawLine(-40, y, -15, y);
			Rectangle2D rect = graphics.getFontMetrics().getStringBounds(labs[i], graphics);
			
			graphics.drawString(labs[i], -45-(int)rect.getWidth(), y-(int)rect.getCenterY());
			
			int x = i*250;
			graphics.drawLine(x, 40, x, 15);
			graphics.drawString(labs[i], x-(int)rect.getCenterX(), 45+(int)rect.getHeight());
		}
		
		Rectangle2D rect = graphics.getFontMetrics().getStringBounds(xlab, graphics);
		graphics.drawString(xlab,500-(int)rect.getCenterX(),100+(int)rect.getHeight());
		
		rect = graphics.getFontMetrics().getStringBounds(ylab, graphics);
		graphics.rotate(-Math.PI/2.0);
		graphics.drawString(ylab,(500-(int)rect.getCenterX()),-110);
		graphics.rotate(Math.PI/2.0);
		
		
		int[] xPoints = new int[dtr.getNumberOfLines()];
		int[] yPoints = new int[dtr.getNumberOfLines()];
		
		for(int i=0;i<dtr.getNumberOfLines();i++) {
			double[] temp = dtr.getLine(i);
			xPoints[i] = (int)(temp[0]*1000);
			yPoints[i] = -(int)(temp[1]*1000);
		}
		
		graphics.setColor(new Color(0.45f,0.65f,0.96f));
		
		graphics.drawPolyline(xPoints, yPoints, xPoints.length);
		
		
	}

}
