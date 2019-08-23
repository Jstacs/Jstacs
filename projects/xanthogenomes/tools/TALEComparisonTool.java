package projects.xanthogenomes.tools;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.algorithms.alignment.Alignment;
import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.PairwiseStringAlignment;
import de.jstacs.algorithms.alignment.cost.SimpleCosts;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.FileManager;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.graphics.GraphicsAdaptor;
import projects.xanthogenomes.TALE;
import projects.xanthogenomes.TALE.Repeat;
import projects.xanthogenomes.TALE.Type;

public class TALEComparisonTool implements JstacsTool {

	
	private static class RepeatPlotGenerator implements PlotGenerator{

		private TALE t1;
		private TALE t2;
		double[][] dist;
		
		
		public RepeatPlotGenerator(StringBuffer xml) throws NonParsableException{
			xml = XMLParser.extractForTag(xml, "RPG");
			t1 = (TALE) XMLParser.extractObjectForTags(xml, "t1");
			t2 = (TALE) XMLParser.extractObjectForTags(xml, "t2");
			dist = (double[][]) XMLParser.extractObjectForTags(xml, "dist");
		}
		
		public RepeatPlotGenerator(TALE t1, TALE t2, double[][] dist){
			this.t1 = t1;
			this.t2 = t2;
			this.dist = dist;
		}
		
		@Override
		public StringBuffer toXML() {
			StringBuffer xml = new StringBuffer();
			XMLParser.appendObjectWithTags(xml, t1, "t1");
			XMLParser.appendObjectWithTags(xml, t2, "t2");
			XMLParser.appendObjectWithTags(xml, dist, "dist");
			XMLParser.addTags(xml, "RPG");
			return xml;
		}

		@Override
		public void generatePlot(GraphicsAdaptor ga) throws Exception {
			
			int margMain = 130;
			int margLeft = 70;
			int margTop = 70;
			int repHeight = 30;
			
			int height = margMain+margTop+(t1.getNumberOfRepeats()+1)*repHeight;
			int width = margLeft+(t2.getNumberOfRepeats()+1)*repHeight;
			
			Graphics2D g = ga.getGraphics(width,height);
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, width, height);
			g.setColor(Color.BLACK);
			
			int offx = margLeft;
			int offy = margMain+margTop;
			
			double maxDist = dist[dist.length-1][ dist[dist.length-1].length-1 ];
			for(int i=0;i<dist.length-1;i++){
				for(int j=0;j<dist[i].length-1;j++){
					if(dist[i][j] > maxDist){
						maxDist = dist[i][j];
					}
				}
			}
			
			
			for(int i=0;i<dist.length;i++){
				for(int j=0;j<dist[i].length;j++){
					if( (i<dist.length-1 && j<dist[i].length-1) || (i==dist.length-1 && j==dist[i].length-1) ){
						Color c = Color.WHITE;
						if(dist[i][j] > 0){
							float mix = (float)(dist[i][j]/maxDist);
							c = new Color(Color.YELLOW.getRed()/255f*(1f-mix) + 1f*(mix), Color.YELLOW.getGreen()/255f*(1f-mix), Color.YELLOW.getBlue()/255f*(1f-mix) );
						}
						g.setColor(c);
						g.fillRect(offx+j*repHeight, offy+i*repHeight, repHeight, repHeight);
						g.setColor(Color.BLACK);
					}else{
						Color c = new Color(0.9f, 0.9f, 0.9f);
						g.setColor(c);
						g.fillRect(offx+j*repHeight, offy+i*repHeight, repHeight, repHeight);
						g.setColor(Color.BLACK);
					}
				}
			}
			
			
			g.setFont( new Font(Font.SANS_SERIF, Font.PLAIN, repHeight/2) );
			
			for(int i=0;i<t1.getNumberOfRepeats();i++){
				String str = (i+1)+"";
				Rectangle2D dim = g.getFontMetrics().getStringBounds( str, g );
				int x = margLeft - margLeft/10 - (int)dim.getWidth();
				int y = offy+i*repHeight + repHeight - (repHeight - (int)dim.getHeight())/2;
				g.drawString(str, x, y);
				
				str = t1.getRepeat(i).getRvd();
				if(i<t1.getNumberOfRepeats()-1){
					if(t1.getRepeat(i).getType() == Type.LONG){
						g.setColor(Color.RED);
					}else if(t1.getRepeat(i).getType() == Type.SHORT){
						g.setColor(Color.GREEN);
					}else if(t1.getRepeat(i).getType() == Type.THIRTYFIVE){
						g.setColor(Color.GRAY);
					}
				}
				dim = g.getFontMetrics().getStringBounds( str, g );
				x = margLeft - margLeft/2 - (int)dim.getWidth();
				y = offy+i*repHeight + repHeight - (repHeight - (int)dim.getHeight())/2;
				g.drawString(str, x, y);
				g.setColor(Color.BLACK);
				
			}
			for(int i=0;i<t2.getNumberOfRepeats();i++){
				String str = (i+1)+"";
				Rectangle2D dim = g.getFontMetrics().getStringBounds( str, g );
				g.drawString(str, offx+i*repHeight + (repHeight - (int)dim.getWidth())/2, margMain+margTop - margTop/10);
				
				str = t2.getRepeat(i).getRvd();
				if(i<t2.getNumberOfRepeats()-1){
					if(t2.getRepeat(i).getType() == Type.LONG){
						g.setColor(Color.RED);
					}else if(t2.getRepeat(i).getType() == Type.SHORT){
						g.setColor(Color.GREEN);
					}else if(t2.getRepeat(i).getType() == Type.THIRTYFIVE){
						g.setColor(Color.GRAY);
					}
				}
				dim = g.getFontMetrics().getStringBounds( str, g );
				g.drawString(str, offx+i*repHeight + (repHeight - (int)dim.getWidth())/2, margMain+margTop - margTop/2);
				g.setColor(Color.BLACK);
			}
			
			String[] temp = null;
			if(t1.getId().equals(t2.getId())){
				temp = new String[]{"Difference of repeats of",t1.getId()};
			}else{
				temp = new String[]{"Difference of repeats of",t1.getId()+" and",t2.getId()};
			}
			
			double sumHeight = 0.0;
			double maxWidth = 0.0;
			for(int i=0;i<temp.length;i++){
				Rectangle2D dim = g.getFontMetrics().getStringBounds( temp[i], g );
				if(dim.getWidth() > maxWidth){
					maxWidth = dim.getWidth();
				}
				sumHeight += dim.getHeight()*1.5;
			}
			
			
			
			int step = (int)((margMain-repHeight)/maxDist);
			
			for(int i=1;i<maxDist;i++){
				
				float mix = (float)(i/maxDist);
				Color c = new Color(Color.YELLOW.getRed()/255f*(1f-mix) + 1f*(mix), Color.YELLOW.getGreen()/255f*(1f-mix), Color.YELLOW.getBlue()/255f*(1f-mix) );
				g.setColor(c);
				g.fillRect(margLeft, repHeight/2+(int)((maxDist-i-1)*step), repHeight, step);
				g.setColor(Color.BLACK);
			}
			g.drawLine((int)(margLeft*0.9), repHeight/2+step/2, margLeft-1, repHeight/2+step/2);
			Rectangle2D dim = g.getFontMetrics().getStringBounds(((int)maxDist)+"", g );
			g.drawString(((int)maxDist)+"", (int)(margLeft*0.8-dim.getWidth()), (int)(repHeight/2+step/2+dim.getHeight()/2) );
			
			g.drawLine((int)(margLeft*0.9), (int)(repHeight/2+maxDist*step-step/2), margLeft-1, (int)(repHeight/2+maxDist*step-step/2));
			dim = g.getFontMetrics().getStringBounds( "0", g );
			g.drawString("0", (int)(margLeft*0.8-dim.getWidth()), (int)(repHeight/2+maxDist*step-step/2+dim.getHeight()/2)  );
			
			g.drawLine(margLeft-1, repHeight/2+step/2, margLeft-1, (int)(repHeight/2+maxDist*step-step/2));
			
			double scaleHeight = margMain/sumHeight;
			
			double scaleWidth = (width-margLeft-2*repHeight)/maxWidth;
			
			double scale = Math.min(2, Math.min(scaleHeight, scaleWidth) );
			
			g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, (int)(g.getFont().getSize()*scale)) );
			
			for(int i=0;i<temp.length;i++){
				dim = g.getFontMetrics().getStringBounds( temp[i], g );
				g.drawString(temp[i], margLeft+repHeight+(float)(width-dim.getWidth()-margLeft-2*repHeight)/2f, (float)(sumHeight*scale*0.9/temp.length*(i+1)) );
			}
			
			
			
			
		}
		
	}
	
	
	@Override
	public ToolParameterSet getToolParameters() {
		
		
		FileParameter tales = new FileParameter( "TALE sequences", "TALE sequences, either as complete DNA or AS sequences (e.g., output of TALE Prediction).", "fasta,fa,fas",true);
		
		return new ToolParameterSet( getShortName(),  tales );
		
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		
		TALE[] input = ClassBuilderTool.readProteinTALEs( ( (FileParameter)parameters.getParameterAt( 0 ) ).getFileContents() , protocol );
		
		TALE[] tales = null;
		
		if(input.length == 1){
			tales = new TALE[]{input[0],input[0]};
		}else{
			tales = input;
		}
		
		
		LinkedList<Result> ress = new LinkedList<Result>();
		
		for(int i=0;i<tales.length;i++){
			for(int j=0;j<=i;j++){
				TALE tal1 = tales[i];
				TALE tal2 = tales[j];
				
				PlotGeneratorResult pgr = compareTALEs(tal1,tal2);
				
				ress.add(pgr);
				
			}
		}
		
		ResultSet set = new ResultSet(ress.toArray(new Result[0]));
		
		return new ToolResult("Result of "+getToolName(), getToolName()+" on \""+((FileParameter)parameters.getParameterAt(0)).getFileContents().getFilename()+"\"", null, set, parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	private PlotGeneratorResult compareTALEs(TALE tal1, TALE tal2) {
		
		Alignment al = new Alignment(new SimpleCosts(0, 1, 1)); 
		
		TALE to1 = tal1;
		TALE to2 = tal2;
		
		if(tal1.getDnaOriginal() != null && tal2.getDnaOriginal() != null){
			tal1 = tal1.getDnaOriginal();
			tal2 = tal2.getDnaOriginal();
			
		}

		double[][] dist = new double[tal1.getNumberOfRepeats()][tal2.getNumberOfRepeats()];
		for(int i=0;i<dist.length;i++){
			for(int j=0;j<dist[i].length;j++){
				Repeat r1 = tal1.getRepeat(i);
				Repeat r2 = tal2.getRepeat(j);

				Sequence s1 = r1.getRepeat();
				Sequence s2 = r2.getRepeat();

				dist[i][j] = al.getAlignment(AlignmentType.GLOBAL, s1, s2).getCost();
			}
		}
		
		
		PlotGenerator gen = new RepeatPlotGenerator(to1, to2, dist);
		
		return new PlotGeneratorResult("Repeat difference of "+tal1.getId() + ( tal1.getId().equals(tal2.getId()) ? "" : " and "+tal2.getId() ), "Pairwise edit distance of TALE repeats", gen, true);
		
	}

	@Override
	public String getToolName() {
		return "TALE Repeat Differences";
	}

	@Override
	public String getToolVersion() {
		return "1.0";
	}

	@Override
	public String getShortName() {
		return "repdiff";
	}

	@Override
	public String getDescription() {
		return "Compares TALE repeats on the level of DNA or AA sequences";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( TALEAnalysisTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/tools/TALEComparisonTool.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}


	@Override
	public ToolResult[] getTestCases(String path) {
		return null;
	}

	@Override
	public void clear() {		
	}

	@Override
	public String[] getReferences() {
		return null;
	}
	
	
}
