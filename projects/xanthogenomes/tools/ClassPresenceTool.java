package projects.xanthogenomes.tools;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.clustering.distances.DistanceMetric;
import de.jstacs.clustering.hierachical.ClusterTree;
import de.jstacs.clustering.hierachical.Hclust;
import de.jstacs.clustering.hierachical.Hclust.Linkage;
import de.jstacs.io.FileManager;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.NiceScale;
import de.jstacs.utils.graphics.GraphicsAdaptor;
import projects.xanthogenomes.TALE;
import projects.xanthogenomes.TALEAligner;
import projects.xanthogenomes.TALEFamilyBuilder;
import projects.xanthogenomes.TALEFamilyBuilder.TALEFamily;

public class ClassPresenceTool implements JstacsTool {

	public ClassPresenceTool() {
	}

	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter builderFile = new FileParameter( "Class builder", "TALE class builder definition", "xml", true );
		FileParameter exclude = new FileParameter( "Excluded strains", "List of the names of strains to be excluded", "txt", true );
		
		FileParameter include = new FileParameter( "Included strains", "List of the names of strains to be included", "txt", true );
		
		
		
		try {
			SimpleParameter cut = new SimpleParameter(DataType.DOUBLE, "Cut height", "The height at which the class tree is cut", true,new NumberValidator<Double>(0.0, Double.MAX_VALUE),10.0);
			
			SelectionParameter sp = new SelectionParameter(DataType.PARAMETERSET, new String[]{"by class tree","alphabetically"}, new ParameterSet[]{new SimpleParameterSet(cut),new SimpleParameterSet()}, "Arrangement of classes", "If the list of classes should be arranged by the order in the class tree or alphabetically", true);

			SelectionParameter inex = new SelectionParameter(DataType.PARAMETERSET, new String[]{"none","include","exclude"}, new ParameterSet[]{new SimpleParameterSet(), new SimpleParameterSet(
					include, new SimpleParameter(DataType.BOOLEAN, "Keep order in plots", "If selected, order of strains in the plot will be kept as in the include file",true,false)
					), new SimpleParameterSet(exclude)}, "Exclude/Include", "Choose wether to exclude or (exclusively) include specific strains", true); 
			
			SimpleParameter bool = new SimpleParameter(DataType.BOOLEAN, "Strains by presence", "If selected, strains will be ordered by presence, otherwise by TALE distance", true, false);
			
			return new ToolParameterSet( getShortName(), builderFile,sp,inex,bool);
		} catch (ParameterException e) {
			return null;
		}
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		progress.setLast(1.0);
		
		progress.setCurrent(0.0);
		protocol.append("Loading class builder.\n");
		
		final TALEFamilyBuilder builder = new TALEFamilyBuilder(new StringBuffer( ((FileParameter)parameters.getParameterAt( 0 )).getFileContents().getContent()));
		progress.setCurrent(0.4);
		
		double height2 = Double.MAX_VALUE;
		if(((SelectionParameter)parameters.getParameterAt(1)).getSelected()==0){
			height2 = (Double) ((ParameterSet)parameters.getParameterAt(1).getValue()).getParameterAt(0).getValue();
		}
		
		final double height = height2;
		
		SelectionParameter inex = (SelectionParameter) parameters.getParameterAt(2);
		
		boolean bool = (Boolean) parameters.getParameterAt(3).getValue();
		String[] strainOrder = null;
		
		HashSet<String> excluded = null;
		HashSet<String> included = null;
		if(inex.getSelected() == 2){
		
			excluded = new HashSet<String>();
			if(((ParameterSet)inex.getValue()).getParameterAt(0).getValue() != null){
				String cont = ((FileParameter)((ParameterSet)inex.getValue()).getParameterAt(0)).getFileContents().getContent();
				String[] parts = cont.split("\n");
				for(int i=0;i<parts.length;i++){
					parts[i] = parts[i].trim();
					if(parts[i].length() > 0){
						excluded.add(parts[i]);
					}
				}
			}
		
		}else if(inex.getSelected() == 1){
			included = new HashSet<String>();
			if(((ParameterSet)inex.getValue()).getParameterAt(0).getValue() != null){
				String cont = ((FileParameter)((ParameterSet)inex.getValue()).getParameterAt(0)).getFileContents().getContent();
				String[] parts = cont.split("\n");
				LinkedList<String> order = new LinkedList<>();
				for(int i=0;i<parts.length;i++){
					parts[i] = parts[i].trim();
					if(parts[i].length() > 0){
						included.add(parts[i]);
						order.add(parts[i]);
					}
				}
				boolean keepOrder = (Boolean) ((SimpleParameter)((ParameterSet)inex.getValue()).getParameterAt(1)).getValue();
				if(keepOrder) {
					strainOrder = order.toArray(new String[0]);
				}
			}
		}
		
		protocol.append("Clustering classes.\n");
		ClusterTree<TALEFamily> famTree = builder.clusterFamilies();
	//	TALEFamily[] fams = builder.getFamilies();
		progress.setCurrent(0.7);
		
		ClusterTree<TALEFamily>[] subtrees = Hclust.cutTree(height, famTree);
		
		
		HashSet<String> allStrains = new HashSet<String>();
		
		TALE[] tales = builder.getAllTALEs();
		for(int i=0;i<tales.length;i++){
			if( (excluded == null && included == null) || (excluded != null && !excluded.contains(tales[i].getStrain())) || ( included != null && included.contains(tales[i].getStrain()) )){
				allStrains.add( tales[i].getStrain() );
			}
		}
		System.out.println(allStrains);
		String[] strains = allStrains.toArray(new String[0]);
		Arrays.sort(strains);
		final HashMap<String,Integer> tempMap = new HashMap<String,Integer>();
		for(int i=0;i<strains.length;i++){
			tempMap.put(strains[i], i);
		}
		
		protocol.append("Collecting class presence in strains.\n");
		final boolean[][][] present2 = new boolean[subtrees.length][][];
		final TALE[][][] presTALEs = new TALE[subtrees.length][][];
		for(int i=0;i<subtrees.length;i++){
			ClusterTree<TALEFamily> curr = subtrees[i];
			TALEFamily[] currFams = curr.getClusterElements();
			if(height==Double.MAX_VALUE){
				Arrays.sort(currFams,new Comparator<TALEFamily>(){

					@Override
					public int compare(TALEFamily o1, TALEFamily o2) {
						return o1.getFamilyId().compareTo(o2.getFamilyId());
					}
					
				});
			}
			present2[i] = new boolean[currFams.length][tempMap.size()];
			presTALEs[i] = new TALE[currFams.length][tempMap.size()];
			for(int j=0;j<currFams.length;j++){
				TALEFamily currFam = currFams[j];
				String famId = currFam.getFamilyId();
				TALE[] famTALEs = currFam.getFamilyMembers();
				
				for(int k=0;k<famTALEs.length;k++){
					String strain = famTALEs[k].getStrain();
					if(tempMap.containsKey(strain)){
						present2[ i ][ j ][ tempMap.get(strain) ] = true;
						presTALEs[ i ][ j ][ tempMap.get(strain) ] = famTALEs[k];
					}
				}
			}			
			
		}
		progress.setCurrent(0.8);
		
		
		protocol.append("Arranging strains by class presence.\n");
		DistanceMetric<String> metric = null;
		if(bool){
			metric = new DistanceMetric<String>() {

				@Override
				public double getDistance(String o1, String o2) throws Exception {
					double dist = 0;
					for(int i=0;i<present2.length;i++){
						for(int j=0;j<present2[i].length;j++){
							if(present2[i][j][tempMap.get(o1)] != present2[i][j][tempMap.get(o2)]){
								dist++;
							}
						}
					}
					return dist;
				}
				
			};
		}else{
			metric = new DistanceMetric<String>() {

				@Override
				public double getDistance(String o1, String o2) throws Exception {
					double dist = 0;
					for(int i=0;i<presTALEs.length;i++){
						for(int j=0;j<presTALEs[i].length;j++){
							double contrib = Double.POSITIVE_INFINITY;
							TALE tal1 = presTALEs[i][j][tempMap.get(o1)];
							TALE tal2 = presTALEs[i][j][tempMap.get(o2)];
							if(tal1 == null && tal2 != null){
								tal1 = findBestMatch(presTALEs,tempMap,o1,o2,tal2);
								//System.out.println("f1");
							}else if(tal1 != null && tal2 == null){
								tal2 = findBestMatch(presTALEs,tempMap,o2,o1,tal1);
								//System.out.println("f2");
							}
							if(tal1 != null && tal2 != null){
								contrib = TALEAligner.align( tal1, tal2, builder.getCosts(), builder.getAlignmentType(), builder.getExtraGapOpening(), builder.getExtraGapExtension() ).getCost();
							}else if(tal1 == null && tal2 == null){
								contrib = 0;
							}
							if(contrib > builder.getCut()*1.2){
								//System.out.println(o1+"\t"+o2+"\t"+i+","+j+"\n"+contrib+"\n"+tal1+"\n"+tal2+"\n+++++++++++++++++++++++++++++++++++++++");
								contrib = builder.getCut()*1.2;
							}
							dist += contrib;
						}
					}
					return dist;
				}

				/**
				 * Find best match for TALE from strain o2 from class (i2,j2) among the TALEs from strain o1 without a match in o2 
				 * @param presTALEs
				 * @param tempMap
				 * @param o1
				 * @param o2
				 * @param i2
				 * @param j2
				 * @return
				 */
				private TALE findBestMatch(TALE[][][] presTALEs, HashMap<String, Integer> tempMap, String o1,
						String o2, TALE origin) {
					double minDist = Double.POSITIVE_INFINITY;
					TALE bestMatch = null;
					for(int i=0;i<presTALEs.length;i++){
						for(int j=0;j<presTALEs[i].length;j++){
							int k2 = tempMap.get(o2);
							int k1 = tempMap.get(o1);
							//TALEs from o1 without a match in o2
							if(presTALEs[i][j][k1] != null && presTALEs[i][j][k2] == null){
								double temp = TALEAligner.align( origin, presTALEs[i][j][k1], builder.getCosts(), builder.getAlignmentType(), builder.getExtraGapOpening(), builder.getExtraGapExtension() ).getCost();
								if(temp < minDist){
									temp = minDist;
									bestMatch = presTALEs[i][j][k1];
								}
							}
						}
					}
					return bestMatch;
				}
				
			};
		}
		
		Hclust<String> clus = new Hclust<String>(metric, Linkage.AVERAGE);
		
		double[][] distMat = DistanceMetric.getPairwiseDistanceMatrix( metric, strains );
		
		ListResult lr2 = getList(distMat,strains);
		
		
		ClusterTree<String> strainTree = clus.cluster(distMat,strains);
		System.out.println(strainTree);
		strainTree.leafOrder(distMat);
		
		String newick = strainTree.toNewick();
		
		TextResult tr = new TextResult("Strain tree", "Strain tree (average linkage) in newick format", new FileParameter.FileRepresentation("", newick), "txt", this.getToolName(), null, true);
		
		strains = strainTree.getClusterElements();
		
		if(strainOrder != null) {
			strains = strainOrder;
		}
		
		progress.setCurrent(0.9);
		
	//	System.out.println(strainTree);
		
		int[] translate = new int[strains.length];
		for(int i=0;i<strains.length;i++){
			translate[ i ] = tempMap.get(strains[i]);
		}
		
		boolean[][][] present = new boolean[subtrees.length][][];
		for(int i=0;i<present.length;i++){
			present[i] = new boolean[present2[i].length][];
			for(int j=0;j<present[i].length;j++){
				present[i][j] = new boolean[present2[i][j].length];
				for(int k=0;k<present[i][j].length;k++){
					present[i][j][k] = present2[i][j][ translate[k] ];
				}
			}
		}
		
		protocol.append("Drawing diagram of class presence.\n");
		
		ClassPresencePlotter cpp = new ClassPresencePlotter(subtrees, present, strains, height, height<Double.MAX_VALUE);
		
		
		PlotGeneratorResult pgr = new PlotGeneratorResult("Class presence plot","",cpp,true);
		
		ListResult lr = getList(subtrees, present, strains);
		
		
		ResultSet rs = new ResultSet(new Result[]{pgr, lr, lr2, tr});
		
		progress.setCurrent(1.0);
		
		return new ToolResult("Result of "+getToolName(), "", null, rs, parameters, getToolName(), new Date(System.currentTimeMillis()));
		
	}

	private ListResult getList(double[][] distMat, String[] strains) {
		LinkedList<ResultSet> ress = new LinkedList<ResultSet>();
		
		for(int i=0;i<distMat.length;i++){
			LinkedList<Result> temp = new LinkedList<Result>();
			temp.add(new CategoricalResult("Strain", "", strains[i]));
			for(int j=0;j<distMat.length;j++){
				if(j<i){
					temp.add( new NumericalResult(strains[j], "", distMat[i][j]) );
				}else if(i < j){
					temp.add( new NumericalResult(strains[j], "", distMat[j][i]) );
				}else{
					temp.add( new NumericalResult(strains[j], "", 0.0) );
				}
			}
			ResultSet res = new ResultSet(temp);
			ress.add(res);
		}
		return new ListResult("Distance matrix", "Distance matrix of strains", null, ress.toArray(new ResultSet[0]));
	}

	private ListResult getList(ClusterTree<TALEFamily>[] subtrees, boolean[][][] present, String[] strains) {
		
		LinkedList<ResultSet> ress = new LinkedList<ResultSet>();
		
		for(int i=0;i<present.length;i++){
			TALEFamily[] fams = subtrees[i].getClusterElements();
			for(int k=0;k<present[i].length;k++){
				boolean isTrue = false;
				for(int j=0;j<present[i][k].length&&!isTrue;j++){
					isTrue |= present[i][k][j];
				}
				if(isTrue){
					String name = fams[k].getFamilyId();
					LinkedList<Result> res = new LinkedList<Result>();
					res.add(new CategoricalResult("TALE class", "", name));
					for(int j=0;j<strains.length;j++){
						res.add(new CategoricalResult(strains[j], "", present[i][k][j] ? "y": "n"));
					}
					ress.add(new ResultSet(res));
				}
			}
		}
		return new ListResult("Class presence table", "", null, ress.toArray(new ResultSet[0]));
		
	}

	@Override
	public String getToolName() {
		return "TALE Class Presence";
	}

	@Override
	public String getToolVersion() {
		return "1.1";
	}

	@Override
	public String getShortName() {
		return "presence";
	}

	@Override
	public String getDescription() {
		return "Inspect class presence in class builder";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( TALEAnalysisTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/tools/ClassPresenceTool.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

	
	
	
	public static class ClassPresencePlotter implements PlotGenerator{

		
		private ClusterTree<TALEFamily>[] subtrees;
		private boolean[][][] present;
		private String[] strains;
		private double maxDist;
		private boolean drawTrees;
		private int treeWidth = 200;
		
		private static int lineHeight = 30;
		protected static DecimalFormat format = new DecimalFormat();
		
		public ClassPresencePlotter(ClusterTree<TALEFamily>[] subtrees, boolean[][][] present, String[] strains, double maxDist, boolean drawTrees){
			this.subtrees = subtrees;
			this.present = present;
			this.strains = strains;
			this.maxDist = maxDist;
			this.drawTrees = drawTrees;
			if(!drawTrees){
				treeWidth=0;
			}
		}
		
		public ClassPresencePlotter(StringBuffer xml) throws NonParsableException {
			xml = XMLParser.extractForTag(xml, "ClassPresencePlotter");
			subtrees = (ClusterTree<TALEFamily>[]) XMLParser.extractObjectForTags(xml, "subtrees");
			present = (boolean[][][]) XMLParser.extractObjectForTags(xml, "present");
			strains = (String[]) XMLParser.extractObjectForTags(xml, "strains");
			maxDist = (Double) XMLParser.extractObjectForTags(xml, "maxDist");
			drawTrees = (Boolean) XMLParser.extractObjectForTags(xml, "drawTrees");
			treeWidth = (Integer) XMLParser.extractObjectForTags(xml, "treeWidth");
		}
		
		@Override
		public StringBuffer toXML() {
			StringBuffer xml = new StringBuffer();
			XMLParser.appendObjectWithTags(xml, subtrees, "subtrees");
			XMLParser.appendObjectWithTags(xml, present, "present");
			XMLParser.appendObjectWithTags(xml, strains, "strains");
			XMLParser.appendObjectWithTags(xml, maxDist, "maxDist");
			XMLParser.appendObjectWithTags(xml, drawTrees, "drawTrees");
			XMLParser.appendObjectWithTags(xml, treeWidth, "treeWidth");
			XMLParser.addTags(xml, "ClassPresencePlotter");
			return xml;
		}

		@Override
		public void generatePlot(GraphicsAdaptor ga) throws Exception {
			
			
			Graphics2D dummy = ga.getGraphics( 10, 10 );
			
			int[] dim = getDimension( dummy, subtrees, strains );
			
			Graphics2D graphics = ga.getGraphics( dim[0], dim[1] );
			
			graphics.setColor( Color.white );
			
			graphics.fillRect( 0, 0, dim[0], dim[1] );
			
			graphics.setColor( Color.black );
			
			graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			graphics.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
			
			graphics.setFont( new Font(Font.SANS_SERIF, Font.PLAIN, (int)(lineHeight*1.5)) );
			
			double[] fontDim = getIDFontDimensions(graphics, subtrees);
			int w = (int)(fontDim[0]*1.2);
			int h = (int)fontDim[1];
			
			
			int yoff = (int)getStrainFontDimensions(graphics, strains)[0] + 2*lineHeight;
			
			for(int i=0;i<strains.length;i++){
				//graphics.drawString(strains[i], yoff, (treeWidth+w));
				drawRotate(graphics, treeWidth+w+lineHeight*2*(i)+lineHeight*1.5+lineHeight/2, yoff-fontDim[1]/2, -50, strains[i]);
				graphics.drawLine(treeWidth+w+lineHeight*2*(i)+lineHeight/2, yoff, treeWidth+w+lineHeight*2*(i)+lineHeight/2, dim[1]-4*lineHeight);
			}
			graphics.drawLine(treeWidth+w+lineHeight*2*strains.length+lineHeight/2, yoff, treeWidth+w+lineHeight*2*strains.length+lineHeight/2, dim[1]-4*lineHeight);
			
			yoff += 2*lineHeight;
			
			for(int i=0;i<subtrees.length;i++){
				yoff = plotSubtree(graphics, subtrees[i], present[i], yoff, w, h, drawTrees);
			}
			
			//yoff += lineHeight;

			graphics.setStroke( new BasicStroke(2f) );
			
			NiceScale scale = new NiceScale( 0, maxDist );
			scale.setMaxTicks( 4 );
			double min = scale.getNiceMin();
			double max = scale.getNiceMax();
			double tick = scale.getTickSpacing();
			//System.out.println(minDist+" "+maxDist+" -> "+min+" "+max+" "+tick);
			
			double rat = (treeWidth)/(maxDist);
			
			
			int minX = lineHeight+(int)Math.round((maxDist - min)*rat);
			int maxX = minX;
			
			for(double i=min;i<=max+tick/2.0;i+=tick){
				if(i<1E-6 && i>-1E-6){
					i=0;
				}
				//System.out.println(i+" "+(i+tick)+" <= "+max+", "+(max+tick/2.0));
				int x = lineHeight+(int)Math.round((maxDist - i)*rat);
				if(x > lineHeight/2){
					graphics.drawLine( x, yoff, x, yoff+lineHeight/2 );
					maxX = x;
					String label = format.format( i );
					int textCenter = (int)Math.round( graphics.getFontMetrics().getStringBounds( label, graphics ).getCenterX() );
					graphics.drawString( label, x-textCenter, yoff+2*lineHeight-5 );
				}
			}
			
			graphics.drawLine(minX, yoff, maxX, yoff);
			
			
		}
		
		public static void drawRotate(Graphics2D g2d, double x, double y, int angle, String text) 
		{    
		    g2d.translate((float)x,(float)y);
		    g2d.rotate(Math.toRadians(angle));
		    g2d.drawString(text,0,0);
		    g2d.rotate(-Math.toRadians(angle));
		    g2d.translate(-(float)x,-(float)y);
		}  
		
		public int[] getDimension(Graphics2D graphics, ClusterTree<TALEFamily>[] trees, String[] strains){
			graphics = (Graphics2D)graphics.create();
			
			graphics.setFont( new Font(Font.SANS_SERIF, Font.PLAIN, (int)(lineHeight*1.5)) );
			
			double[] idDim = getIDFontDimensions( graphics, trees );
			
			int fh = (int)getStrainFontDimensions(graphics, strains)[0];
			int height = fh;
			for(int i=0;i<trees.length;i++){
				//height += trees[i].getNumberOfElements()*lineHeight*2;
				TALEFamily[] fams = trees[i].getClusterElements();
				for(int j=0;j<fams.length;j++){
					boolean included = false;
					TALE[] tales = fams[j].getFamilyMembers();
					outerloop:
					for(int k=0;k<tales.length;k++){
						for(int l=0;l<strains.length;l++){
							if(strains[l].equals( tales[k].getStrain() )){
								included = true;
								break outerloop;
							}
						}
					}
					if(included){
						height += lineHeight*2;
					}
					
				}
			}
			height += lineHeight*2*2;
			
			height += lineHeight*2;
			
			int tw = treeWidth + (int)(idDim[0]*1.2) + (strains.length+1)*lineHeight*2 + (int)(fh*0.6);
			
			
			return new int[]{tw, height};
		}
		
		private double[] getIDFontDimensions(Graphics2D graphics, ClusterTree<TALEFamily>[] trees){
			double fontHeight = 0;
			double fontWidth = 0;
			
			for(int j=0;j<trees.length;j++){
				TALEFamily[] members = trees[j].getClusterElements();
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
			}
			return new double[]{fontWidth,fontHeight};
		}
		
		private double[] getStrainFontDimensions(Graphics2D graphics, String[] strains){
			double fontHeight = 0;
			double fontWidth = 0;

			for(int i=0;i<strains.length;i++){
				String id = strains[i];
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
		
		private int plotSubtree(Graphics2D graphics, ClusterTree<TALEFamily> subtree, boolean[][] present, int yoff, int idWidth, int idHeight, boolean drawTree){
			
			graphics = (Graphics2D) graphics.create();
			
			int xoff = treeWidth+lineHeight/2; 
			
			int yoffBack = yoff;
			
			TALEFamily[] fams = subtree.getClusterElements();
			xoff += idWidth;
			
			yoff = yoffBack;
			
			for(int i=0;i<present.length;i++){
				boolean isTrue = false;
				for(int j=0;j<present[i].length&&!isTrue;j++){
					isTrue |= present[i][j];
				}
				if(isTrue){
					graphics.drawString( fams[i].getFamilyId(), xoff-idWidth+(int)(lineHeight*0.75), yoff- (lineHeight*2-idHeight) );
					graphics.drawLine(xoff, yoff-2*lineHeight, xoff+present[i].length*lineHeight*2, yoff-2*lineHeight);
					for(int j=0;j<present[i].length;j++){
						if(present[i][j]){
							graphics.setColor(new Color(0f, 0f, 0.5f, 0.6f));
							graphics.fillRect(xoff+j*lineHeight*2+2, yoff-2*lineHeight+2, 2*lineHeight-4, 2*lineHeight-4);
							graphics.setColor(Color.BLACK);
						}
					}
					yoff += 2*lineHeight;
				}
			}
			graphics.drawLine(xoff, yoff-2*lineHeight, xoff+present[0].length*lineHeight*2, yoff-2*lineHeight);
			
			graphics.setStroke( new BasicStroke(2f) );
			if(drawTree){
				drawTree(graphics,subtree,present,0,yoffBack-2*lineHeight,0);
			}
			
			return yoff;
			
		}

		private int drawTree(Graphics2D graphics, ClusterTree<TALEFamily> subtree, boolean[][] present, int off, int yoff, int level) {
			
			int n=0;
			for(int i=off;i<off+subtree.getNumberOfElements();i++){
				boolean isTrue = false;
				for(int j=0;j<present[i].length&&!isTrue;j++){
					isTrue |= present[i][j];
				}
				if(isTrue){
					n++;
				}
			}
			
			int xoff = lineHeight + (int) ( (1.0-subtree.getDistance()/maxDist)*treeWidth );
			
			if(n>1){
				
				ClusterTree<TALEFamily> left = subtree.getSubTrees()[0];
				ClusterTree<TALEFamily> right = subtree.getSubTrees()[1];
				
				
				int numLeft = drawTree(graphics, left,present,off,yoff,level+1);
				int numRight = drawTree(graphics, right,present,off+left.getNumberOfElements(),yoff + numLeft*2*lineHeight,level+1);
				
				
				graphics.drawLine(xoff, yoff+numLeft*lineHeight, xoff, yoff+numLeft*2*lineHeight + numRight*lineHeight);
				
				int xoff2 = treeWidth+lineHeight;
				if(numLeft>1){
					xoff2 = lineHeight + (int) ( (1.0-left.getDistance()/maxDist)*treeWidth ); 
				}
				graphics.drawLine(xoff, yoff+numLeft*lineHeight, xoff2 , yoff+numLeft*lineHeight);
				//System.out.println("1: "+xoff+" "+xoff2+" "+(yoff+numLeft*lineHeight));
				
				xoff2 = treeWidth+lineHeight;
				if(numRight>1){
					xoff2 = lineHeight + (int) ( (1.0-right.getDistance()/maxDist)*treeWidth ); 
				}
				graphics.drawLine(xoff, yoff+numLeft*2*lineHeight + numRight*lineHeight, xoff2 , yoff+numLeft*2*lineHeight + numRight*lineHeight);
				//System.out.println("2: "+xoff+" "+xoff2+" "+(yoff+numLeft*2*lineHeight + numRight*lineHeight));
				
			}else if(n>0 && level>0){
				//int xoff = (int) ( (1.0-subtree.getDistance()/maxDist)*treeWidth );
				//graphics.drawLine(prevXoff, yoff+lineHeight, treeWidth-lineHeight/5, yoff+lineHeight);
			}
			
			
			return n;
		}
		
		
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
		return new String[] {"@article{grau16annotale,\n" + 
				"	title = {{AnnoTALE}: bioinformatics tools for identification, annotation, and nomenclature of {TALEs} from \\emph{Xanthomonas} genomic sequences},\n" + 
				"	author = {Grau, Jan and Reschke, Maik and Erkes, Annett and Streubel, Jana and Morgan, Richard D. and Wilson, Geoffrey G. and Koebnik, Ralf and Boch, Jens},\n" + 
				"	journal = {Scientific Reports},\n" + 
				"	year = {2016},\n" + 
				"	volume = {6},\n" + 
				"	pages = {21077},\n" + 
				"	doi = {10.1038/srep21077}\n" + 
				"	}\n",
				"@article{erkes17evolution,\n" + 
				"	title = {Evolution of Transcription Activator-Like Effectors in \\emph{Xanthomonas oryzae}},\n" + 
				"	author = {Erkes, Annett and Reschke, Maik and Boch, Jens and Grau, Jan},\n" + 
				"	journal = {Genome Biology and Evolution},\n" + 
				"	year = {2017},\n" + 
				"	volume = {9},\n" + 
				"	number = {6},\n" + 
				"	pages = {1599-1615},\n" + 
				"	doi = {10.1093/gbe/evx108}\n" + 
				"	}\n"};
	}

	
	
}
