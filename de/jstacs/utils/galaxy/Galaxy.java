package de.jstacs.utils.galaxy;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.cli.CLI;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.DataSetResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.utils.Pair;
import de.jstacs.utils.galaxy.GalaxyAdaptor.FileResult;
import de.jstacs.utils.galaxy.GalaxyAdaptor.LinkedImageResult;
import de.jstacs.utils.galaxy.GalaxyAdaptor.Protocol;
import de.jstacs.utils.graphics.GraphicsAdaptor;
import de.jstacs.utils.graphics.GraphicsAdaptorFactory;
import de.jstacs.utils.graphics.PDFAdaptor;
import de.jstacs.utils.graphics.RasterizedAdaptor;
import de.jstacs.utils.graphics.GraphicsAdaptorFactory.OutputFormat;


public class Galaxy {

	private JstacsTool[] tools;
	private ParameterSet[] toolParameters;
	private boolean[][] addLine;
	private String vmargs;
	private boolean configThreads;
	
	public Galaxy(String vmargs, boolean configThreads, JstacsTool... tools){
		this.tools = tools;
		this.vmargs = vmargs;
		if(this.vmargs == null){
			this.vmargs ="";
		}else if(this.vmargs.length()!=0){
			this.vmargs = " "+this.vmargs;
		}
		this.configThreads = configThreads;
		toolParameters = new ParameterSet[tools.length];
		this.addLine = new boolean[tools.length][];
		for(int i=0;i<tools.length;i++){
			toolParameters[i] = tools[i].getToolParameters();
			addLine[i] = new boolean[toolParameters[i].getNumberOfParameters()];
		}
	}
	
	private int getToolIndex(String shortName){
		for(int i=0;i<tools.length;i++){
			if(shortName.equals( tools[i].getShortName() )){
				return i;
			}
		}
		return -1;
	}
	
	public void run(String[] args) throws Exception{
		
		File jarfile = new File(Galaxy.class.getProtectionDomain().getCodeSource().getLocation().toURI().getPath());
		String jar = jarfile.getAbsolutePath();
		
		if("--create".equals( args[0]) ){
			
			for(int i=0;i<tools.length;i++){
				String name = tools[i].getShortName();
				
				String command = "java"+vmargs+" -jar "+jar+" "+name;
				
				GalaxyAdaptor ga = new GalaxyAdaptor(toolParameters[i], addLine[i], tools[i].getToolName(), tools[i].getDescription(), "1.0", command, "jobname");
				
				ga.setHelp( tools[i].getHelpText() );
				
				ga.parse( new String[]{"--create",name+".xml"}, configThreads );
				
			}
			
		}else{
			
			String toolname = args[0];
			int idx = getToolIndex( toolname );
			
			String[] args2 = new String[args.length-1];
			System.arraycopy( args, 1, args2, 0, args2.length );
			
			String name = tools[idx].getShortName();
			String command = "java"+vmargs+" -jar "+jar+" "+name;
			
			GalaxyAdaptor ga = new GalaxyAdaptor(toolParameters[idx], addLine[idx], tools[idx].getToolName(), tools[idx].getDescription(), "1.0", command, "jobname");
			
			ga.setHelp( tools[idx].getHelpText() );
			
			ga.parse( args2, configThreads );
			
			Protocol protocol = ga.getProtocol( false );
			
			ProgressUpdater progress = new ProgressUpdater();
			
			ResultSet ress = tools[idx].run( toolParameters[idx], protocol, progress ).getRawResult()[0];			
			
			Pair<Result,boolean[]>[] temp = flatten(ress);
			
			HashSet<String> names = new HashSet<String>();
			
			for(int i=0;i<temp.length;i++){
				
				Result res = temp[i].getFirstElement();
				boolean export = temp[i].getSecondElement()[0];
				boolean includeInSummary = temp[i].getSecondElement()[1];
				
				String exportExtension = null;
				if(res instanceof TextResult){
					String mime = ( (TextResult)res ).getMime();
					if(mime != null){
						String[] exts = mime.split( "\\," );
						exportExtension = exts[0];
					}
				}
				if(res instanceof PlotGeneratorResult){
					
					PlotGenerator pg = ((PlotGeneratorResult)res).getValue();
					RasterizedAdaptor rast = new RasterizedAdaptor( "png" );
					pg.generatePlot( rast );
					BufferedImage bi = rast.getImage();
					
					PDFAdaptor pdf = new PDFAdaptor();
					pg.generatePlot( pdf );
					
					String filename = res.getName().replaceAll( "[\\s\\:\\/]", "_" );
					String temp2 = filename;
					int j=1;
					while(names.contains( temp2 )){
						temp2 = filename+"_"+j;
						j++;
					}
					filename = temp2;
					names.add(filename);
					
					filename = "./"+filename+".pdf";
					
					pdf.generateOutput( filename );
					
					FileResult link = new FileResult( res.getName(), "PDF", filename );
					
					LinkedImageResult lir = new LinkedImageResult( res.getName(), res.getComment(), bi, link );
					
					res = lir;
					
				}
				
				ga.addResult( res, export, includeInSummary, exportExtension );
				
			}
			
			ga.writeOutput();
			
		}
		
	}
	
	
	private static Pair<Result,boolean[]>[] flatten(ResultSet ress){
		LinkedList<Pair<Result,boolean[]>> all = new LinkedList<Pair<Result,boolean[]>>();
		
		for(int i=0;i<ress.getNumberOfResults();i++){
			Result res = ress.getResultAt( i );
			
			boolean export = res instanceof FileResult || res instanceof TextResult || res instanceof DataSetResult;
			if(res instanceof TextResult){
				export = ( (TextResult)res ).getExport();
			}
			boolean includeInSummary = !export;
			
			//TODO PlotGeneratorResult
			if(res instanceof ResultSetResult){
				Pair<Result,boolean[]>[] temp = flatten( ( (ResultSetResult)res ).getRawResult()[0]);
				all.add(new Pair<Result, boolean[]>(new GalaxyAdaptor.HeadResult(res.getName(), res.getComment()), new boolean[]{false,true}));
				Collections.addAll( all, temp );
			}else{
				all.add( new Pair<Result, boolean[]>( res, new boolean[]{export,includeInSummary} ) );
			}
		}
		return all.toArray( new Pair[0] );
		
	}
	
	
	
	
}