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

package de.jstacs.tools.ui.galaxy;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.parameters.ParameterSet;
import de.jstacs.results.DataSetResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;
import de.jstacs.results.StorableResult;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.JstacsTool.ResultEntry;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor.FileResult;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor.LinkedImageResult;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor.Protocol;
import de.jstacs.utils.Pair;
import de.jstacs.utils.graphics.PDFAdaptor;
import de.jstacs.utils.graphics.RasterizedAdaptor;

/**
 * Class for generating a generic Galaxy interface for a set of {@link JstacsTool}s.
 * Internally, it uses much of the infrastructure provided by {@link GalaxyAdaptor}.
 * 
 * This class belongs to the trinity of {@link Galaxy}, {@link CLI}, and JavaFX interface that
 * may all be created from a set of {@link JstacsTool}s.
 * 
 * @author Jan Grau
 *
 */
public class Galaxy {

	private JstacsTool[] tools;
	private ParameterSet[] toolParameters;
	private ResultEntry[][] defaultResults;
	private boolean[][] addLine;
	private String vmargs;
	private boolean configThreads;
	
	/**
	 * Creates a new Galaxy interface from a set of {@link JstacsTool}s.
	 * @param vmargs the arguments supplied to the JavaVM when running within Galaxy
	 * @param configThreads if threads should be configured
	 * @param tools the tools that should be displayed
	 * @see GalaxyAdaptor#parse(String[], boolean)
	 */
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
		this.defaultResults = new ResultEntry[tools.length][];
		for(int i=0;i<tools.length;i++){
			toolParameters[i] = tools[i].getToolParameters();
			this.defaultResults[i] = tools[i].getDefaultResultInfos();
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
	
	private GalaxyAdaptor getGalaxyAdaptor( int i, String jar, String vmargs, String[] args ) throws Exception {
		String name = tools[i].getShortName();
		
		String command = "java"+vmargs+" -jar "+jar+" "+name;
		
		GalaxyAdaptor ga = new GalaxyAdaptor(toolParameters[i], defaultResults[i], addLine[i], tools[i].getToolName(), tools[i].getDescription(), tools[i].getToolVersion(), command, "jobname");
		
		ga.setHelp( tools[i].getHelpText() );
		
		ga.parse( args, configThreads );
		return ga;
	}
	
	/**
	 * Runs this Galaxy interface with the supplied arguments.
	 * If the first argument equals "--create", the Galaxy config files for the tools are created.
	 * Otherwise, the tools selected by the supplied arguments is run with these arguments.
	 * @param args the arguments
	 * @throws Exception if the configuration could not be created or the tool threw an exception
	 */
	public void run(String[] args) throws Exception{
		
		File jarfile = new File(Galaxy.class.getProtectionDomain().getCodeSource().getLocation().toURI().getPath());
		String jar = jarfile.getAbsolutePath();
		
		if("--create".equals( args[0]) ){
			
			if( args.length == 1 ) 
			{
				for(int i=0;i<tools.length;i++){
					String name = tools[i].getShortName();
					getGalaxyAdaptor(i, jar, vmargs, new String[]{"--create",name+".xml"} );
				}
			} else {
				String myVMArgs;
				if( args.length == 2 ) {
					myVMArgs = vmargs;
				} else {
					myVMArgs = "";
					for( int i = 3; i <args.length; i++ ) {
						myVMArgs += " " + args[i];
					}
				}
				getGalaxyAdaptor(getToolIndex(args[1]), jar, myVMArgs, new String[]{"--create",args[1]+".xml"} );
			}
		}else{
			
			String toolname = args[0];
			int idx = getToolIndex( toolname );
			
			String[] args2 = new String[args.length-1];
			System.arraycopy( args, 1, args2, 0, args2.length );
			
			GalaxyAdaptor ga = getGalaxyAdaptor(idx, jar, "", args2);//do not use vmargs here
			
			Protocol protocol = ga.getProtocol( false );
			
			ProgressUpdater progress = new ProgressUpdater();
			
			ResultSet ress = tools[idx].run( toolParameters[idx], protocol, progress, ga.getThreads() ).getRawResult()[0];			
			
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
			
			boolean export = res instanceof FileResult || res instanceof DataSetResult || res instanceof StorableResult;
			if(res instanceof TextResult){
				export = ( (TextResult)res ).getExport();
			}
			if(res instanceof ListResult){
				export = ( (ListResult)res ).getExport();
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
