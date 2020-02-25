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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.parameters.AbstractSelectionParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
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
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
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
	
	private String vmargs;
	private boolean[] configThreads;
	private boolean useDiscoverDatasets;
	
	public Galaxy(String vmargs, boolean configThreads, boolean useDiscoverDatasets, JstacsTool... tools){
		this( vmargs, new boolean[]{configThreads}, useDiscoverDatasets, tools );
	}
	
	/**
	 * Creates a new Galaxy interface from a set of {@link JstacsTool}s.
	 * @param vmargs the arguments supplied to the JavaVM when running within Galaxy
	 * @param configThreads if threads should be configured
	 * @param tools the tools that should be displayed
	 * @see GalaxyAdaptor#parse(String[], boolean)
	 */
	public Galaxy(String vmargs, boolean[] configThreads, boolean useDiscoverDatasets, JstacsTool... tools){
		this.tools = tools;
		this.vmargs = vmargs;
		if(this.vmargs == null){
			this.vmargs ="";
		}else if(this.vmargs.length()!=0){
			this.vmargs = " "+this.vmargs;
		}
		if( configThreads == null || configThreads.length==1 ) {
			this.configThreads = new boolean[tools.length];
			Arrays.fill(this.configThreads, configThreads==null ? false : configThreads[0] );
		} else if( configThreads.length == tools.length ) {
			this.configThreads = configThreads.clone();
		} else {
			throw new IllegalArgumentException("Check the length of the configureThreads array.");
		}
		this.useDiscoverDatasets = useDiscoverDatasets;
		//TODO remove this.configThreads = configThreads;
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
		GalaxyAdaptor ga = new GalaxyAdaptor( tools[i], "java"+vmargs+" -jar "+jar+" "+name, "jobname", new File(jar).getParentFile().getAbsolutePath(), useDiscoverDatasets );
		ga.setHelp( tools[i].getHelpText() );
		ga.parse( args, configThreads[i] );
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
			//System.out.println("Creating " + tools.length);
			if( args.length == 1 ) 
			{
				for(int i=0;i<tools.length;i++){
					String name = tools[i].getShortName();
					//System.out.println(name);
					getGalaxyAdaptor(i, jar, vmargs, new String[]{"--create",name+".xml"} );
				}
			} else {
				String myVMArgs;
				if( args.length == 2 ) {
					myVMArgs = vmargs;
				} else {
					myVMArgs = "";
					for( int i = 2; i <args.length; i++ ) {
						myVMArgs += " " + args[i];
					}
				}
				//System.out.println(args[1] + "\t" + myVMArgs);
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
			
			/*out
			ParameterSet ps = toolParameters[idx];
			System.out.println( "Parameters of tool \""+tools[idx].getToolName()+"\" ("+tools[idx].getShortName()+", version: " + tools[idx].getToolVersion() + "):" );
			print( ps, "" );
			System.out.println( "The number of threads used for the tool, defaults to 1\t= "+ga.getThreads() );
			/**/
			ToolResult tr = tools[idx].run( (ToolParameterSet) ga.parameters, protocol, progress, ga.getThreads() );
			ResultSet[] ress = tr==null ? new ResultSet[0] : tr.getRawResult();
			if( ress.length > 0 ) {
				Pair<Result,boolean[]>[] temp = flatten(ress[0]);
				
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
			}
			ga.writeOutput();
			
		}
		
	}
	
	@Deprecated
	private void print(ParameterSet parameters, String tabPrefix){
		boolean isExp = parameters instanceof ExpandableParameterSet;
		ExpandableParameterSet exp = null;
		//TODO parameter names
		if( isExp ) {
			exp = (ExpandableParameterSet) parameters;
			parameters=(ParameterSet) parameters.getParameterAt(0).getValue();
			
			System.out.println( tabPrefix+"This parameter can be used multiple times:" );//TODO
			ParameterSet template = (ParameterSet) exp.getParameterAt(0).getValue();
			int n = exp.getNumberOfParameters();
			for(int k=0;k<n;k++){
				ParameterSet ps2 = (ParameterSet) exp.getParameterAt(k).getValue();
				for(int j=0;j<template.getNumberOfParameters();j++){
					System.out.println( tabPrefix+"\t" + (n>1?"("+(k+1)+") ":"") + ps2.getParameterAt(j).toString() );
				}
			}
		} else {		
			for(int i=0;i<parameters.getNumberOfParameters();i++){
				Parameter par = parameters.getParameterAt( i );
				if(par.getDatatype() != DataType.PARAMETERSET){
					System.out.println( tabPrefix + par.toString() );
				}else{
					if(par instanceof AbstractSelectionParameter){
						System.out.println( tabPrefix+par.toString() );
						String offset = tabPrefix;
						int l = par.getName().length();
						for( int k = 0; k<=l; k++ ) {
							offset += " ";
						}
						ParameterSet incoll = ( (AbstractSelectionParameter)par ).getParametersInCollection();
						for(int j=0;j<incoll.getNumberOfParameters();j++){
							ParameterSetContainer cont = (ParameterSetContainer)incoll.getParameterAt( j );
							if( cont.getValue().getNumberOfParameters()>0 ) {
								System.out.println( offset+"Parameters for selection \""+cont.getName()+"\":" );
								print(cont.getValue(),offset+"\t");
							} else {
								System.out.println( offset+"No parameters for selection \""+cont.getName()+"\"" );
							}
						}
					} else {
						ParameterSet ps = (ParameterSet)par.getValue();
						print(ps,tabPrefix+"\t");
					}
				}
			}
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
