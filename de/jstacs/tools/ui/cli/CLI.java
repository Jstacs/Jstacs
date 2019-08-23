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

package de.jstacs.tools.ui.cli;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.DataType;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.AbstractSelectionParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.MultiSelectionParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.results.savers.ResultSaver;
import de.jstacs.results.savers.ResultSaverLibrary;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SafeOutputStream;

/**
 * Class that allows for building generic command line interface (CLI) applications based on the {@link JstacsTool} interface.
 * The generic description of command line parameters is automatically inferred from the corresponding {@link Parameter}s and {@link ParameterSet}s of
 * the {@link JstacsTool}s.
 * 
 * If this class is used for a single {@link JstacsTool}, the parameters may be supplied directly using a <code>key=value</code> interface. In case of multiple tools, each specific tool may be addressed by its short name ({@link JstacsTool#getShortName()}).
 * 
 * @author Jan Grau, Jens Keilwagen
 *
 */
public class CLI {

	/**
	 * Class for a {@link Protocol} that prints messages to {@link System#out} and warnings to {@link System#err}
	 * and additionally hold a log of all messages in a local {@link StringBuffer} that may be, e.g., stored to a file later. 
	 * @author Jan Grau
	 *
	 */
	public static class SysProtocol implements Protocol {

		private StringBuffer log;
		protected PrintStream out, err;
		/**
		 * Creates a new, empty protocol.
		 */
		public SysProtocol() {
			this.log = new StringBuffer();
			out = System.out;
			err = System.err;
		}
		
		@Override
		public void append( String str ) {
			out.print(str);
			log.append( str );
		}

		@Override
		public void appendHeading( String heading ) {
			out.print("* "+heading);	
			log.append( "* "+heading );
		}

		@Override
		public void appendWarning( String warning ) {
			err.print(warning);
			log.append( warning );
		}
		
		
		/**
		 * Returns the {@link StringBuffer} containing all messages since the creation of this
		 * object.
		 * @return the messages
		 */
		public StringBuffer getLog(){
			return log;
		}

		@Override
		public void appendThrowable(Throwable th) {
			StringWriter str = new StringWriter();
			
			th.printStackTrace(new PrintWriter(str));
			
			String strstr = str.toString();
			err.println(strstr);
			log.append(strstr);
		}
		
		public void appendVerbatim(String verbatim) {
			append(verbatim);
		}

		public void flush() {
			err.flush();
			out.flush();
		}
	}
	
	/**
	 * A dummy class that writes no protocol.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class QuietSysProtocol extends SysProtocol {
		public QuietSysProtocol() {
			out = new PrintStream( SafeOutputStream.getSafeOutputStream(null) );
		}
	}
	
	private String opt;
	private JstacsTool[] tools;
	private boolean[] configureThreads;
	
	private ToolParameterSet[] toolParameters;
	private HashMap<String,String>[] keyMap;
	private String description;
	
	/**
	 * Creates a new command line interface for the tools provided.
	 * @param tools the tools
	 * @throws CloneNotSupportedException 
	 */
	public CLI(JstacsTool... tools) throws CloneNotSupportedException {
		this(null,tools);
	}
	
	/**
	 * Creates a new command line interface from a set of Jstacs tools.
	 * @param configureThreads if the tool at the corresponding index should be configured to use multiple threads, which also means that a parameter <code>threads</code> is displayed for this tool
	 * @param tools (an array of) tool(s) that can be run via the command line interface
	 * @throws CloneNotSupportedException 
	 */
	public CLI(boolean[] configureThreads, JstacsTool... tools) throws CloneNotSupportedException {
		this("",configureThreads,tools);
	}
	
	/**
	 * Creates a new command line interface for the tools provided, where for
	 * each tool multi-threading may be configured. Additional, required arguments when calling the tools
	 * may be specified.
	 * @param opt the additional arguments, may be an empty {@link String}
	 * @param configureThreads for each tool, if multi-threading option should be added to the list of parameters
	 * @param tools the tools
	 * @throws CloneNotSupportedException 
	 */
	public CLI( String opt, boolean[] configureThreads, JstacsTool... tools) throws CloneNotSupportedException {
		this(null, opt, configureThreads, tools);
	}
	
	/**
	 * Creates a new command line interface for the tools provided, where for
	 * each tool multi-threading may be configured. Additional, required arguments when calling the tools
	 * may be specified.
	 * @param description a description of this {@link CLI} displayed on the command line
	 * @param opt the additional arguments, may be an empty {@link String}
	 * @param configureThreads for each tool, if multi-threading option should be added to the list of parameters
	 * @param tools the tools
	 * @throws CloneNotSupportedException 
	 */
	public CLI( String description, String opt, boolean[] configureThreads, JstacsTool... tools) throws CloneNotSupportedException {
		if(configureThreads == null){
			this.configureThreads = new boolean[tools.length];
		}else{
			this.configureThreads = configureThreads;
		}
		this.tools = tools;
		this.toolParameters = new ToolParameterSet[tools.length];
		this.keyMap = new HashMap[tools.length];
		for(int i=0;i<tools.length;i++){
			toolParameters[i] = tools[i].getToolParameters();
			keyMap[i] = new HashMap<String, String>();
			//System.out.println(tools[i].getToolName() + "===================================================");
			addToKeyMap("","",keyMap[i],toolParameters[i]);
		}
		this.opt = opt.trim();
		if( this.opt.indexOf(' ') >= 0 ) {
			throw new IllegalArgumentException(opt);
		} else if( this.opt.length()>0 ) {
			this.opt = " " +this.opt;
		}
		this.description=description;
	}
	
	
	
	
	private static void addToKeyMap( String pathPrefix, String keyPrefix, HashMap<String, String> hashMap, ParameterSet parameterSet ) throws CloneNotSupportedException {
		boolean isExp = parameterSet instanceof ExpandableParameterSet;
		for(int i=0;i<(isExp?1:parameterSet.getNumberOfParameters());i++){
			Parameter par;
			if( isExp ) {
				ExpandableParameterSet e = (ExpandableParameterSet) parameterSet;
				boolean intermediate = parameterSet.getNumberOfParameters()==0;
				if( intermediate ) {
System.out.println("before extension " + e.getNumberOfParameters());
					e.addParameterToSet();
System.out.println("after extension " + e.getNumberOfParameters());
				}
				par = parameterSet.getParameterAt( 0 );
				if( intermediate ) {
					e.removeParameterFromSet();
System.out.println("after deletion " + e.getNumberOfParameters());
				}
				
			} else {
				par = parameterSet.getParameterAt( i );
			}
			
			String parKey = pathPrefix+":"+par.getName();
			String add = isExp ? "?" : i + "";
			if(par.getDatatype() != DataType.PARAMETERSET){
				String key = hashMap.get(parKey); 
				if(key == null){	
					key = getKey(keyPrefix,hashMap,par);
					hashMap.put( parKey, key );
				}
			}else{
				if(par instanceof AbstractSelectionParameter){
					if(par instanceof MultiSelectionParameter){
						throw new RuntimeException("Not implemented for "+par.getClass() );
					}
					
					String key = hashMap.get(parKey);
					if(key == null){
						key = getKey(keyPrefix,hashMap,par);
						hashMap.put( parKey, key );
					}
					ParameterSet incoll = ( (AbstractSelectionParameter)par ).getParametersInCollection();
					HashMap<String,String> copyTemp = (HashMap)hashMap.clone();
					for(int j=0;j<incoll.getNumberOfParameters();j++){
						
						HashMap<String,String> copy = (HashMap) copyTemp.clone();
						
						ParameterSetContainer cont = (ParameterSetContainer)incoll.getParameterAt( j );
						ParameterSet ps = cont.getValue();
						if( ps instanceof ToolParameterSet ) {
							addToKeyMap( pathPrefix+":"+add+"-"+j, keyPrefix + ((ToolParameterSet)ps).getToolName() + ".", copy, ps );
						} else {
							addToKeyMap( pathPrefix+":"+add+"-"+j, keyPrefix, copy, ps );
						}
						
						
						hashMap.putAll(copy);
						
					}
				}else{
					ParameterSet ps = (ParameterSet)par.getValue();
					if( ps instanceof ToolParameterSet ) {
						addToKeyMap( pathPrefix+":"+add, keyPrefix + ((ToolParameterSet)ps).getToolName() + ".", hashMap, ps );
					} else {
						addToKeyMap( pathPrefix+":"+add, keyPrefix, hashMap, ps );
					}
				}
			}
		}
	}



	private static String getKey( String keyPrefix, HashMap<String, String> hashMap, Parameter par ) {
		Collection<String> valueSet2 = hashMap.values();
		LinkedList<String> valueSet = new LinkedList<String>( valueSet2 );
		valueSet.add( "outdir" );
		valueSet.add( "info" );
		valueSet.add( "threads" );
		String parName = par.getName();
		String key = keyPrefix+(parName.charAt( 0 )+"").toLowerCase();
		if(!valueSet.contains( key )){
			return key;
		}
		
		String[] temp = parName.split( "\\s" );
		key = "";
		for(int i=0;i<temp.length;i++){
			key += temp[i].charAt( 0 );
		}
		key = keyPrefix+key.toLowerCase();
		if(!valueSet.contains( key )){
			return key;
		}
		
		key = keyPrefix+parName.replaceAll( "[\\s=]", "" );//TODO .toLowerCase();
		int k=1;
		String temp2 = key;
		while(valueSet.contains( temp2 )){
			temp2 = key+k;
			k++;
		}
		return key;
	}


	private int getToolIndex(String shortName){
		for(int i=0;i<tools.length;i++){
			if(shortName.equals( tools[i].getShortName() )){
				return i;
			}
		}
		return -1;
	}
	
	/**
	 * Runs this command line application with the specified arguments in <code>args</code>. This method would typically be called directly from the <code>main</code>
	 * method.
	 * @param args the arguments as supplied on the command line
	 * @throws Exception in case of an error
	 */
	public void run(String[] args) throws Exception {
		SysProtocol protocol = new SysProtocol();
		
		String outdir = ".";
		int toolIndex = tools.length == 1 ? 0 : ( args.length>0 ? getToolIndex( args[0] ) : -1);
		
		if( toolIndex < 0 ){
			//show summary
			String jar = "<name>.jar";
			
			try{
				File jarfile = new File(CLI.class.getProtectionDomain().getCodeSource().getLocation().toURI().getPath());
				jar = jarfile.getName();
			}catch(Exception e){ }
			
			if( description!=null && description.length()>0) {
				System.err.println( description + "\n" );
			}
			System.err.println( "Available tools:\n" );
			for(int i=0;i<tools.length;i++){
				System.err.println("\t"+tools[i].getShortName()+" - "+tools[i].getToolName());
			}
			System.err.println();
			System.err.println("Syntax: java -jar "+jar+opt+" <toolname> [<parameter=value> ...]\n");
			System.err.println("Further info about the tools is given with\n\tjava -jar "+jar+opt+" <toolname> info\n");
			System.err.println("For tests of individual tools:\n\tjava -jar "+jar+opt+" <toolname> test [<verbose>]\n");
			System.err.println("Tool parameters are listed with\n\tjava -jar "+jar+opt+" <toolname>\n");
		} else if( (args.length > 1 && args[1].equals( "info" ) ) || (tools.length == 1 && args.length>0 && args[0].equals("info"))) {
			//show info
			System.err.println("\n"+parse(tools[toolIndex].getHelpText()));
			String[] ref = tools[toolIndex].getReferences();
			int n = ref==null?0:ref.length;
			if( n>0 ) {
				System.err.println("If you use this tool, please cite");
				for( int i = 0; i<n; i++  ) {
					System.err.print("\n"+parse(ref[i]));
				}
			}
		} else if( (args.length > 1 && args[1].equals( "test" ) ) || (tools.length == 1 && args.length>0 && args[0].equals("test"))) {
			//run tests
			System.err.println("\nRun tests for "+tools[toolIndex].getShortName() + ":");
			boolean verbose = false;
			int i = (args.length > 1 && args[1].equals( "test" )) ? 2 : 3;
			if( args.length >= i+1 ) {
				verbose = Boolean.parseBoolean(args[i]);
			}
			double val = JstacsTool.test(tools[toolIndex],verbose);
			System.err.println("\nSummary: " + val + "\n" );
		} else {
			//run tool
			Pair<String,Integer> pair = null;
			int threads = 1;
			try{
				 pair = setToolParameters(configureThreads[toolIndex], tools.length == 1 ? 0 : 1, toolParameters[toolIndex],keyMap[toolIndex],args, protocol);
				 outdir = pair.getFirstElement();
				 threads = pair.getSecondElement();
			}catch(IllegalArgumentException e){
				printToolParameters(toolIndex,protocol,outdir,1);
				
				System.err.println("\n"+e.getMessage());
				return;
			}
			
			
			
			
			if(!toolParameters[toolIndex].hasDefaultOrIsSet()){
				printToolParameters(toolIndex,protocol,outdir,threads);
				System.err.println("\nAt least one parameter has not been set (correctly).");
				String msg = toolParameters[toolIndex].getErrorMessage();
				if( msg != null ) System.err.println(msg);
				
				return;
			}else{
				printToolParameters(toolIndex,protocol,outdir,threads);
			}
			protocol.flush();
			
			ToolResult results = tools[toolIndex].run( toolParameters[toolIndex], protocol, new ProgressUpdater(), threads );
			
			writeToolResults( results, protocol, outdir, tools[toolIndex], toolParameters[toolIndex] );
		}		
	}
	
	public static void writeToolResults( ToolResult results, SysProtocol protocol, String outdir, JstacsTool tool, ToolParameterSet toolParameters ) throws IOException {
		Result[] r;
		if( results != null ) {
			ResultSet[] rs = results.getRawResult();
			r = new Result[rs[0].getNumberOfResults()+1];
			for( int i = 0; i+1 < r.length; i ++ ) {
				r[i] = rs[0].getResultAt(i);
			}
		} else {
			r = new Result[1];
		}
		r[r.length-1] = new TextResult( "protocol " + tool.getShortName(), "Result", new FileRepresentation("CLI-protocol.txt",protocol.getLog().toString()), "txt", tool.getToolName(), null, true );
		
		results = new ToolResult(
				results==null?"":results.getName(),
				results==null?"":results.getComment(),
				results==null?null:results.getAnnotation(), 
				new ResultSet( r ),
				toolParameters,
				results==null?"":results.getToolName(),
				results==null?new Date():results.getFinishedDate() );
		
		ResultSaver saver = ResultSaverLibrary.getSaver( results.getClass() );	
		saver.writeOutput( results, new File(outdir) );

/*
		String prefix = outdir+File.separator+"protocol_" + tools[toolIndex].getShortName();//new
		File protout = new File( prefix + ".txt");
		int k=1;
		while(protout.exists()){
			protout = new File(prefix + "_"+k+".txt");
			k++;
		}
		
		FileManager.writeFile( protout, protocol.getLog() );
*/
	}


	private Pair<String,Integer> setToolParameters( boolean configureThreads, int off, ParameterSet parameterSet, HashMap<String, String> hashMap, String[] args, Protocol protocol ) throws IllegalValueException, CloneNotSupportedException {
		HashMap<String, LinkedList<String>> valueMap = new HashMap<String, LinkedList<String>>();
		String outdir = ".";
		int threads = 1;
		boolean newLine=false;
		for(int i=off;i<args.length;i++){
			int idx = args[i].indexOf("=");
			if(idx < 0){
				throw new IllegalValueException("Parameter mis-specified in: "+args[i]);
			}
			String[] temp = new String[]{args[i].substring(0, idx),args[i].substring(idx+1)};

			if("outdir".equals( temp[0] ) ){
				outdir = temp[1];
			}else if(configureThreads && "threads".equals( temp[0] ) ){
				threads = Integer.parseInt(temp[1]);
			}else{
				LinkedList<String> v = valueMap.get(temp[0]);
				if(v == null){
					valueMap.put(temp[0], new LinkedList<String>());
				}
				/*if( v!=null ) {
					protocol.appendWarning( "Overwriting parameter: " + temp[0]+"=" +v +"\n");
					newLine=true;
				}*/
				valueMap.get(temp[0]).add(temp[1]);
			}
		}

		if( newLine ) {
			protocol.append("\n");
		}
		set("",parameterSet,hashMap,valueMap,protocol,0);
		if( valueMap.size() > 0 ) {
			Iterator<String> it = valueMap.keySet().iterator();
			HashSet<String> h = new HashSet<String>();
			h.addAll( hashMap.values() );
			while( it.hasNext() ) {
				String k = it.next();
				if( h.contains(k) ) {
					
					throw new IllegalArgumentException( "Parameter "+k+" specified multiple times or not applicable for current selection(s). Could not use value(s): "+ valueMap.get(k) );
				}
			}
			throw new IllegalArgumentException("Unknown parameters: "+ valueMap );
		}
		return new Pair<String, Integer>(outdir, threads);
	}

	private void set( String pathPrefix, ParameterSet parameters, HashMap<String, String> hashMap, HashMap<String, LinkedList<String>> valueMap, Protocol protocol, int exp ) throws IllegalValueException, CloneNotSupportedException {
		boolean isExp = parameters instanceof ExpandableParameterSet;
		ParameterSet template = null;
		boolean intermediate = false;
		ExpandableParameterSet eps = null;
		if( isExp ) {
			eps = (ExpandableParameterSet) parameters;
			intermediate=eps.getNumberOfParameters() == 0;
			if( intermediate ) {
				eps.addParameterToSet();
			}
			template=(ParameterSet) parameters.getParameterAt(0).getValue();
			exp++;
			if( exp > 1 ) {
				throw new RuntimeException("Nested ExpandableParameterSets not implemented.");//TODO
			}
			pathPrefix += ":?";
		}
		int n = (isExp ? template : parameters).getNumberOfParameters();
		for(int i=0;i<n;i++){
			Parameter par = (isExp ? template : parameters).getParameterAt( i );
			String parKey = pathPrefix+":"+par.getName();
			//System.out.println(parKey+" -> "+par.getName());
			
			if(par.getDatatype() != DataType.PARAMETERSET || par instanceof AbstractSelectionParameter){
				String key = hashMap.get( parKey );
if( key == null ) {
	throw new IllegalArgumentException("Could not find a key for parameter: " + par.getName() + "\t" + parKey);
}
				LinkedList<String> value = valueMap.get( key );
				
				if(value != null){
					if( isExp ) {
						for(int k=0;value.size()>0;k++){
							if(k >= parameters.getNumberOfParameters()){
								try {
									eps.addParameterToSet();
								} catch (CloneNotSupportedException doesnothappen) {
									doesnothappen.printStackTrace();
								}
							}
							Parameter par2 = ((ParameterSet) parameters.getParameterAt(k).getValue()).getParameterAt(i);
							par2.setValue(value.removeFirst());
							
							if( par2.getDatatype() == DataType.PARAMETERSET ) {
								set(pathPrefix+":"+i+(par2 instanceof SelectionParameter ? "-"+((SelectionParameter)par2).getSelected() : ""),(ParameterSet)par2.getValue(),hashMap,valueMap,protocol,exp);
							}
						}						
					} else {
						par.setValue( value.removeFirst() );
						if( par.getDatatype() == DataType.PARAMETERSET ) {
							set(pathPrefix+":"+i+(par instanceof SelectionParameter ? "-"+((SelectionParameter)par).getSelected() : ""),(ParameterSet)par.getValue(),hashMap,valueMap,protocol, exp );
						}
					}
				} else if( par.getDatatype() == DataType.PARAMETERSET ) {
					set(pathPrefix+":"+i+(par instanceof SelectionParameter ? "-"+((SelectionParameter)par).getSelected() : ""),(ParameterSet)par.getValue(),hashMap,valueMap,protocol, exp );
				}
				
				if( value != null && value.size() == 0 ) {
					valueMap.remove(key);
				}
			} else {
				set(pathPrefix+":"+i,(ParameterSet) par.getValue(), hashMap, valueMap, protocol, exp );
			}
		}
		if( intermediate ) {
			eps.removeParameterFromSet();
		}
	}




	public static void print(String pathPrefix, HashMap<String, String> keyMap, ParameterSet parameters, String tabPrefix, Protocol protocol, String add ) throws CloneNotSupportedException{
		boolean isExp = parameters instanceof ExpandableParameterSet;
		ExpandableParameterSet exp = null;
		boolean intermediate = false;
		if( isExp ) {
			exp = (ExpandableParameterSet) parameters;
			intermediate = exp.getNumberOfParameters()==0;
			if( intermediate ) exp.addParameterToSet();
			parameters=(ParameterSet) parameters.getParameterAt(0).getValue();
			
			protocol.appendWarning( tabPrefix+"The following parameter(s) can be used " + (intermediate?"zero or ":"") + "multiple times:\n" );//TODO
			int n = exp.getNumberOfParameters();
			for(int k=0;k<n;k++){
				ParameterSet ps2 = (ParameterSet) exp.getParameterAt(k).getValue();
				print(pathPrefix+":?", keyMap, ps2, tabPrefix, protocol, n>1?" ("+(k+1)+")":"" );
			}
		} else {		
			for(int i=0;i<parameters.getNumberOfParameters();i++){
				Parameter par = parameters.getParameterAt( i );
				String parKey = pathPrefix+":"+par.getName();
				String k = keyMap==null ? "" : keyMap.get( parKey );
				if(par.getDatatype() != DataType.PARAMETERSET){
					protocol.appendWarning( tabPrefix+k+add+" - "+par.toString()+"\n" );
				}else{
					if(par instanceof AbstractSelectionParameter){
if( k == null ) {
	throw new IllegalArgumentException("Could not find a key for parameter: " + par.getName() + "\t" + parKey);
}
						protocol.appendWarning( tabPrefix+k+add+" - "+par.toString()+"\n" );
						ParameterSet incoll = ( (AbstractSelectionParameter)par ).getParametersInCollection();
						char[] array = new char[k.length()+3];
						Arrays.fill(array,' ');
						String off = tabPrefix+(new String(array));
						for(int j=0;j<incoll.getNumberOfParameters();j++){
							ParameterSetContainer cont = (ParameterSetContainer)incoll.getParameterAt( j );
							if( cont.getValue().getNumberOfParameters()>0 ) {
								protocol.appendWarning( off+"Parameters for selection \""+cont.getName()+"\":\n" );
								print(pathPrefix+":"+i+"-"+j, keyMap,cont.getValue(),off+"\t",protocol,add);
							} else {
								protocol.appendWarning( off+"No parameters for selection \""+cont.getName()+"\"\n" );
							}
						}
					} else {
						ParameterSet ps = (ParameterSet)par.getValue();
						boolean b = ps instanceof ToolParameterSet;
						if( b ) protocol.appendWarning("\n" + ((ToolParameterSet)ps).getToolName() + " parameters" + "\n" );
						print(pathPrefix+":"+i,keyMap,ps,tabPrefix+(!b?"\t":""),protocol,add);
						if( b ) protocol.appendWarning("\n");
					}
				}
			}
		}
		if( intermediate ) {
			exp.removeParameterFromSet();
		}
	}
	
	private static void printTable(String prefix,HashMap<String, String> keyMap, ParameterSet parameters, StringBuffer out) throws CloneNotSupportedException{
		boolean isExp = parameters instanceof ExpandableParameterSet;
		boolean intermediate = false;
		ExpandableParameterSet e = null;
		if( isExp ) {
			e = (ExpandableParameterSet) parameters;
			intermediate = e.getNumberOfParameters()==0;
			if( intermediate ) e.addParameterToSet();
			out.append( "<tr><td colspan=3>The following parameter(s) can be used " + (intermediate?"zero or ":"") + "multiple times:</td></tr>\n"
					+"<tr><td></td><td colspan=2><table border=0 cellpadding=0 align=\"center\" width=\"100%\">\n");
		}
		for(int i=0;i<parameters.getNumberOfParameters();i++){
			Parameter par = parameters.getParameterAt( i );
			String parKey = prefix+":"+par.getName();
			String add = isExp ? "?" : i + "";
			if(par.getDatatype() == DataType.PARAMETERSET && ! (par instanceof AbstractSelectionParameter) ){
				ParameterSet ps = (ParameterSet)par.getValue();
				printTable(prefix+":"+add, keyMap, ps, out);
			} else {
				out.append( "<tr style=\"vertical-align:top\">\n<td><font color=\"green\">" + keyMap.get( parKey )+ "</font></td>\n" );
				String s = par.toString();
				out.append( "<td>"+s.substring(0,s.lastIndexOf(")\t= ")+1) + "</td>\n" );
				if( par.getDatatype() != DataType.PARAMETERSET ) {
					out.append( "<td style=\"width:100px;\">"+par.getDatatype() + "</td>\n</tr>\n" );
				} else {
					out.append( "<td style=\"width:100px;\"></td></tr>");
					out.append( "<tr><td></td><td colspan=2><table border=0 cellpadding=0 align=\"center\" width=\"100%\">\n" );
					ParameterSet incoll = ( (AbstractSelectionParameter)par ).getParametersInCollection();
					for(int j=0;j<incoll.getNumberOfParameters();j++){
						ParameterSetContainer cont = (ParameterSetContainer)incoll.getParameterAt( j );
						if( cont.getValue().getNumberOfParameters()>0 ) {
							out.append( "<tr><td colspan=3><b>Parameters for selection &quot;"+cont.getName()+"&quot;:</b></td></tr>\n" );
							printTable(prefix+":"+add+"-"+j, keyMap, cont.getValue(), out);
						} else {
							out.append( "<tr><td colspan=3><b>No parameters for selection &quot;"+cont.getName()+"&quot;</b></td></tr>\n" );
						}
					}
					out.append( "</table></td></tr>\n" );
				}
			}
		}
		if( isExp ) {
			out.append( "</table>\n</td></tr>\n");
			if( intermediate ) e.removeParameterFromSet();
		}
	}

	private void printToolParameters( int toolIndex, Protocol protocol, String outdir, int threads ) throws CloneNotSupportedException {
		ParameterSet ps = toolParameters[toolIndex];
		protocol.appendWarning( "Parameters of tool \""+tools[toolIndex].getToolName()+"\" ("+tools[toolIndex].getShortName()+", version: " + tools[toolIndex].getToolVersion() + "):\n" );
		print( "",keyMap[toolIndex], ps, "", protocol, "" );
		protocol.appendWarning( "outdir - The output directory, defaults to the current working directory (.)\t= "+outdir+"\n" );
		if(configureThreads[toolIndex]){
			protocol.appendWarning( "threads - The number of threads used for the tool, defaults to 1\t= "+threads+"\n" );
		}
	}
	
	
	/**
	 * Creates the HTML code for wiki entries for the supplied tools.
	 */
	public void wiki() throws IOException {
		wiki("");
	}
	

	/**
	 * Creates the HTML code for wiki entries for the supplied tools.
	 * 
	 * @param home the home directory of the files
	 */
	public void wiki( String home ) {
		File h = new File( home+"/" );
		h.mkdirs();
		for( int toolIndex=0; toolIndex<keyMap.length; toolIndex++ ) {
			ParameterSet ps = toolParameters[toolIndex];
			try {
				StringBuffer old = null;
				File f = new File(home+"/"+tools[toolIndex].getShortName()+".txt");
				if( f.exists() ) {
					old = FileManager.readFile(f);
				}
				StringBuffer current = new StringBuffer();
				current.append( "<table border=0 cellpadding=10 align=\"center\" width=\"100%\">\n<tr>\n<td>name</td>\n<td>comment</td>\n<td>type</td>\n</tr>\n<tr><td colspan=3><hr></td></tr>\n" );
				printTable("", keyMap[toolIndex], ps, current);
				current.append( "<tr style=\"vertical-align:top\">\n<td><font color=\"green\">outdir</font></td>\n" );
				current.append( "<td>The output directory, defaults to the current working directory (.)</td>\n" );
				current.append( "<td>STRING</td>\n</tr>\n" );
				if(configureThreads[toolIndex]){
					current.append( "<tr style=\"vertical-align:top\">\n<td><font color=\"green\">threads</font></td>\n" );
					current.append( "<td>The number of threads used for the tool, defaults to 1</td>\n" );
					current.append( "<td>INT</td>\n</tr>\n" );
				}				
				current.append( "</table>" );
				if( old == null || !current.toString().equals(old.toString()) ) {
					FileManager.writeFile(f, current);
					if( old != null ) FileManager.writeFile(home+"/"+tools[toolIndex].getShortName()+".bak", old);
				}
			} catch( Exception ex ) {
				//nothing
			}
		}
	}
	
	
	
	private static String parse(String restruct){
		String[] lines = restruct.split( "\n" );
		
		Pattern bold = Pattern.compile( "\\*\\*(.+?)\\*\\*" );
		Pattern italics = Pattern.compile( "\\*(.+?)\\*" );
		Pattern tt = Pattern.compile( "\\`\\`(.+?)\\`\\`" );
		
		Pattern link = Pattern.compile( "^\\.\\.\\s+\\_(.*?)\\s*\\:\\s*(.*)$" );
		
		HashMap<Pattern, String> linkTargets = new HashMap<Pattern, String>();
		
		for(int i=0;i<lines.length;i++){
			
			Matcher m = bold.matcher( lines[i] );
			lines[i] = m.replaceAll( "§$1§" );
			
			m = italics.matcher( lines[i] );
			lines[i] = m.replaceAll( "'$1'" );
			
			lines[i] = lines[i].replaceAll( "§", "*" );
			
			m = tt.matcher( lines[i] );
			lines[i] = m.replaceAll( "$1" );	
			
			m = link.matcher( lines[i] );
			
			if(m.matches()){
				String key = m.group( 1 );
				String target = m.group( 2 );
				linkTargets.put( Pattern.compile( "\\`?("+key+")\\`?\\_" ), target );
				lines[i] = "";
			}
			
		}
		
		
		
		for(int i=0;i<lines.length;i++){
			Set<Pattern> pats = linkTargets.keySet();
			Iterator<Pattern> it = pats.iterator();
			while( it.hasNext() ){
				Pattern pat = it.next();
				Matcher m = pat.matcher( lines[i] );
				lines[i] = m.replaceAll( "$1 ("+linkTargets.get( pat )+")" );
			}
		}
		
		StringBuffer sb = new StringBuffer();
		for(int i=0;i<lines.length;i++){
			sb.append( lines[i] );
			sb.append( "\n" );
		}
		return sb.toString();
	}

}
