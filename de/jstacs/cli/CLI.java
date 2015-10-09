package de.jstacs.cli;

import java.io.File;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.DataType;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.AbstractSelectionParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;
import de.jstacs.results.savers.ResultSaver;
import de.jstacs.results.savers.ResultSaverLibrary;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;


public class CLI {

	public static class SysProtocol implements Protocol{

		private StringBuffer log;
		
		public SysProtocol() {
			this.log = new StringBuffer();
		}
		
		@Override
		public void append( String str ) {
			System.out.print(str);
			log.append( str );
		}

		@Override
		public void appendHeading( String heading ) {
			System.out.print("* "+heading);	
			log.append( "* "+heading );
		}

		@Override
		public void appendWarning( String warning ) {
			System.err.print(warning);
			log.append( warning );
		}
		
		
		
		public StringBuffer getLog(){
			return log;
		}

		@Override
		public void appendThrowable(Throwable th) {
			StringWriter str = new StringWriter();
			
			th.printStackTrace(new PrintWriter(str));
			
			String strstr = str.toString();
			System.err.println(strstr);
			log.append(strstr);
		}
		
	}
	
	
	private JstacsTool[] tools;
	
	private ParameterSet[] toolParameters;
	private HashMap<Parameter,String>[] keyMap;
	
	public CLI(JstacsTool... tools) {
		this.tools = tools;
		this.toolParameters = new ParameterSet[tools.length];
		this.keyMap = new HashMap[tools.length];
		for(int i=0;i<tools.length;i++){
			toolParameters[i] = tools[i].getToolParameters();
			keyMap[i] = new HashMap<Parameter, String>();
			addToKeyMap(keyMap[i],toolParameters[i]);
		}
		
	}
	
	
	
	
	private void addToKeyMap( HashMap<Parameter, String> hashMap, ParameterSet parameterSet ) {
		for(int i=0;i<parameterSet.getNumberOfParameters();i++){
			Parameter par = parameterSet.getParameterAt( i );
			if(par.getDatatype() != DataType.PARAMETERSET){
				String key = getKey(hashMap,par);
				hashMap.put( par, key );
			}else{
				if(par instanceof AbstractSelectionParameter){
					String key = getKey(hashMap,par);
					hashMap.put( par, key );
					ParameterSet incoll = ( (AbstractSelectionParameter)par ).getParametersInCollection();
					for(int j=0;j<incoll.getNumberOfParameters();j++){
						ParameterSetContainer cont = (ParameterSetContainer)incoll.getParameterAt( j );
						addToKeyMap( hashMap, cont.getValue() );
					}
				}else{
					ParameterSet ps = (ParameterSet)par.getValue();
					addToKeyMap( hashMap, ps );
				}
			}
		}		
	}



	private String getKey( HashMap<Parameter, String> hashMap, Parameter par ) {
		Collection<String> valueSet2 = hashMap.values();
		LinkedList<String> valueSet = new LinkedList<String>( valueSet2 );
		valueSet.add( "outdir" );
		valueSet.add( "info" );
		String parName = par.getName();
		String key = (parName.charAt( 0 )+"").toLowerCase();
		if(!valueSet.contains( key )){
			return key;
		}
		
		String[] temp = parName.split( "\\s" );
		key = "";
		for(int i=0;i<temp.length;i++){
			key += temp[i].charAt( 0 );
		}
		key = key.toLowerCase();
		if(!valueSet.contains( key )){
			return key;
		}
		
		key = parName.replaceAll( "[\\s=]", "" );
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

	
	
	
	public void run(String[] args) throws Exception {
		SysProtocol protocol = new SysProtocol();
		
		String jar = "<name>.jar";
		
		try{
			File jarfile = new File(CLI.class.getProtectionDomain().getCodeSource().getLocation().toURI().getPath());
			jar = jarfile.getName();
		}catch(Exception e){ }
		
		String outdir = ".";
		
		if(args.length == 0 || (tools.length > 1 && getToolIndex( args[0] ) < 0)){
			System.err.println( "Available tools:\n" );
			for(int i=0;i<tools.length;i++){
				System.err.println("\t"+tools[i].getShortName()+" - "+tools[i].getToolName());
			}
			System.err.println();
			System.err.println("Syntax: java -jar "+jar+" <toolname> [<parameter=value> ...]\n");
			System.err.println("Further info about the tools is given with\n\tjava -jar "+jar+" <toolname> info\n");
			System.err.println("Tool parameters are listed with\n\tjava -jar "+jar+" <toolname>\n");
			return;
		}/*else if( ( tools.length == 1 && args.length==0) || args.length == 1){
			int toolIndex = getToolIndex( args[0] );
			printToolParameters(toolIndex,protocol,".");
			System.err.println( "Syntax: java -jar "+jar+" [<parameter=value> ...]" );
			return;
		}*/else if(args.length > 1 && args[1].equals( "info" )){
			int toolIndex = getToolIndex( args[0] );
			System.err.println("\n"+parse(tools[toolIndex].getHelpText()));
		}else{
			int toolIndex = getToolIndex( args[0] );
			outdir = setToolParameters(tools.length == 1 ? 0 : 1, toolParameters[toolIndex],keyMap[toolIndex],args);
			
			if(!toolParameters[toolIndex].hasDefaultOrIsSet()){
				System.err.println("At least one parameter has not been set (correctly):\n");
				printToolParameters(toolIndex,protocol,outdir);
				return;
			}else{
				printToolParameters(toolIndex,protocol,outdir);
			}
			
			
			ToolResult results = tools[toolIndex].run( toolParameters[toolIndex], protocol, new ProgressUpdater() );
			
		//	ResultSetResult res = new ResultSetResult( "Tool results", "", null, results );
			
			ResultSaver saver = ResultSaverLibrary.getSaver( results );
			
			saver.writeOutput( results, new File(outdir) );
			
			File protout = new File(outdir+File.separator+"protocol.txt");
			int k=1;
			while(protout.exists()){
				protout = new File(outdir+File.separator+"protocol_"+k+".txt");
				k++;
			}
			
			FileManager.writeFile( protout, protocol.getLog() );
			
		}
		
		
	}


	private String setToolParameters( int off, ParameterSet parameterSet, HashMap<Parameter, String> hashMap, String[] args ) throws IllegalValueException {
		HashMap<String, String> valueMap = new HashMap<String, String>();
		String outdir = ".";
		for(int i=off;i<args.length;i++){
			int idx = args[i].indexOf("=");
			if(idx < 0){
				throw new IllegalValueException("Parameter mis-specified in: "+args[i]);
			}
			String[] temp = new String[]{args[i].substring(0, idx),args[i].substring(idx+1)};
			
			if("outdir".equals( temp[0] ) ){
				outdir = temp[1];
			}else{
				valueMap.put( temp[0], temp[1] );
			}
		}
		
		set(parameterSet,hashMap,valueMap);
		
		return outdir;
	}




	private void set( ParameterSet parameters, HashMap<Parameter, String> hashMap, HashMap<String, String> valueMap ) throws IllegalValueException {
		for(int i=0;i<parameters.getNumberOfParameters();i++){
			Parameter par = parameters.getParameterAt( i );
			if(par.getDatatype() != DataType.PARAMETERSET){
				String key = hashMap.get( par );
				String value = valueMap.get( key );
				if(value != null){
					par.setValue( value );
				}
			}else{
				if(par instanceof AbstractSelectionParameter){
					String key = hashMap.get( par );
					String value = valueMap.get( key );
					if(value != null){
						par.setValue( value );
						ParameterSet set = (ParameterSet)par.getValue();
						set(set,hashMap,valueMap);
					}
				}else{
					ParameterSet ps = (ParameterSet)par.getValue();
					set(ps,hashMap,valueMap);
				}
			}
		}
	}




	private void print(HashMap<Parameter, String> keyMap, ParameterSet parameters, String tabPrefix, Protocol protocol){
		for(int i=0;i<parameters.getNumberOfParameters();i++){
			Parameter par = parameters.getParameterAt( i );
			if(par.getDatatype() != DataType.PARAMETERSET){
				protocol.appendWarning( tabPrefix+keyMap.get( par )+" - "+par.toString()+"\n" );
			}else{
				if(par instanceof AbstractSelectionParameter){
					protocol.appendWarning( tabPrefix+keyMap.get( par )+" - "+par.toString()+"\n" );
					ParameterSet incoll = ( (AbstractSelectionParameter)par ).getParametersInCollection();
					for(int j=0;j<incoll.getNumberOfParameters();j++){
						ParameterSetContainer cont = (ParameterSetContainer)incoll.getParameterAt( j );
						protocol.appendWarning( "\n"+tabPrefix+"\tParameters for selection \""+cont.getName()+"\":\n" );
						print(keyMap,cont.getValue(),tabPrefix+"\t",protocol);
					}
				}else{
					ParameterSet ps = (ParameterSet)par.getValue();
					print(keyMap,ps,tabPrefix+"\t",protocol);
				}
			}
		}
	}

	private void printToolParameters( int toolIndex, Protocol protocol, String outdir ) {
		ParameterSet ps = toolParameters[toolIndex];
		protocol.appendWarning( "Parameters of tool \""+tools[toolIndex].getToolName()+"\" ("+tools[toolIndex].getShortName()+"):\n" );
		print( keyMap[toolIndex], ps, "", protocol );
		protocol.appendWarning( "outdir - The output directory, defaults to the current working directory (.)\t= "+outdir+"\n" );
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
