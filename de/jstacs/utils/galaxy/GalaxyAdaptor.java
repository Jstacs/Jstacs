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

package de.jstacs.utils.galaxy;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.LinkedList;

import javax.imageio.ImageIO;

import de.jstacs.DataType;
import de.jstacs.classifiers.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.io.FileManager;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.DataSetResult;
import de.jstacs.results.ImageResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.MeanResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.SimpleResult;
import de.jstacs.results.StorableResult;

/**
 * Adaptor class between the parameter representation of Jstacs in {@link de.jstacs.parameters.Parameter}s and {@link ParameterSet}s and the parameter representation
 * in <a href="http://galaxy.psu.edu/">Galaxy</a>.
 * A {@link GalaxyAdaptor} can be created from a {@link ParameterSet} containing all parameters that are necessary for
 * the execution of some program that shall be included in a Galaxy installation.
 * 
 * The conversion is done by the method {@link GalaxyAdaptor#parse(String[])}. If this method is called with first argument equal to &quot;--create&quot;
 * and second argument a filename, then the Galaxy-representation of the current {@link ParameterSet} is stored to that file. Afterwards, the application
 * and this file must be added to <code>tool_conf.xml</code> (see <a href="https://bitbucket.org/galaxy/galaxy-central/wiki/AddToolTutorial">Galaxy tutorial</a> for details).
 * 
 * During the execution of the main program, {@link GalaxyAdaptor#addResult(Result, boolean, boolean)} and {@link GalaxyAdaptor#addResultSet(ResultSet, boolean, boolean)}
 * can be used to add results of the computation to the results displayed in Galaxy, to export these and to add these to a summary page. If a
 * protocol of the program run shall be written, a {@link Protocol} object can be obtained by the {@link GalaxyAdaptor#getProtocol(boolean)} method.
 * 
 * After all results have been added and the protocol has been written, the method {@link GalaxyAdaptor#writeOutput()} must be called to create
 * the appropriate output files within Galaxy.
 * 
 * @author Jan Grau
 *
 */
public class GalaxyAdaptor {

	private ParameterSet parameters;
	private boolean[] addLine;
	private String toolname;
	private String description;
	private String help;
	private String command;
	private String version;
	
	private NumberFormat format;
	private NumberFormat expFormat;
	
	private String labelName;
	
	private Protocol protocol;
	private boolean exportProtocol;
	
	private LinkedList<OutputElement> list = new LinkedList<OutputElement>();
	private String outfile;
	private String outfileId;
	private String newFilePath;
	private String htmlFilesPath;
	/**
	 * The stylesheet used for the Galaxy HTML output.
	 */
	public static String stylesheet = "<style type=\"text/css\">" +
			"body{font-family:sans-serif;font-size:10pt}\n" +
			"table{font-size:10pt;border-spacing:0px}\n" +
			"th{font-weight:bold}\n" +
			"div.head{font-size:15px;line-height:24px;padding:5px 10px;background:#ebd9b2;border-bottom:solid #d8b365 1px;font-weight:bold}\n" +
			"div.comment{color:grey}\n" +
			"h2{text-align:center}" +
			"</style>";
	private static int htmlId = 0;
	
	private static String[] colors = new String[]{"#99FFFF","#CCCCFF", "#99FFCC", "#CCFF99", "#FFCC99"};
	
	/**
	 * Returns the color for a specified depth within the parameter hierarchy.
	 * @param depth the depth
	 * @return the color
	 */
	public static String getColor(int depth){
		while(depth >= colors.length){depth -= colors.length;}
		return colors[depth];
	}
	
	/**
	 * Gets a unique id that can be used, e.g. to create unique filenames.
	 * @return the id
	 */
	public static int getHtmlId(){
		return htmlId++;
	}
	
	/**
	 * Creates a new {@link GalaxyAdaptor} from a given {@link ParameterSet} containing all parameters
	 * that are necessary for a program is shall be included in a Galaxy installation. Besides the
	 * parameters, a name, a description, and a version number of the program must be provided.
	 * Additionally, the user must provide the command to run this program. For instance, if the program is bundled into
	 * a jar <code>MyJar.jar</code>, this command could be <code>java -jar MyJar.jar</code>.
	 * @param parameters the parameters of the program
	 * @param addLine indicates for each parameter in <code>parameters</code> if a line is displayed before its name, ignored if <code>parameters</code> is not a {@link SimpleParameterSet}. May be <code>null</code>.
	 * @param toolname the name of the program
	 * @param description a description of the program, may be supplemented by additional help provided by {@link GalaxyAdaptor#setHelp(File)} or {@link GalaxyAdaptor#setHelp(String)}
	 * @param version the version of the program
	 * @param command the command to run the program
	 * @param labelName if <code>null</code>, the default names of Galaxy are used, otherwise a field &quot;Job name&quot; 
	 *     is added as first parameter of the tool with internal name <code>labelName</code>. The internal name must not collide with the name of any other parameter.
	 */
	public GalaxyAdaptor(ParameterSet parameters, boolean[] addLine, String toolname, String description, String version, String command, String labelName){
		this.parameters = parameters;
		this.addLine = addLine == null ? null : addLine.clone();
		this.toolname = toolname;
		this.description = description;
		this.version = version;
		this.help = description;
		this.command = command;
		this.labelName = labelName;
		this.format = NumberFormat.getNumberInstance();
		this.format.setMaximumFractionDigits( 3 );
		this.format.setMinimumFractionDigits( 3 );
		this.expFormat = new DecimalFormat("0.00E0");
	}

	/**
	 * Sets the help, i.e., a more detailed description of the program
	 * to <code>help</code>.
	 * @param help the help text
	 */
	public void setHelp(String help){
		this.help = help;
	}
	
	/**
	 * Sets the help, i.e., a more detailed description of the program
	 * to the contents of <code>helpfile</code>.
	 * @param helpfile the file containing the help text
	 * @throws IOException if <code>helpfile</code> could not be read
	 */
	public void setHelp(File helpfile) throws IOException {
		BufferedReader read = new BufferedReader( new FileReader( helpfile ) );
		String temp = null;
		StringBuffer helpb = new StringBuffer();
		while( (temp=read.readLine()) != null ){
			helpb.append( temp );
			helpb.append( "\n" );
		}
		read.close();
		this.help = helpb.toString();
	}
	
	/**
	 * Returns an object for writing a protocol of a program run
	 * @param exportProtocol if <code>true</code> the protocol will be exported and accessible 
	 * 					as extra result within Galaxy
	 * @return the protocol object
	 */
	public Protocol getProtocol(boolean exportProtocol){
		this.exportProtocol = exportProtocol;
		protocol = new Protocol();
		return protocol;
	}
	
	/**
	 * Creates the contents of a Galaxy configuration file from all the information
	 * provided to this {@link GalaxyAdaptor}.
	 * @return the configuration
	 * @throws Exception if any of the parameters could not be converted
	 */
	public String toGalaxyConfig() throws Exception{
		
		StringBuffer allBuffer = new StringBuffer();
		XMLParser.appendObjectWithTagsAndAttributes( allBuffer, description, "description", null, false );
		allBuffer.append( "\n" );
		XMLParser.appendObjectWithTagsAndAttributes( allBuffer, command+" --run $script_file $summary $summary.id $__new_file_path__ $summary.extra_files_path", "command", null, false );
		allBuffer.append( "\n" );
		
		StringBuffer descBuffer = new StringBuffer();
		StringBuffer confBuffer = new StringBuffer();
		
		if(labelName != null){
			XMLParser.addTagsAndAttributes( descBuffer, "param", "type=\"text\" size=\"40\" name=\""+getLegalName( toolname )+"_"+labelName+"\" label=\"Job name\" value=\"\" optional=\"true\" help=\"Please enter a name for your job that should be used in the history (optional)\"" );
		}
		
		if(parameters instanceof SimpleParameterSet){
			((SimpleParameterSet)parameters).toGalaxy( getLegalName( toolname ) , "", 0, descBuffer, confBuffer, addLine );
		}else{
			parameters.toGalaxy( getLegalName( toolname ) , "", 0, descBuffer, confBuffer, false );
		}
		XMLParser.addTags( descBuffer, "inputs" );
		
		confBuffer = escape( confBuffer );
		XMLParser.addTagsAndAttributes( confBuffer, "configfile", "name=\"script_file\"");
		XMLParser.addTags( confBuffer, "configfiles" );
		allBuffer.append( descBuffer );
		allBuffer.append( confBuffer );
		
		StringBuffer outBuf = new StringBuffer();
		if(labelName != null){
			XMLParser.addTagsAndAttributes( outBuf, "data", "format=\"html\" name=\"summary\" label=\"#if str($"+getLegalName( toolname )+"_"+labelName+") == '' then $tool.name + ' on ' + $on_string else $"+getLegalName( toolname )+"_"+labelName+"#\"" );
		}else{
			XMLParser.addTagsAndAttributes( outBuf, "data", "format=\"html\" name=\"summary\"" );
		}
		XMLParser.addTags( outBuf, "outputs" );
		
		allBuffer.append( outBuf );
		
		StringBuffer helpBuf = new StringBuffer();
		helpBuf.append( help );
		XMLParser.addTags( helpBuf, "help" );
		
		allBuffer.append( helpBuf );
		
		XMLParser.addTagsAndAttributes( allBuffer, "tool", "id=\""+getLegalName( toolname )+"\" name=\""+toolname+"\" version=\""+version+"\" force_history_refresh=\"true\"" );
		return allBuffer.toString();
	}
	
	/**
	 * Parses the values of the parameters from a galaxy script file
	 * @param filename the name of the script file
	 * @throws Exception if the parameter values could not be parsed
	 */
	public void fromGalaxyConfig(String filename) throws Exception{
		StringBuffer buf = new StringBuffer();
		BufferedReader read = new BufferedReader( new FileReader( filename ) );
		String tmp = null;
		while( (tmp = read.readLine()) != null){
			buf.append( tmp );
		}
		read.close();
		parameters.fromGalaxy( getLegalName( toolname ), buf );
	}
	
	private String getLROutput(ListResult lr) throws IOException{
		StringBuffer all = new StringBuffer();
		Result r;
		if (lr.getAnnotation() != null) {
			ResultSet annotation = lr.getAnnotation();
			DataType d;
			for (int i = 0; i < annotation.getNumberOfResults(); i++) {
				r = annotation.getResultAt(i);
				d = r.getDatatype();
				if (d != DataType.PNG && d != DataType.HTML
						&& d != DataType.LIST && d != DataType.STORABLE) {
					StringBuffer sb = new StringBuffer();
					sb.append("<label>"+r.getName()+":</label>");
					sb.append(r.getValue().toString());
					XMLParser.addTags( sb, "div" );
					all.append( sb );
				}
			}
		}
		StringBuffer list = new StringBuffer();
			ResultSet[] res = lr.getValue();
			boolean newNames;
			int i = 0, j, k;
			for (; i < res.length; i++) {
				newNames = i == 0;
				k = res[i].getNumberOfResults() - 1;
				if (!newNames) {
					if (k + 1 != res[i - 1].getNumberOfResults()) {
						newNames = true;
					} else {
						for (j = 0; j <= k; j++) {
							if (!res[i].getResultAt(j).getName().equals(
									res[i - 1].getResultAt(j).getName())) {
								newNames = true;
								break;
							}
						}
					}
				}
				if (newNames) {
					if(i != 0){
						list.append( "</table>" );
					}
					list.append( "<table border=\"1\">" );
					list.append( "<tr>" );
					for (j = 0; j <= k; j++) {
						list.append( "<th>"+res[i].getResultAt(j).getName()+"</th>" );
					}
					list.append( "</tr>" );
				}
				// write results
				list.append( "<tr>" );
				for (j = 0; j <= k; j++) {
					if(res[i].getResultAt(j) instanceof SimpleResult){
						if(res[i].getResultAt( j ).getDatatype() == DataType.DOUBLE || res[i].getResultAt( j ).getDatatype() == DataType.FLOAT){
							double d = ((Number)res[i].getResultAt( j ).getValue()).doubleValue();
							if(Math.abs( d )<0.01 && d != 0){
								list.append( "<td>"+expFormat.format( d )+"</td>"  );
							}else{	
								list.append( "<td>"+format.format( d )+"</td>" );
							}
						}else{
							list.append("<td>"+res[i].getResultAt(j).getValue()+"</td>");
						}
					}else{
						list.append("<td>"+getOutput( res[i].getResultAt(j) )+"</td>");
					}
				}
				list.append( "</tr>" );
			}
			list.append( "</table>" );
			XMLParser.addTags( list, "div" );
			all.append( list );
			return all.toString();
	}
	
	private String getOutput(Result res) throws IOException{
		StringBuffer buf = new StringBuffer();
		StringBuffer temp2 = new StringBuffer();
		temp2.append( res.getName() );
		XMLParser.addTagsAndAttributes( temp2, "div", "class=\"head\"" );
		temp2.append( "<br />" );
		if(res instanceof SimpleResult){
			StringBuffer temp = new StringBuffer();
			temp.append( res.getValue().toString().replaceAll( "\\n", "<br />" ) );
			XMLParser.addTags( temp, "div" );
			buf.append( temp );
		}else if(res instanceof ListResult){
			buf.append( getLROutput( ((ListResult) res ) ) );
		}else if(res instanceof DataSetResult){
			buf.append( getDataSetOutput((DataSetResult)res) );
		}else if(res instanceof StorableResult){
			buf.append( getStorableOutput((StorableResult)res) );
		}else if(res instanceof DoubleTableResult){
			buf.append( getDTROutput((DoubleTableResult)res) );
		}else if(res instanceof ImageResult){
			buf.append( getIROutput((ImageResult)res) );
		}else if(res instanceof FileResult){
			buf.append( getFileOutput( (FileResult)res ) );
		}
		StringBuffer temp = new StringBuffer();
		temp.append( res.getComment() );
		if(res instanceof LinkedImageResult){
			temp.append("<br />Obtain &quot;"+((LinkedImageResult)res).getLink().getName()+"&quot; ("+((LinkedImageResult)res).getLink().getComment()+") by clicking on the image");
		}
		XMLParser.addTagsAndAttributes( temp, "div", "class=\"comment\"" );
		temp.append( "<br />" );
		buf.append( temp );
		//XMLParser.addTagsAndAttributes( buf, "div", "class=\"form-row\"" );
		//XMLParser.addTagsAndAttributes( buf, "div", "class=\"toolFormBody\"" );
		temp2.append( buf );
	//	XMLParser.addTagsAndAttributes( temp2, "div", "class=\"toolForm\"" );
		return temp2.toString();
	}

	private String getDTROutput( DoubleTableResult res ) {
		double[][] r = res.getValue();
		StringBuffer sb = new StringBuffer();
		sb.append( "<table border=\"1\">" );
		for(int i=0;i<r.length;i++){
			sb.append( "<tr>" );
			for(int j=0;j<r[i].length;j++){
				sb.append( "<td>"+r[i][j]+"</td>" );
			}
			sb.append( "</tr>" );
		}
		sb.append( "</table>" );
		return sb.toString();
	}

	private String getIROutput( ImageResult res ) throws IOException {
		String name = getLegalName( res.getName() )+getHtmlId()+".";
		
		String filename = htmlFilesPath+System.getProperty( "file.separator" )+name;
		File f = new File(filename+"png");
		f.getParentFile().mkdirs();
		BufferedImage img = ((ImageResult)res).getValue();
		ImageIO.write( img, "png", f );
		String ext = "png";		
		
		StringBuffer sb = new StringBuffer();
		sb.append( "<img src=\""+name+ext+"\" alt=\""+res.getName()+"\" width=\""+img.getWidth()+"\" height=\""+img.getHeight()+"\"/>" );
		if(res instanceof LinkedImageResult){
			FileResult fr = ((LinkedImageResult)res).getLink();
			XMLParser.addTagsAndAttributes( sb, "a", "href=\""+fr.getFilename()+"."+fr.getExtension()+"\"" );
		}
		return sb.toString();
	}
	
	private String getFileOutput( FileResult res ){
		return "<a href=\""+res.getFilename()+"."+res.getExtension()+"\">"+res.getName()+"</a>";
	}
	
	private String getStorableOutput( StorableResult res ) throws IOException {
		String name = getLegalName( res.getName() )+getHtmlId()+".";
		String ext = export( htmlFilesPath+System.getProperty( "file.separator" )+name, res, null );
		return "<a href=\""+name+ext+"\">"+res.getName()+"</a>";
	}

	private String getDataSetOutput( DataSetResult res ) throws IOException {
		DataSet data = res.getValue();
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		if(res.getParser() == null){
			data.save( baos,'>', new SplitSequenceAnnotationParser( ":", ";" ) );
		}else{
			data.save( baos,'>', res.getParser() );
		}
		return baos.toString().replaceAll( "\\n", "<br />" );
	}

	private String getOutput(ResultSet res) throws IOException{
		Result[] ress = null;
		if(res instanceof MeanResultSet){
			ress = ((MeanResultSet)res).getStatistics().getResults();
		}else{
			ress = res.getResults();
		}
		StringBuffer sb = new StringBuffer();
		for(int i=0;i<ress.length;i++){
			sb.append( getOutput( ress[i] ) );
		}
		return sb.toString();
	}
	
	/**
	 * Exports a specified {@link Result} of a program execution
	 * to a file provided by <code>filename</code> and returns the
	 * corresponding Galaxy data type.
	 * @param filename the filename
	 * @param res the result
	 * @param exportExtension the extension used for the exported file
	 * @return the data type
	 * @throws IOException if the contents of <code>res</code> could not be written to the file
	 */
	public String export(String filename, Result res, String exportExtension) throws IOException{
		String ee = exportExtension;
		if(res instanceof SimpleResult){
			if(ee == null){
				ee = "txt";
			}
			File f = new File(filename+ee);
			f.getParentFile().mkdirs();
			PrintWriter pw = new PrintWriter( filename+ee );
			pw.println(res.toString());
			pw.close();
			return ee;
		}else if(res instanceof ListResult){
			if(ee == null){
				ee = "tabular";
			}
			File f = new File(filename+ee);
			f.getParentFile().mkdirs();
			PrintWriter pw = new PrintWriter( filename+ee );
			if(ee.equalsIgnoreCase( "gff3" )){
				pw.println("##gff-version 3");
			}else if(ee.equalsIgnoreCase( "gff" )){
				pw.println("##gff-version 2");
			}
			((ListResult)res).print( pw );
			pw.close();
			return ee;
		}else if(res instanceof DataSetResult){
			if(ee == null){
				ee = "fasta";
			}
			File f = new File(filename+ee);
			f.getParentFile().mkdirs();
			FileOutputStream fos = new FileOutputStream( filename+ee );
			if( ((DataSetResult)res).getParser() == null ){
				((DataSetResult)res).getValue().save(fos,'>',new SplitSequenceAnnotationParser( ":", ";" ) );
			}else{
				((DataSetResult)res).getValue().save(fos,'>', ((DataSetResult)res).getParser() );
			}
			fos.close();
			return ee;
		}else if(res instanceof StorableResult){
			if(ee == null){
				ee = "xml";
			}
			File f = new File(filename+ee);
			f.getParentFile().mkdirs();
			PrintWriter pw = new PrintWriter( filename+ee );
			pw.println(((StorableResult)res).getValue());
			pw.close();
			return ee;
		}else if(res instanceof DoubleTableResult){
			if(ee == null){
				ee = "tabular";
			}
			File f = new File(filename+ee);
			f.getParentFile().mkdirs();
			PrintWriter pw = new PrintWriter( filename+ee );
			double[][] tab = ((DoubleTableResult)res).getValue();
			for(int i=0;i<tab.length;i++){
				for(int j=0;j<tab[i].length-1;j++){
					pw.print( tab[i][j]+"\t" );
				}
				if(tab[i].length > 0){
					pw.println(tab[i][tab[i].length-1]);
				}else{
					pw.println();
				}
			}
			pw.close();
			return ee;
		}else if(res instanceof LinkedImageResult){
			return export( filename, ((LinkedImageResult)res).getLink(), exportExtension );
		}else if(res instanceof ImageResult){
			if(ee == null){
				ee = "png";
			}
			File f = new File(filename+ee);
			f.getParentFile().mkdirs();
			BufferedImage img = ((ImageResult)res).getValue();
			ImageIO.write( img, ee, new File(filename+ee) );
			return ee;
		}else if(res instanceof FileResult){
			String ext = ((FileResult)res).getExtension();
			File f = new File(filename+ext);
			f.getParentFile().mkdirs();
			FileManager.copy( ((FileResult)res).getValue().getAbsolutePath(), f.getAbsolutePath() );
			//System.out.println("exported "+f.getAbsolutePath());
			return ext;
		}
		return null;
	}
	
	/**
	 * Writes all output files of one program execution. This includes
	 * the summary of the execution, exported results and protocols, and images.
	 * @throws IOException if any of the files could not be created or written
	 */
	public void writeOutput() throws IOException{
		StringBuffer summary = new StringBuffer();
		
		int i=0;
		boolean exported = false;
		for(OutputElement el : list){
			boolean export = el.export;
			Object res = el.result;
			boolean include = el.includeInSummary;
			
			String str = null;
			if(include){
				if(res instanceof Result){
					str = getOutput( (Result) res );
				}else if(res instanceof ResultSet){
					str = getOutput( (ResultSet) res );
				}else{
					str = res.toString().replaceAll( "\\n", "<br />" );
				}
				summary.append( str );
			}
			if(export){
				exported = true;
				if(res instanceof Result){
					i++;
					String name = i+": "+((Result)res).getName();
					String ext = export( newFilePath+System.getProperty( "file.separator" )+"primary_"+outfileId+"_"+name+"_visible_", (Result)res, el.exportExtension );
				}else{
					ResultSet rs = (ResultSet)res;
					for(int j=0;j<rs.getNumberOfResults();j++){
						i++;
						String name = i+": "+rs.getResultAt( j ).getName();
						String ext = export( newFilePath+System.getProperty( "file.separator" )+"primary_"+outfileId+"_"+name+"_visible_", rs.getResultAt( j ), el.exportExtension );
					}
				}
			}
		}
		
		if(protocol != null){
			CategoricalResult prot = new CategoricalResult( "Protocol", "The protocol of this "+toolname+" run", protocol.toString() );
			summary.append( getOutput( prot ) );
			if(exportProtocol){
				i++;
				export( newFilePath+System.getProperty( "file.separator" )+"primary_"+outfileId+"_"+i+"_visible_", prot, null );
			}
		}
		XMLParser.addTags( summary, "body" );
		
		StringBuffer all = new StringBuffer();
		StringBuffer head = new StringBuffer();
		head.append( "Summary of "+toolname+" results" );
		XMLParser.addTags( head, "title" );
		all.append( head );
		//head = new StringBuffer();
		//XMLParser.addTagsAndAttributes( head, "link", "href=\""+stylesheetURL +"\" rel=\"stylesheet\" type=\"text/css\"" );
		//all.append( head );
		all.append( stylesheet );
		XMLParser.addTags( all, "head" );
		all.append( "<h2>"+head+"</h2>" );
		all.append( summary );
		
		XMLParser.addTags( all, "html" );
		
		PrintWriter wr = new PrintWriter( outfile );
		wr.println(all);
		
		wr.close();
		
	}
	
	/**
	 * Parses the command line. If the first argument is equal to &quot;--create&quot;, then
	 * the configuration file is written to the file with filename provided in the second argument.
	 * Otherwise the first argument must be &quot;--run&quot; to parse the parameters of a program
	 * run within Galaxy.
	 * @param args the arguments
	 * @return <code>true</code> if this execution should be a program run (as opposed to writing a configuration file)
	 * @throws Exception if the arguments could not be parsed or the Galaxy configuration file could not be created
	 */
	public boolean parse(String[] args) throws Exception{
		if("--create".equals( args[0] )){
			String str = toGalaxyConfig();
			PrintWriter wr = new PrintWriter( args[1] );
			wr.println(str);
			wr.close();
			return false;
		}else if("--run".equals( args[0] )){
			fromGalaxyConfig( args[1] );
			outfile = args[2];
			outfileId = args[3];
			newFilePath = args[4];
			htmlFilesPath = args[5];
			return true;
		}
		return false;
	}
	
	/**
	 * Adds a result to the results of a program run.
	 * @param res the result
	 * @param export if <code>true</code> the result is exported to its own Galaxy result, e.g. for
	 * 				evaluation in other application within Galaxy
	 * @param includeInSummary if <code>true</code> the result is shown on the summary page
	 */
	public void addResult(Result res, boolean export, boolean includeInSummary){
		list.add( new OutputElement( res, export, includeInSummary, null ) );
	}
	
	/**
	 * Adds a result to the results of a program run.
	 * @param res the result
	 * @param export if <code>true</code> the result is exported to its own Galaxy result, e.g. for
	 * 				evaluation in other application within Galaxy
	 * @param includeInSummary if <code>true</code> the result is shown on the summary page
	 * @param exportExtension the file extension used for the exported file
	 */
	public void addResult(Result res, boolean export, boolean includeInSummary, String exportExtension){
		list.add( new OutputElement( res, export, includeInSummary, exportExtension ) );
	}
	
	/**
	 * Adds a set of results to the results of a program run.
	 * @param res the results
	 * @param exportAll if <code>true</code> all results in this set are exported to their own Galaxy result, e.g. for
	 * 				evaluation in other application within Galaxy
	 * @param includeInSummary if <code>true</code> the results are shown on the summary page
	 */
	public void addResultSet(ResultSet res, boolean exportAll, boolean includeInSummary){
		list.add( new OutputElement( res, exportAll, includeInSummary, null ) );
	}
	
	/**
	 * Returns a legal variable name in Galaxy
	 * @param name the original name
	 * @return the legalized name
	 */
	public static String getLegalName(String name){
		return name.replaceAll( "[\\s,-:]+", "_" );
	}
	
	private static StringBuffer escape(StringBuffer str){
		return new StringBuffer( str.toString().replaceAll( "<", "&lt;" ).replaceAll( ">", "&gt;" ));
	}
	
	
	/**
	 * Class for an {@link ImageResult} that is linked to a file that can be downloaded.
	 * @author Jan Grau
	 *
	 */
	public static class LinkedImageResult extends ImageResult{

		private FileResult link;
		
		/**
		 * Create a new {@link ImageResult} with linked {@link FileResult} <code>link</code>
		 * @param name the name of the result
		 * @param comment the comment for the result
		 * @param image the image
		 * @param link the linked file
		 */
		public LinkedImageResult( String name, String comment, BufferedImage image, FileResult link ) {
			super( name, comment, image );
			this.link = link;
		}

		/**
		 * Creates a new {@link LinkedImageResult} from its XML-representation
		 * @param xml the representation
		 * @throws NonParsableException if xml could not be parsed
		 */
		public LinkedImageResult( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}

		/*
		 * (non-Javadoc)
		 * @see de.jstacs.results.ImageResult#getXMLTag()
		 */
		@Override
		public String getXMLTag() {
			return getClass().getSimpleName();
		}

		/*
		 * (non-Javadoc)
		 * @see de.jstacs.results.ImageResult#appendFurtherInfos(java.lang.StringBuffer)
		 */
		@Override
		protected void appendFurtherInfos( StringBuffer sb ) {
			XMLParser.appendObjectWithTags( sb, link, "link" );
		}

		/*
		 * (non-Javadoc)
		 * @see de.jstacs.results.ImageResult#extractFurtherInfos(java.lang.StringBuffer)
		 */
		@Override
		protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
			link = XMLParser.extractObjectForTags( representation, "link", FileResult.class );
		}
		
		/**
		 * Returns the linked file
		 * @return the file
		 */
		public FileResult getLink(){
			return link;
		}
		
	}
	
	/**
	 * {@link Result} for files that are results of some computation. Also used to link
	 * e.g. PDFs of images to {@link LinkedImageResult}s.
	 * @author Jan Grau
	 *
	 */
	public static class FileResult extends Result{

		private String path;
		private String filename;
		private String extension;
		
		/**
		 * Creates a new {@link FileResult} with name, comment, and path to the file.
		 * @param name the name of the result
		 * @param comment the comment for the result
		 * @param fullPath the path to the file
		 */
		public FileResult(String name, String comment, String fullPath){
			super(name, comment, DataType.FILE);
			int idx = fullPath.lastIndexOf( System.getProperty( "file.separator" ) );
			int extIdx = fullPath.lastIndexOf( "." );
			this.path = fullPath.substring( 0, idx );
			this.filename = fullPath.substring( idx+1, extIdx );
			this.extension = fullPath.substring( extIdx+1 );
			this.getValue().getParentFile().mkdirs();
		}
		
		/**
		 * Creates a new {@link FileResult} with name, comment, path to the file, filename and extension.
		 * @param name the name of the result
		 * @param comment the comment for the result
		 * @param path the path to the directory containing the file
		 * @param filename the filename without extension
		 * @param extension the filename extension
		 */
		public FileResult(String name, String comment, String path, String filename, String extension){
			super(name, comment, DataType.FILE);
			this.path = path;
			this.filename = filename;
			this.extension = extension;
			this.getValue().getParentFile().mkdirs();
		}
		
		
		
		/**
		 * Creates a new {@link FileResult} from its XML-representation
		 * @param rep the representation
		 * @throws NonParsableException if pre could not be parsed
		 */
		public FileResult( StringBuffer rep ) throws NonParsableException {
			super( rep );
		}
		
		/*
		 * (non-Javadoc)
		 * @see de.jstacs.results.ImageResult#getXMLTag()
		 */
		@Override
		public String getXMLTag() {
			return getClass().getSimpleName();
		}

		/*
		 * (non-Javadoc)
		 * @see de.jstacs.AnnotatedEntity#appendFurtherInfos(java.lang.StringBuffer)
		 */
		@Override
		protected void appendFurtherInfos( StringBuffer sb ) {
			XMLParser.appendObjectWithTags( sb, path, "path" );
			XMLParser.appendObjectWithTags( sb, filename, "filename" );
			XMLParser.appendObjectWithTags( sb, extension, "extension" );
		}

		/*
		 * (non-Javadoc)
		 * @see de.jstacs.AnnotatedEntity#extractFurtherInfos(java.lang.StringBuffer)
		 */
		@Override
		protected void extractFurtherInfos( StringBuffer rep ) throws NonParsableException {
			path = XMLParser.extractObjectForTags( rep, "path", String.class );
			filename = XMLParser.extractObjectForTags( rep, "filename", String.class );
			extension = XMLParser.extractObjectForTags( rep, "extension", String.class );
		}

		@Override
		public File getValue() {
			String sep = System.getProperty( "file.separator" );
			return new File(path+sep+filename+"."+extension);
		}

		/**
		 * Returns the path of the directory containing the file
		 * @return the path
		 */
		public String getPath() {
			return path;
		}

		/**
		 * Sets the path of the directory containing the file to <code>path</code>
		 * @param path the new path
		 */
		public void setPath( String path ) {
			this.path = path;
		}

		/**
		 * Returns the filename. 
		 * @return the filename
		 */
		public String getFilename() {
			return filename;
		}

		/**
		 * Sets the file
		 * @param filename the new filename
		 */
		public void setFilename( String filename ) {
			this.filename = filename;
		}

		/**
		 * Returns the filename extension
		 * @return the extension
		 */
		public String getExtension() {
			return extension;
		}

		/**
		 * Sets the filename extension
		 * @param extension the new extension
		 */
		public void setExtension( String extension ) {
			this.extension = extension;
		}
		
		
		
	}
	
	/**
	 * Returns the path where files, e.g. images, that shall be linked from the HTML summary
	 * shall be stored in.
	 * @return the path
	 */
	public String getHtmlFilesPath() {
		return htmlFilesPath;
	}
	
	/**
	 * Class for a Protocol writer.
	 * @author Jan Grau
	 *
	 */
	public static class Protocol{

		private ByteArrayOutputStream baos;
		private PrintWriter wr;
		
		/**
		 * @param out
		 */
		private Protocol( ) {
			baos = new ByteArrayOutputStream();
			wr = new PrintWriter( baos );
		}
		
		/**
		 * Appends <code>str</code> to the protocol.
		 * @param str the string to be appended
		 */
		public void append(String str){
			wr.append( str );
		}
		
		/**
		 * Returns the {@link PrintWriter} of this protocol
		 * @return the writer
		 */
		public PrintWriter getWriter(){
			return wr;
		}
		
		/**
		 * Returns the {@link ByteArrayOutputStream} of this protocol
		 * @return the stream
		 */
		public ByteArrayOutputStream getOutputStream(){
			return baos;
		}
		
		/**
		 * Append a heading to the protocol
		 * @param str the title of the heading
		 */
		public void appendHeading(String str){
			wr.append( "<strong>"+str+"</strong>\n" );
			wr.flush();
		}
		
		/**
		 * Appends a warning to the protocol
		 * @param str the warning
		 */
		public void appendWarning(String str){
			wr.append( "<em>"+str+"</em>\n" );
			wr.flush();
		}
		
		/**
		 * Returns the current version of the protocol as {@link String}.
		 * @return the protocol as {@link String}
		 */
		public String toString(){
			wr.flush();
			return baos.toString().replaceAll( "\n", "<br />" );
		}
		
	}
	
	private static class OutputElement{
		private Object result;
		private boolean export;
		private boolean includeInSummary;
		private String exportExtension;
		
		/**
		 * @param result
		 * @param export
		 * @param includeInSummary
		 */
		private OutputElement( Object result, boolean export, boolean includeInSummary, String exportExtension ) {
			this.result = result;
			this.export = export;
			this.includeInSummary = includeInSummary;
			this.exportExtension = exportExtension;
		}
		
		
		
	}
	
}
