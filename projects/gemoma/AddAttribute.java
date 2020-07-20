package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;

import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.FileExistsValidator;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;

/**
 * Allows to add an attribute to features of an annotation file.
 * 
 * @author Jens Keilwagen
 */
public class AddAttribute extends GeMoMaModule {
	public File addAttribute(ToolParameterSet parameters, Protocol protocol, String tag, String annotation, String extraFileName, int idColumn, int valColumn, String attribute ) throws Exception {
		BufferedReader r = Tools.openGzOrPlain(extraFileName);
		HashMap<String,HashSet<String>> ids = new HashMap<String,HashSet<String>>();
		String line;
		HashSet<String> dummy = new HashSet<String>(); 
		while( (line=r.readLine()) != null ) {
			String[] split = line.split("\t");
			String key = split[idColumn];
			HashSet<String> vals = ids.get(key);
			if( vals == null ) {
				vals = (valColumn<0 ? dummy : new HashSet<String>());
				ids.put(key, vals);
			}
			if( valColumn >= 0 ) {
				vals.add(split[valColumn]);
			}
		}
		r.close();

		File output = Tools.createTempFile("AddAttribute");
		BufferedWriter w = new BufferedWriter( new FileWriter(output) );
		r = Tools.openGzOrPlain(annotation);
		int[] a = new int[2];
		String[] empty = new String[0];
		boolean first=true;
		while( (line=r.readLine()) != null ) {
			if( line.length()==0 || line.charAt(0) == '#' ) {
				w.write(line);
			} else {
				if( first ) {
					first=false;
					w.append(GeMoMa.INFO + getShortName() + " " + getToolVersion() + "; ");
					String info = JstacsTool.getSimpleParameterInfo(parameters);
					if( info != null ) {
						w.append("SIMPLE PARAMETERS: " + info );
					}
					w.newLine();
				}
				String[] split = line.split("\t");
				if( split[2].equals(tag) ) {
					String temp = split[split.length-1];
					int idx1 = temp.indexOf("ID=")+3;
					int idx2 = temp.indexOf(';',idx1);
					String id = temp.substring(idx1,idx2);
					HashSet<String> c = ids.get(id);
					temp += (temp.charAt(temp.length()-1)==';'?"":";");
					if( valColumn<0 ) {
						temp += attribute + "=" + (c!=null);
						a[c!=null?0:1]++;
					} else {
						if( c!=null ) {
							temp += attribute + "=";
							String[] v = c.toArray(empty);
							if( v.length > 1 ) temp+="\"";
							for( int i = 0; i < v.length; i++ ) temp += (i==0?"":", ") + v[i];
							if( v.length > 1 ) temp+="\"";
							a[0]++;
						} else {
							a[1]++;
						}
					}
					split[split.length-1] = temp;
					for( int i = 0; i < split.length; i++ ) {
						w.write((i==0?"":"\t") + split[i]);
					}
				} else {
					w.write(line);
				}
			}
			w.newLine();
		}
		r.close();
		w.close();
		
		protocol.append( "items with additional attribute " + attribute + " :\t" + a[0] + "\n");
		protocol.append( "items without additional attribute " + attribute + " :\t" + a[1] + "\n");
		return output;
	}

	@Override
	public ToolParameterSet getToolParameters() {
		try {
			return  new ToolParameterSet( getShortName(),
					new FileParameter( "annotation", "annotation file", "gff,gff3", true, new FileExistsValidator(), true ),
					new SimpleParameter( DataType.STRING, "feature", "a feature of the annotation, e.g., gene, transcript or mRNA", true, GeMoMa.TAG ),
					new SimpleParameter( DataType.STRING, "attribute", "the name of the attribute that is added to the annotation", true ),
					new FileParameter( "table", "a tab-delimited file containing IDs and additional attribute", "tabular", true, new FileExistsValidator(), false ),
					new SimpleParameter( DataType.INT, "ID column", "the ID column in the tab-delimited file", true, new NumberValidator<Integer>(0, Integer.MAX_VALUE) ),
					new SelectionParameter( DataType.PARAMETERSET, 
							new String[]{"VALUES","BINARY"},
							new Object[]{
									new SimpleParameterSet(
											new SimpleParameter( DataType.INT, "attribute column", "the attribute column in the tab-delimited file", true, new NumberValidator<Integer>(0, Integer.MAX_VALUE) )
									),
									new SimpleParameterSet()
							},
							"type", "type of addition attribute", true)
			);
		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		SimpleParameterSet sps = (SimpleParameterSet) parameters.getParameterForName("type").getValue();
		File output = addAttribute( parameters, protocol,
				(String) parameters.getParameterForName("feature").getValue(),
				(String) parameters.getParameterForName("annotation").getValue(),
				(String) parameters.getParameterForName("table").getValue(),
				(Integer) parameters.getParameterForName("ID column").getValue(),
				sps.getNumberOfParameters()==0 ? -1 : (Integer) sps.getParameterAt(0).getValue(),
				(String) parameters.getParameterForName("attribute").getValue()
		);
		
		return new ToolResult("", "", null, new ResultSet(new TextResult("extended annotation", "Result", new FileRepresentation(output.getAbsolutePath()), "gff", getToolName(), null, true)), parameters, getToolName(), new Date());
	}

	@Override
	public String getToolName() {
		return "AddAttribute";
	}

	@Override
	public String getShortName() {
		return getToolName();
	}

	@Override
	public String getDescription() {
		return "adds an attribute to an annotation file";
	}

	@Override
	public String getHelpText() {
		return 
			"This tool allows to add an additional attribute to specific features of an annotation.\n\n"
			+ "Those additional attributes might be used in **GAF** for filtering or sorting or might be displayed in genome browsers like IGV or WebApollo. "
			+ "The user can choose binary attributes (true or false) or attributes with values according to given tab-delimited table."
			+ MORE;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", AnnotationFinalizer.defResult ),
		};
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		// TODO Auto-generated method stub
		return null;
	}	
}