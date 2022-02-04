package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;

import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.FileExistsValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;

/**
 * This class allows to a table of attributes from a specific class of features of a GFF.
 * 
 * @author Jens Keilwagen
 */
public class GFFAttributes extends GeMoMaModule {

	public GFFAttributes() {}

	@Override
	public ToolParameterSet getToolParameters() {
		try{
			return
				new ToolParameterSet( getShortName(),
					new FileParameter( "annotation", "GFF file containing the gene annotations", "gff,gff3",  true, new FileExistsValidator(), true ),
					new SimpleParameter(DataType.STRING,"feature","the feature which is used to parse the attributes",true,GeMoMa.TAG),
					new SimpleParameter(DataType.STRING,"missing","the value used for missing attributes of a feature",true,"")
				);
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public String getToolName() {
		return getShortName();
	}

	@Override
	public String getShortName() {
		return "GFFAttributes";
	}

	@Override
	public String getDescription() {
		return "extracts attributes from a feature class of GFF";
	}

	@Override
	public String getHelpText() {
		return "Annotations that are build with **GeMoMaPipeline** or augmented with **AnnotationEvidence** have lots of attributes that might be intersting for the user."
				+ " This module allows to create a simple table that can easily be parsed and used for visualization of statistics."
				+ " However, the module could also be used for annotations that are not created of modified with GeMoMa modules.";
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp) throws Exception {
		//read and parse annotation
		String missing = parameters.getParameterForName("missing").getValue().toString();
		String tag = parameters.getParameterForName("feature").getValue().toString();
		String fName = parameters.getParameterForName("annotation").getValue().toString();
		
		HashMap<String, Integer> keys = new HashMap<String,Integer>();
		HashMap<String, String[]> feature = new HashMap<String, String[]>();
		BufferedReader r = Tools.openGzOrPlain(fName);
		String line;
		while( (line=r.readLine()) != null ) {
			if( line.equalsIgnoreCase("##FASTA") ) break; //http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
			if( line.length() == 0 || line.startsWith("#") ) continue; 
			
			String[] split = line.split("\t");
			if( split[2].equalsIgnoreCase(tag) ) {
				//parse and extend keys if necessary
				String[] att = split[8].split(";");
				String[] v = new String[att.length];
				int[] index = new int[att.length];
				for( int i = 0; i < att.length; i++ ) {
					int j = att[i].indexOf('=');
					if( j < 0 ) throw new IllegalArgumentException("Could not parse attribute: " + att[i]);
					String a = att[i].substring(0,j);
					v[i] = att[i].substring(j+1);
					if( !a.equals("ID") ) {
						Integer ind = keys.get(a);
						if( ind == null ) {
							ind = keys.size();
							keys.put(a, ind );
						}
						index[i]=ind;
					} else {
						index[i]=-1;
					}
				}
				
				//create entry
				String[] filled = new String[keys.size()];
				Arrays.fill(filled, missing);
				String id = null;
				for( int i = 0; i < att.length; i++ ) {
					if( index[i]<0 ) {
						id=v[i];
					} else {
						filled[index[i]] = v[i];
					}
				}
				feature.put(id, filled);
			}
		}
		r.close();
		
		//output
		String[] k = keys.keySet().toArray(new String[0]);
		Arrays.sort(k);
		String[] id = feature.keySet().toArray(new String[0]);
		Arrays.sort(id);
		
		File out = Tools.createTempFile("GFF-attributes", temp );
		BufferedWriter w = new BufferedWriter( new FileWriter(out) );
		w.append("ID");
		for( String g: k ) {
			w.append( "\t" + g );
		}
		w.newLine();
		for( String i: id ) {
			String[] attributes = feature.get(i);
			
			w.append(i);
			for( String g: k ) {
				int idx = keys.get(g);
				w.append( "\t" + (idx < attributes.length? attributes[idx] : missing) );
			}
			w.newLine();
		}
		w.close();
		
		return new ToolResult("", "", null, new ResultSet(new TextResult("GFF attributes", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "tabular", getToolName(), null, true)), parameters, getToolName(), new Date());
	}
	
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "tabular", "GFF attributes"),
		};
	}
}