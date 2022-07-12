package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
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

	static String[] pos = {"0chromosome","1start position","2end postion"};
	static int[] ind = {0, 3, 4};
	
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp) throws Exception {
		//read and parse annotation
		String missing = parameters.getParameterForName("missing").getValue().toString();
		String tag = parameters.getParameterForName("feature").getValue().toString();
		String fName = parameters.getParameterForName("annotation").getValue().toString();
		
		HashMap<String, Integer> keys = new HashMap<String,Integer>();
		for( int i = 0; i < pos.length; i++ ) {
			keys.put(pos[i], i);
		}
		ArrayList<String[]> feature = new ArrayList<String[]>();
		BufferedReader r = Tools.openGzOrPlain(fName);
		String line;
		while( (line=r.readLine()) != null ) {
			if( line.equalsIgnoreCase("##FASTA") ) break; //http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
			if( line.length() == 0 || line.startsWith("#") ) continue; 
			
			String[] split = line.split("\t");
			if( split[2].equalsIgnoreCase(tag) ) {
				//parse and extend keys if necessary
				String[] att = split[8].split(";");
				String[] v = new String[att.length+3];
				int[] index = new int[att.length+3];
				for( int i = 0; i < att.length; i++ ) {
					int j = att[i].indexOf('=');
					if( j < 0 ) throw new IllegalArgumentException("Could not parse attribute: " + att[i]);
					String a = att[i].substring(0,j);
					v[i] = att[i].substring(j+1);
					
					Integer ind = keys.get(a);
					if( ind == null ) {
						ind = keys.size();
						keys.put(a, ind );
					}
					index[i]=ind;
				}
				for( int i = 0; i < pos.length; i++ ) {
					index[att.length+i]=keys.get(pos[i]);
					v[att.length+i]=split[ind[i]];
				}
				
				//create entry
				String[] filled = new String[keys.size()];
				Arrays.fill(filled, missing);
				for( int i = 0; i < v.length; i++ ) {
					filled[index[i]] = v[i];
				}
				feature.add(filled);
			}
		}
		r.close();
		
		//output
		String[] k = keys.keySet().toArray(new String[0]);
		Arrays.sort(k);
		
		Collections.sort(feature,new StringArrayComparator() );
		
		File out = Tools.createTempFile("GFF-attributes", temp );
		BufferedWriter w = new BufferedWriter( new FileWriter(out) );
		for( int i = 0; i <k.length; i++ ) {
			String h =k[i];
			if( i<pos.length) {
				h=h.substring(1);
			}
			w.append( (i==0?"":"\t") + h );
		}
		w.newLine();
		for( int i = 0; i < feature.size(); i++ ) {
			String[] attributes = feature.get(i);
			
			for( int j = 0; j <k.length; j++ ) {
				int idx = keys.get(k[j]);
				w.append( (j==0?"":"\t") + (idx < attributes.length? attributes[idx] : missing) );
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
	
	public static class StringArrayComparator implements Comparator<String[]> {

		@Override
		public int compare(String[] a, String[] b) {
			int res = 0;
			if( a == b ) {}
			else if( a==null ) {
				res=1;
			} else if( b==null ) {
				res=-1;
			} else {
				int i = 0, n = Math.min(a.length,b.length);
				while( res == 0 && i < n ) {
					if( i==2 || i==1 ) {
						res = Integer.compare( Integer.parseInt(a[i]), Integer.parseInt(b[i]) );
					} else {
						res = a[i].compareTo(b[i]);
					}
					i++;
				}
			}
			return res;
		}
		
	}
}