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

package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;

import de.jstacs.DataType;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.FileExistsValidator;
import de.jstacs.parameters.validation.RegExpValidator;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import projects.gemoma.GeMoMaAnnotationFilter.Prediction;

/**
 * This class allows get attribute for all references.
 * 
 * @author Jens Keilwagen
 * 
 * @see GeMoMaAnnotationFilter
 */
public class Attribute2Table extends GeMoMaModule {
		
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String tempD) throws Exception {
		String tag = parameters.getParameterForName("tag").getValue().toString();
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterForName("raw annotation").getValue();
		boolean addPrefix = (Boolean) parameters.getParameterForName("add").getValue();
		String attribute = parameters.getParameterForName("attribute").getValue().toString();
		
		protocol.append( "\nused attribute: " + attribute + "\n\n" );
		
		String line, t;
		String[] split;
				
		Prediction current = null;
		HashMap<String,int[]> del = new HashMap<String, int[]>();
		
		//read genes in GeMoMa gff format!
		BufferedReader r;
		int MAX = eps.getNumberOfParameters();
		String[] prefix = new String[MAX];
		String[] sampleInfo = new String[MAX];
		HashMap<String,ArrayList<Prediction>> raw = new HashMap<String,ArrayList<Prediction>>();
		for( int k = 0; k < MAX; k++ ) {
			SimpleParameterSet sps = ((SimpleParameterSet)eps.getParameterAt(k).getValue());
			String h = (String) sps.getParameterAt(0).getValue();
			prefix[k] = h==null?"":h;
			String pref = addPrefix ? prefix[k] : "";
			sampleInfo[k] = k + (prefix[k]==null||prefix[k].length()==0?"":(" (" + prefix[k] + ")"));
			protocol.append("species " + sampleInfo[k] + "\n");
					
			//gff needs to be "clustered" (all features of one transcript should be in one block)
			String fName = sps.getParameterAt(1).getValue().toString();
			//System.out.println(fName);
			r = new BufferedReader(new FileReader(fName));
			while( (line=r.readLine()) != null ) {
				if( line.length() == 0 ) continue;
				if( line.charAt(0)=='#' ) {
					continue;
				}
				
				split = line.split("\t");
				
				t = split[2];
				if( t.equals("gene") ) {
					//TODO WARNING
				} else if( t.equals(tag) ) { //tag: prediction/mRNA/transcript
					current = new Prediction(split, MAX, k, pref, del );
					String key = current.hash.get("ref-gene"); //TODO add prefix?
					ArrayList<Prediction> list = raw.get( key );
					if( list == null ) {
						list = new ArrayList<Prediction>();
						raw.put(key, list);
					}
					list.add(current);
				} else {
					if( current==null ) {
						r.close();
						throw new NullPointerException("There is no gene model. Please check parameter \"tag\" and the order within your annotation file "+sampleInfo+": " + fName );
					}
					if( split[8].contains("Parent="+current.oldId) ) {
						if( t.equals( "CDS" ) ) current.addCDS( line );
					} else {
						//TODO WARNING
						r.close();
						throw new IllegalArgumentException("The GFF has to be clustered, i.e., all features of a transcript must be adjacent lines:\n" + line);
					}
				}
			}
			r.close();
		}
		
		File out = Tools.createTempFile("Attribute", tempD);
		BufferedWriter outW = new BufferedWriter( new FileWriter(out) );
		outW.append( "#geneID\tchr\tstart\tend" );
		for( int i = 0; i < sampleInfo.length; i++ ) {
			outW.append("\t" + sampleInfo[i] );
		}
		outW.newLine();
		
		String fName = parameters.getParameterForName("final gene annotation file").getValue().toString();
System.out.println( "\"" + fName + "\"\t" + ((new File(fName).exists())) );
		r = new BufferedReader(new FileReader(fName));
		double[] best = new double[MAX];
		current = null;
		while( (line=r.readLine()) != null ) {
			if( line.length() == 0 ) continue;
			if( line.charAt(0)=='#' ) {
				continue;
			}
			
			split = line.split("\t");
			
			t = split[2];
			if( t.equals("gene") ) {
				//TODO WARNING
			} else if( t.equals(tag) ) { //tag: prediction/mRNA/transcript
				if( current != null ) createLine( current, raw, attribute, outW, best );
				current = new Prediction(split, MAX, 0, "", del);
			} else {
				if( current==null ) {
					r.close();
					throw new NullPointerException("There is no gene model. Please check parameter \"tag\" and the order within your annotation file "+sampleInfo+": " + fName );
				}
				if( split[8].contains("Parent="+current.oldId) ) {
					if( t.equals( "CDS" ) ) current.addCDS( line );
				} else {
					//TODO WARNING
					r.close();
					throw new IllegalArgumentException("The GFF has to be clustered, i.e., all features of a transcript must be adjacent lines:\n" + line);
				}
			}
		}
		r.close();
		if( current != null ) createLine( current, raw, attribute, outW, best );
		outW.close();
		
		ArrayList<Result> res = new ArrayList<Result>();
		res.add( new TextResult("best attribute", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "tabular", getToolName(), null, true) );
		return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
	}
	
	void createLine( Prediction current, HashMap<String,ArrayList<Prediction>> raw, String attribute, BufferedWriter w, double[] best ) throws IOException {
		Arrays.fill(best, -1000);
		set( current, raw.get(current.hash.get("ref-gene")), attribute, best );
		String alt = current.hash.get("alternative");
		if( alt!= null && alt.length()>0 ) {
			if( alt.charAt(0)=='"' && alt.charAt(alt.length()-1)=='#' ) {
				alt = alt.substring(1,alt.length()-1);
			}
			String[] a = alt.split(",");
			for( String alter: a ) {
				set( current, raw.get(alter), attribute, best );
			}
		}
		w.append( current.id + "\t" + current.split[0] + "\t" + current.split[3] + "\t" + current.split[4] );
		for( int i = 0; i < best.length; i++ ) {
			w.append("\t" + ( best[i]<-100 ? "NA" : best[i] ) );
		}
		w.newLine();
	}
	
	void set( Prediction current, ArrayList<Prediction> p, String attribute, double[] best ) {
		if( p != null ) {
			for( Prediction b: p ) {
				if( current.compareTo(b) == 0 ) {
					int index = b.getIndex();
					double val = Double.parseDouble( b.hash.get(attribute) );
					if( val > best[index] ) best[index]=val;
				}
			}
		}
	}
	
	@Override
	public ToolParameterSet getToolParameters() {
		try{
			
			return
				new ToolParameterSet( getShortName(),
					new SimpleParameter(DataType.STRING,"tag","the tag used to read the GeMoMa annotations",true,GeMoMa.TAG),
					new FileParameter( "final gene annotation file", "GFF file containing the gene annotations (predicted by GeMoMa)", "gff,gff3",  true, new FileExistsValidator(), true ),
					new ParameterSetContainer( "raw annotation", "", new ExpandableParameterSet( new SimpleParameterSet(		
							new SimpleParameter(DataType.STRING,"prefix","the prefix can be used to distinguish predictions from different input files", false, new RegExpValidator("\\w*") ),
							new FileParameter( "raw gene annotation file", "GFF file containing the gene annotations (predicted by GeMoMa)", "gff,gff3",  true, new FileExistsValidator(), true )
					), "gene annotations", "", 1 ) ),
					new SimpleParameter(DataType.BOOLEAN,"add", "add the prefix to the gene ID", true, true ),
					new SimpleParameter(DataType.STRING,"attribute", "the attribute to be checked", true, "iAA" )
				);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public String getToolName() {
		return "Attribute2Table";
	}

	@Override
	public String getShortName() {
		return getToolName();
	}

	@Override
	public String getDescription() {
		return "returns a table of best attribute per predicted final annotation";
	}

	@Override
	public String getHelpText() {
		return	"This tool returns a table of best attribute per predicted final annotation."
				
				+ MORE;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}
	
	@Override
	public ToolResult[] getTestCases( String path ) {
		return null;
	}
}