package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import de.jstacs.DataType;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.FileExistsValidator;
import de.jstacs.parameters.validation.RegExpValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
//import projects.gemoma.old.GeMoMa;

public class GAFComparison extends GeMoMaModule {

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads,
			String temp) throws IOException {
		String tag = parameters.getParameterForName("tag").getValue().toString();
		boolean splitPrefix = (Boolean)parameters.getParameterForName("split prefix").getValue();
		boolean differences = (Boolean)parameters.getParameterForName("differences").getValue();
		
		HashMap<String,HashMap<String,int[]>[]> hash = new HashMap<String,HashMap<String,int[]>[]>();
		
		//read genes in GAF gff format!
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterForName("predicted annotation").getValue();
		
		BufferedReader r;
		int MAX = eps.getNumberOfParameters();
		String[] name = new String[MAX];
		for( int k = 0; k < MAX; k++ ) {
			SimpleParameterSet sps = ((SimpleParameterSet)eps.getParameterAt(k).getValue());
			name[k] = (String) sps.getParameterForName("name").getValue();
			protocol.append("reading gene annotation: " + name[k] + "\n");
			r = Tools.openGzOrPlain( (String) sps.getParameterForName("gene annotation file").getValue());
			String line;
			while( (line=r.readLine()) != null ) {
				if( line.length()== 0 || line.charAt(0)=='#') continue;
				String[] split = line.split("\t");
				if( split[2].equals(tag) ) {
					split = split[8].split(";");
					String pa=null;
					for( int j = 0; j < split.length; j++ ) {
						if( split[j].startsWith("Parent=") ) {
							pa=split[j].substring(7);
						}
					}
					for( int j = 0; j < split.length; j++ ) {
						if( split[j].startsWith("ref-gene=") ) {
							add(hash,pa,split[j].substring(9),name.length,k);
						} else if ( split[j].startsWith("alternative=\"") ) {
							String[] s = split[j].substring(13,split[j].length()-1).split(",");
							for( int h = 0; h < s.length; h++ ) {
								add(hash,pa,s[h],name.length,k);
							}
						}
					}					
				}
			}
			r.close();
		}

		String[] genes = hash.keySet().toArray(new String[0]);
		Arrays.sort(genes);

		File out = Tools.createTempFile("GAFComparison", temp);
		
		BufferedWriter w = new BufferedWriter( new FileWriter(out) );
		w.append((splitPrefix?"\t":"") + "\t\t");
		for( int i = 0; i < name.length; i++ ) {
			w.append( "\t" + name[i] + "\t" );
		}
		w.newLine();
		if( splitPrefix ) System.out.print("prefix\t");
		w.append("gene-ID\tmin #genes\tmax #genes");
		for( int i = 0; i < name.length; i++ ) {
			w.append( "\t#genes\t#transcripts" );
		}
		w.newLine();
		for( String g : genes ) {
			HashMap<String,int[]>[] val = hash.get(g);
			int[] v = summary(val, name.length);
			
			int max = 0, min = Integer.MAX_VALUE;
			for( int i = 0; i < name.length; i++ ) {
				if( v[2*i]>max )  max = v[2*i];
				if( v[2*i]<min )  min = v[2*i];
			}
			
			if( !differences || min < max ) {
				int idx = g.indexOf('_');
				if( idx>=0 && splitPrefix ) {
					w.append(g.substring(0, idx) + "\t" + g.substring(idx+1) ); 
				} else {
					w.append(g);
				}
				w.append( "\t" + min );
				w.append( "\t" + max );
				for( int i = 0; i < v.length; i++ ) {
					w.append( "\t" + v[i] );
				}
				w.newLine();
			}
		}
		w.close();
		
		return new ToolResult("", "", null, new ResultSet(new TextResult("GAF comparison", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "tabular", getToolName(), null, true)), parameters, getToolName(), new Date());
	}
	
	static void add( HashMap<String, HashMap<String,int[]>[]> hash, String parent, String rGene, int max, int current ) {
		rGene=rGene.toUpperCase();
		
		HashMap<String,int[]>[] val = hash.get(rGene);
		if( val == null ) {
			val = new HashMap[max];
			for( int i = 0; i < max; i++ ) {
				val[i] = new HashMap<String,int[]>();
			}
			hash.put(rGene, val);
		}
		int[] v = val[current].get(parent);
		if( v == null ) {
			v = new int[1];
			val[current].put( parent, v );
		}
		v[0]++;
	}
	
	static int[] summary( HashMap<String,int[]>[] val, int max ) {
		int[] res = new int[2*max];
		for( int i = 0; i < max; i++ ) {
			Iterator<Entry<String,int[]>> it = val[i].entrySet().iterator();
			while( it.hasNext() ) {
				Entry<String,int[]> e = it.next();
				res[2*i]++;
				res[2*i+1]+=e.getValue()[0];
			}
		}
		return res;
	}

	@Override
	public ToolParameterSet getToolParameters() {
		try{
			return
				new ToolParameterSet( getShortName(),
					new SimpleParameter(DataType.STRING,"tag","the tag used to read the GAF annotations",true,GeMoMa.TAG),
					new ParameterSetContainer( "predicted annotation", "", new ExpandableParameterSet( new SimpleParameterSet(
							new SimpleParameter(DataType.STRING,"name","a simple name for the organism", true, new RegExpValidator(".+") ),
							new FileParameter( "gene annotation file", "GFF file containing the gene annotations (predicted by GAF)", "gff,gff3,gff.gz.gff3.gz",  true, new FileExistsValidator(), true )
						), "gene annotations", "", 1 ) ),
					new SimpleParameter(DataType.BOOLEAN, "split prefix", "a switch to decide whether the prefix should be split and writen in a separat column", true, false),
					new SimpleParameter(DataType.BOOLEAN, "differences", "a switch to decide whether only genes with difference should be returned", true, true)
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
		return "GAFComparison";
	}

	@Override
	public String getDescription() {
		return "creates stats for the GFF attributes ref-gene and alternative";
	}

	@Override
	public String getHelpText() {
		return
				"This tool allows to compare results from GAF based on the attributed ref-gene and alternative.\n"			
				+ "Hence, you can compare the annotation of different genomes or the effect of different parameters on the annotation of one genome."
				+ MORE;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "tabular", "GAF comparison"),
		};
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		// TODO Auto-generated method stub
		return null;
	}
}
