package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import de.jstacs.DataType;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
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

/**
 * 
 * @author Jens Keilwagen
 */
public class SyntenyChecker extends GeMoMaModule {


	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp) throws Exception {
			String tag = (String) parameters.getParameterForName("tag").getValue();

			//read assignment
			ArrayList<Assignment> assign = new ArrayList<Assignment>();
			ExpandableParameterSet eps = (ExpandableParameterSet) (parameters.getParameterForName("references").getValue());
			for( int i = 0; i < eps.getNumberOfParameters(); i++ ) {
				ParameterSet ps = (ParameterSet) eps.getParameterAt(i).getValue();
				String prefix = (String) ps.getParameterForName("prefix").getValue();
				if( prefix != null && !prefix.endsWith("_") ) {
					prefix += "_";
				}
				String assignFile = (String) ps.getParameterForName("assignment").getValue();
				assign.add( new Assignment( assignFile, prefix ) );
			}
			
			//read annotation
			String annot = (String) parameters.getParameterForName("gene annotation file").getValue();
			BufferedReader r = new BufferedReader( new FileReader( annot ) );
			HashMap<String, Gene> annotation = new HashMap<String, Gene>();
			String line;
			HashMap<String,String> a = new HashMap<String,String>();
			while( (line=r.readLine()) != null ) {
				if( line.length()==0 || line.charAt(0)=='#' ) continue;
				String[] split = line.split("\t");
				if( split[2].equals(tag) ) {
					String[] att = split[8].split(";");
					a.clear();
					for( int i = 0; i < att.length; i++ ) {
						int idx = att[i].indexOf('=');
						if( idx<0 ) continue;//TODO throw new RuntimeException("Could not parse attributes: " + att[i] + " from " + split[8]);
						String key = att[i].substring(0,idx);
						String value = att[i].substring(idx+1);
						if( value.charAt(0)=='"' ) {
							value = value.substring(1,value.length()-1);
						}
						a.put(key, value);
					}
					String parent = a.get("Parent");
					Gene g = annotation.get(parent);
					if( g == null ) {
						g = new Gene( parent );
						g.rg = new HashSet<String>();
						annotation.put(parent, g);
					}
					g.extend( split[0], split[6].charAt(0)=='+'?1:-1, Integer.parseInt(split[3]), Integer.parseInt(split[4]) );
	
					String refG = a.get("ref-gene");
					g.rg.add(refG);
					String altG = a.get("alternative");
					if( altG != null ) {
						String[] al = altG.split(",");
						for( int i = 0; i < al.length; i++ ) {
							g.rg.add(al[i]);
						}
					}
				}
			}
			r.close();
			
			//sort
			ArrayList<Gene> anno = new ArrayList<Gene>();
			Iterator<Gene> it = annotation.values().iterator();
			while( it.hasNext() ) {
				Gene g = it.next();
				anno.add(g);
			}
			Collections.sort(anno);
			
			//output
			String oldContig = null;
			int j=0;
			File output = Tools.createTempFile("synteny",temp);
			BufferedWriter w= new BufferedWriter( new FileWriter( output ) );
			w.append("contig\tmiddle position\tstrand\tord\tname\tref-genes");
			for( int as = 0; as < assign.size(); as++ ) {
				Assignment ass = assign.get(as);
				w.append( "\t" + (ass.prefix.length()>0?ass.prefix.substring(0,ass.prefix.length()-1):"") );
			}
			w.newLine();
			for( int i = 0; i < anno.size(); i++ ) {
				Gene g = anno.get(i);
				if( !g.chr.equals(oldContig) ) {
					j=0;
					oldContig=g.chr;
				}
				String target = g.chr +"\t" + g.middle;
				w.append(g.chr +"\t" + g.middle + "\t" + g.strand + "\t" + j++ + "\t" + g.name + "\t" + g.rg.size() );
				for( int as = 0; as < assign.size(); as++ ) {
					w.append( "\t" + assign.get(as).check( g, target ) );
				}
				w.newLine();
			}
			w.close();
			
			String out = output.getAbsolutePath();
			ArrayList<TextResult> res = new ArrayList<TextResult>();
			res.add( new TextResult("reference gene table", "Result", new FileParameter.FileRepresentation(out), "tabular", getToolName(), null, true) );
			return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
	}
	
	/**
	 * A class representing the assignment of a reference organism.
	 * 
	 * @author Jens Keilwagen
	 */
	static class Assignment {
		String prefix;
		HashMap<String,Gene> hash;
		
		public Assignment( String fName, String pref ) throws IOException {
			if( pref!=null ) {
				this.prefix = pref;
			} else {
				this.prefix = "";
			}
			hash = new HashMap<String, Gene>();

			//read reference (assignment)
			BufferedReader r = new BufferedReader( new FileReader( fName ) );
			r.readLine(); //skip header
			String line;
			while( (line=r.readLine()) != null ) {
				add(line);				
			}
			r.close();
			
			//sort and index
			ArrayList<Gene> ref = new ArrayList<Gene>();
			Iterator<Gene> it = hash.values().iterator();
			while( it.hasNext() ) {
				Gene g = it.next();
				ref.add(g);
			}
			Collections.sort(ref);
			
			String oldContig=null;
			int num = 0;
			for( int i = 0; i < ref.size(); i++ ) {
				Gene g = ref.get(i);
				if( !g.chr.equals(oldContig) ) {
					num=0;
					oldContig = g.chr;
				}
				g.num=num++;
			}
		}

		public String check(Gene g, String target) throws IOException {
			Iterator<String> rgs = g.rg.iterator();
			StringBuffer sb = new StringBuffer();
			while( rgs.hasNext() ) {
				String current = rgs.next();
				if( current.startsWith(prefix) ) {
					Gene ref = hash.get(current);
					if( sb.length()>0 ) sb.append("; ");
					sb.append(ref.chr+","+ref.middle+","+ref.strand+","+ref.num+","+ref.name);
				}
			}
			return sb.toString();			
		}

		public void add( String line ) {
			String[] split = line.split("\t");
			split[0]=prefix+split[0];
			Gene g = hash.get(split[0]);
			if( g==null ) {
				g = new Gene( split[0] );
				hash.put( split[0], g );
			}
			g.extend( split[4], Integer.parseInt(split[5]), Integer.parseInt(split[6]), Integer.parseInt(split[7]) );
		}
	}
	
	/**
	 * A class representing a gene.
	 * 
	 * @author Jens Keilwagen
	 */
	static class Gene implements Comparable<Gene>{
		String name;
		String chr;
		int strand, start, end, middle, num=-1;
		HashSet<String> rg;
		
		Gene( String name ) {
			this.name=name;
			chr=null;
			strand=0;
			middle=-1;
			rg=null;
		} 
		
		void extend( String chr, int strand, int start, int end ) {
			if( this.chr == null ) {
				this.chr=chr;
				this.strand=strand;
				this.start=start;
				this.end=end;
			} else {
				if( this.chr.equals(chr) && this.strand == strand ) {
					this.start = Math.min( start, this.start );
					this.end = Math.max( end, this.end );
				} else {
					throw new RuntimeException("Check chromsome/contig and strand for gene " + name);
				}
			}
			middle = (end+start)/2;
		}
		
		public int compareTo( Gene g ) {
			int d = chr.compareTo(g.chr);
			if( d == 0 ) {
				d = Integer.compare(middle,  g.middle );
			}
			return d;
		}
	}

	@Override
	public ToolParameterSet getToolParameters() {
		try{
			return
				new ToolParameterSet( getShortName(),
					new SimpleParameter(DataType.STRING,"tag","the tag used to read the GeMoMa annotations",true,GeMoMa.TAG),
					new ParameterSetContainer( "references", "", new ExpandableParameterSet( new SimpleParameterSet(		
							new SimpleParameter(DataType.STRING,"prefix","the prefix can be used to distinguish predictions from different input files (=reference organisms)", false, new RegExpValidator("\\w*")),
							new FileParameter( "assignment", "the assignment file of this reference organism, which combines parts of the CDS to transcripts", "tabular", true, new FileExistsValidator() )
					), "reference", "", 1 ) ),
					new FileParameter( "gene annotation file", "GFF file containing the gene annotations predicted by GAF", "gff,gff3", true, new FileExistsValidator(), true )
				);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}

	}

	@Override
	public String getToolName() {
		return "Synteny checker";
	}

	@Override
	public String getShortName() {
		return "SyntenyChecker";
	}

	@Override
	public String getDescription() {
		return "extracts the reference genes from the annotation";
	}

	@Override
	public String getHelpText() {
		return 
				"This tool can be used to determine syntenic regions between target organism and reference organism based on similiarity of genes. "
				+ "The tool returns a table of reference genes per predicted gene. This table can be easily visualized with an R script that is included in the GeMoMa package."
				+ MORE;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "tabular", "reference gene table")
		};
	}
	
	@Override
	public ToolResult[] getTestCases( String path ) {
		// TODO missing test cases
		return null;
	}
}