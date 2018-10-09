package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import projects.gemoma.GeMoMa.IntArrayComparator;

/**
 * This class finalizes an annotation. It predicts UTRs for a given set of CDS annotations using RNA-seq evidence and allows to generate generic gene and transcript names.
 * 
 * @author Jens Keilwagen
 */
public class AnnotationFinalizer implements JstacsTool {
	
	static boolean rename;
	
	/*public static void main( String[] args ) throws Exception {
		CLI cli = new CLI(new AnnotationFinalizer());
		cli.run(args);
	}*/
	
	static interface Feature {
		Feature addFeature( String[] split );
	}
	
	static class Gene implements Feature, Comparable<Gene> {
		ArrayList<Transcript> list;
		String[] split;
		int min=-100, max=-100;
		
		Gene( String[] split ) {
			this.split = split;
			list = new ArrayList<Transcript>();
		}
		
		public int strand() {
			return split[6].charAt(0)=='+' ? 1 : -1;
		}
		
		public Transcript addFeature( String[] split ){
			Transcript t = new Transcript(split);
			list.add(t);
			return t;
		}

		void determine() {
			min = Integer.MAX_VALUE;
			for( int i = 0; i < list.size(); i++ ) {
				min = Math.min( list.get(i).determine(), min);
			}
			max = Integer.MIN_VALUE;
			for( int i = 0; i < list.size(); i++ ) {
				max = Math.max( Integer.parseInt(list.get(i).split[4]), max);
			}			
		}
		
		void set() {
			determine();
			split[3] = "" + min;
			split[4] = "" + max;
			Collections.sort(list);
		}
		
		@Override
		public int compareTo(Gene g) {
			if( min == -100 ) {
				determine();
				//Collections.sort(list);
			}
			if( g.min == -100 ) {
				g.determine();
				//Collections.sort(g.list);
			}
			return Integer.compare( min, g.min );
		}
		
		public String toString() {
			StringBuffer all = new StringBuffer();
			all.append(toString(split));
			for( int i = 0; i < list.size(); i++ ) {
				all.append(list.get(i).toString());
			}
			return all.toString();
		}
		
		static String toString( String[] split ) {
			String res = "";
			for( int i = 0; i < split.length; i++ ) {
				res += split[i] + (i+1==split.length?"\n":"\t");
			}
			return res;
		}

		public void genericName(String name) {
			split[8]="Name=" + name + ";" + split[8];
			for( int i = 0; i < list.size(); i++ ) {
				Transcript t = list.get(i);
				t.split[8]="Name=" + name + "." + (i+1)+";" + t.split[8];
			}			
		}
	}
	
	static class Transcript implements Feature, Comparable<Transcript> {
		ArrayList<CDS> list;
		ArrayList<UTR> upUTR;
		ArrayList<UTR> downUTR;
		String[] split;
		String id;
		double score;
		
		Transcript( String[] split ) {
			this.split = split;
			String[] attr = split[8].split(";");
			int l = 0;
			while( l < attr.length && !attr[l].startsWith("ID=") ) {
				l++;
			}
			id=attr[l].substring(3);
			
			if( rename ) {
				//GeMoMa specific!
				l = 0;
				while( l < attr.length && !attr[l].startsWith("score=") ) {
					l++;
				}
				score = Double.parseDouble(attr[l].substring(6));
			} else {
				score = 0;
			}
			
			list = new ArrayList<CDS>();
			upUTR = new ArrayList<UTR>();
			downUTR = new ArrayList<UTR>();
		}
		
		public int strand() {
			return split[6].charAt(0)=='+' ? 1 : -1;
		}

		public int determine() {
			int min = getMin();
			int max = getMax();
			split[3] = "" + min;
			split[4] = "" + max;
			return min;
		}
		
		public int getMin() {
			int min = Integer.MAX_VALUE;
			for( int i = 0; i < upUTR.size(); i++ ) {
				min = Math.min(min, upUTR.get(i).getMin() );
			}
			for( int i = 0; i < list.size(); i++ ) {
				min = Math.min(min, list.get(i).getMin() );
			}
			for( int i = 0; i < downUTR.size(); i++ ) {
				min = Math.min(min, downUTR.get(i).getMin() );
			}
			return min;
		}
		
		public int getMax() {
			int max = Integer.MIN_VALUE;
			for( int i = 0; i < upUTR.size(); i++ ) {
				max = Math.max(max, upUTR.get(i).getMax() );
			}
			for( int i = 0; i < list.size(); i++ ) {
				max = Math.max(max, list.get(i).getMax() );
			}
			for( int i = 0; i < downUTR.size(); i++ ) {
				max = Math.max(max, downUTR.get(i).getMax() );
			}
			return max;
		}

		@Override
		public CDS addFeature(String[] split) {
			CDS c = new CDS(split);
			list.add(c);
			return c;
		}
		
		public void addUTR( boolean up, String[] split ) {
			UTR utr = new UTR(split);
			if( up ) {
				upUTR.add( utr );
			} else {
				downUTR.add( utr );
			}
		}
		
		public String toString() {
			StringBuffer all = new StringBuffer();
			all.append(Gene.toString(split));
			for( int i = 0; i < upUTR.size(); i++ ) {
				all.append(upUTR.get(i).toString());
			}
			for( int i = 0; i < list.size(); i++ ) {
				all.append(list.get(i).toString());
			}
			for( int i = 0; i < downUTR.size(); i++ ) {
				all.append(downUTR.get(i).toString());
			}
			return all.toString();
		}

		@Override
		public int compareTo(Transcript t) {
			return Double.compare(score, t.score);
		}
	}
	
	static class CDS implements Feature {
		String[] split;
		
		CDS( String[] split ) {
			this.split = split;
		}

		@Override
		public Feature addFeature(String[] split) {
			throw new RuntimeException("not supported");
		}
		
		public int getMin() {
			return Integer.parseInt(split[3]);
		}
		
		public int getMax() {
			return Integer.parseInt(split[4]);
		}
				
		public String toString() {
			return Gene.toString(split);
		}
	}
	
	static class UTR implements Feature {
		String[] split;
		
		UTR( String[] split ) {
			this.split=split;
		}
		
		@Override
		public Feature addFeature(String[] split) {
			throw new RuntimeException("not supported");
		}
		
		public int getMin() {
			return Integer.parseInt(split[3]);
		}
		
		public int getMax() {
			return Integer.parseInt(split[4]);
		}
				
		public String toString() {
			return Gene.toString(split);
		}		
	}
	
	static StringBuffer begin = new StringBuffer();
	
	static HashMap<String, HashMap<String,Gene>> read( String fName, Protocol protocol, String tag ) throws IOException {
		String[] tags = {"gene", tag, "CDS"};//GeMoMa specific
		//read annotation
		BufferedReader r = new BufferedReader( new FileReader( fName ) );
		ArrayList<String[]>[] list = new ArrayList[tags.length];
		HashMap<String,Integer> tagIndex = new HashMap<String, Integer>();
		for( int i = 0; i < list.length; i++ ) {
			list[i] = new ArrayList<String[]>();
			tagIndex.put( tags[i],  i );
		}
		
		String line;
		boolean gff=true;//TODO
		String[] split;
		boolean start = true;
		while( (line=r.readLine()) != null ) {
			if( gff && line.equalsIgnoreCase("##FASTA") ) {//http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
				//protocol.append("Stop reading the annotation file because of '##FASTA'\n"); 
				break;  
			}
			if( line.length() == 0 ) continue;
			if( line.startsWith("#") ) {
				if( start ) {
					begin.append(line+"\n");
				}
				continue;
			}
			start=false;
			if( !gff ) {
				int index = line.indexOf('#');
				if( index > 0 ) {
					line = line.substring(0, index);
				}
			}
			
			split = line.split("\t");
			Integer i = tagIndex.get( split[2] );
			if( i!= null ) {
				list[i].add(split);
			}
		}		
		r.close();
		
		//parse annotation
		HashMap<String,HashMap<String,Gene>> res = new HashMap<String, HashMap<String,Gene>>();
		HashMap<String, Feature> all = new HashMap<String,Feature>();
		int old=0;
		for( int i = 0; i < list.length; i++ ) {
			for( int j =0; j < list[i].size(); j++ ) {
				split = list[i].get(j);
				
				if( i == 0 ) {
					String chr = split[0];
					HashMap<String, Gene> h = res.get(chr);
					if( h == null ) {
						h = new HashMap<String, Gene>();
						res.put(chr, h);
					}
					String[] attr = split[8].split(";");
					int k = 0;
					while( k < attr.length && !attr[k].startsWith("ID=") ) {
						k++;
					}
					if( k < attr.length ) {
						String id = attr[k].substring(3);
						Gene g = new Gene(split);
						h.put(id, g );
						all.put(id, g );
					} else {
						System.out.println("WARNING1: " + Arrays.toString(split) );
					}
				} else {
					String[] attr = split[8].split(";");
					int k = 0;
					while( k < attr.length && !attr[k].startsWith("Parent=") ) {
						k++;
					}
					int l = 0;
					while( l < attr.length && !attr[l].startsWith("ID=") ) {
						l++;
					}
					if( k < attr.length && l < attr.length ) {
						String parent = attr[k].substring(7);
						Feature f = all.get(parent);
						Feature n = f.addFeature(split);
						String id = attr[l].substring(3);
						
						//System.out.println(id + "\t" + parent );
						
						all.put( id, n );
					} else {
						System.out.println("WARNING2: " + Arrays.toString(split) );
					}
				}
			}
			protocol.append( tags[i] + "s: " + (all.size()-old) + "\n" );
			old=all.size();
		}
		return res;
	}
	
	static class SequenceIDComparator implements Comparator<String> {
	    public int compare(String o1, String o2) {
	        int diff = extractInt(o1) - extractInt(o2);
	        if( diff == 0 ) {
	        	diff = o1.compareTo(o2);
	        }
	        return diff;
	    }

	    int extractInt(String s) {
	        String num = s.replaceAll("\\D", "");
	        // return 0 if no digits found
	        return num.isEmpty() ? 0 : Integer.parseInt(num);
	    }
	}
	
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {
		progress.setIndeterminate();
		
		if( GeMoMa.seqs == null ) {
			GeMoMa.fill(protocol, false, 0, parameters.getParameterForName("genome").getValue().toString(), null, 
					(Integer) parameters.getParameterForName("reads").getValue(), 
					(ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(3)).getValue(),
					(ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(5)).getValue() );
		}
				
		SimpleParameterSet renamePS = (SimpleParameterSet) parameters.getParameterForName("rename").getValue();
		int n = renamePS.getNumberOfParameters();
		int digits=-100;
		String prefix, infix=null, suffix = "";
		String[] del = null;
		rename = n>0;
		if( rename ) {
			prefix = renamePS.getParameterForName("prefix").getValue().toString();
			digits = (Integer) renamePS.getParameterForName("digits").getValue();
			if( n>2 ) {
				infix = renamePS.getParameterForName("infix").getValue().toString();
				suffix = renamePS.getParameterForName("suffix").getValue().toString();
				del = renamePS.getParameterForName("delete infix").getValue().toString().split(",");
				for( int i = 0; i < del.length; i++ ) {
					del[i] = del[i].trim();
				}
			}
		} else {
			prefix=null;
		}
		
		protocol.append("\n");
		
		//annotation
		String tag = parameters.getParameterForName("tag").getValue().toString();
		HashMap<String, HashMap<String,Gene>> annotation = 
				//Extractor.read( parameters.getParameterForName("annotation").getValue().toString(), null, protocol);
				read( parameters.getParameterForName("annotation").getValue().toString(), protocol, tag );

		
		//iterate over all chromosomes/contigs
		File gffFile = GeMoMa.createTempFile(getShortName()+"-");
		BufferedWriter w = new BufferedWriter(new FileWriter(gffFile));
		w.append(begin);
		w.append(GeMoMa.INFO + getShortName() + " " + getToolVersion() + "; ");
		String info = JstacsTool.getSimpleParameterInfo(parameters);
		if( info != null ) {
			w.append("SIMPLE PARAMETERS: " + info );
		}
		w.newLine();
		
		int num = 1;
		String[] keys = new String[annotation.size()];
		annotation.keySet().toArray(keys);
		Arrays.sort( keys, new SequenceIDComparator() );
		for( String c : keys ) {
			String newC = c;
			String pref = null;
			if( rename ) {
				if( n==2 ) {
					pref = prefix;
				} else {
					for( int i = 0; i < del.length; i++ ) {
						newC = newC.replaceAll(del[i], "");
					}
					pref = prefix + newC + infix;
					num=1;
				}
			}
			
			int[][] fwdCov = GeMoMa.coverage[0].get(c);
			int[][] revCov = GeMoMa.coverage[1].get(c);
			
			Collection<Gene> gg = annotation.get(c).values();
			ArrayList<Gene> genes = new ArrayList<>(gg);
			Collections.sort(genes);
			
			int[][][] acc = GeMoMa.acceptorSites.get(c);
			int[][][] don = GeMoMa.donorSites.get(c);
			
			//System.out.println(c + "\t" + newC + "\t" + fwdCov + "\t" + revCov + "\t" + acc + "\t" + don );
			
			if( fwdCov!= null || revCov != null ) {
				for( Gene g : genes ) {					
					//get coverage for current strand
					int strand = g.strand();
					int[][] a = acc==null?null:acc[strand==1?0:1];
					int[][] d = don==null?null:don[strand==1?0:1];
					int[][] cov = strand == 1 ? fwdCov : revCov;
					if( cov != null ) {
						for( Transcript t : g.list ) {
							boolean extend;
							extend = extendUTR( t, -1, cov, a, 1 ); //upstream
							extend |= extendUTR( t, 1, cov, d, 0 ); //downstream
							
							//if( extend ) System.exit(0);
						}
					}
					g.set();
				}
			}
			
			Collections.sort(genes);
			for( Gene g : genes ) {
				if( rename ) {
					g.genericName( pref + String.format("%0"+digits+"d", num++) + suffix );
				}
				w.append(g.toString());
			}
		}
		w.close();
		
		ArrayList<TextResult> res = new ArrayList<TextResult>();
		res.add( new TextResult(defResult, "Result", new FileParameter.FileRepresentation(gffFile.getAbsolutePath()), "gff3", getToolName(), null, true) );
				
		return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());

	}
	
	private static boolean extendUTR( Transcript t, int d, int[][] cov, int[][] splice, int sIdx ) { //XXX include introns?
		int strand = t.strand();
		int min = t.getMin();
		int max = t.getMax();
		int dir = d*strand;
		
		String utr=(d==1?"three_prime_UTR":"five_prime_UTR");
		
		LinkedList<String[]> utrs = new LinkedList<String[]>();
		
		int start = dir==-1 ? min : max;
		boolean extend = true;
		while( extend ) {
			int firstBase = start+dir;
			extend=false;
			
			//find stretch that is has some positive coverage
			int h=start;
			if( cov != null ) {
				int idx = Arrays.binarySearch(cov, new int[]{firstBase}, IntArrayComparator.comparator[2] );
				if( idx < 0 ) {
					idx = -(idx+1);
					idx = Math.max(0, idx-1);//has to be reduced since we did not find a match
				}
				if( cov[idx][0] <= firstBase && firstBase <= cov[idx][1] ) {
					//System.out.println(dir);
					do {
						h = cov[idx][dir==-1?0:1];
						idx+=dir;
					} while( idx >= 0 && idx < cov.length && cov[idx][dir==-1?1:0] == h+dir );
				}
			}
			int lastBase = h+dir;
			if( dir==-1 ) {
				int xxx = lastBase;
				lastBase = firstBase;
				firstBase = xxx;
			}
			
			//find splice site in interval [firstBase, lastBase]
			if( splice != null && splice[sIdx]!= null) {
				int idx = Arrays.binarySearch(splice[sIdx], firstBase );
				if( idx < 0 ) {
					idx = -(idx+1);
					//we dont need this here: idx = Math.max(0, idx-1);//has to be reduced since we did not find a match
				}
				int maxSplitReads = -1, maxIdx=-1;
				while( idx < splice[sIdx].length && splice[sIdx][idx] <= lastBase ) {
					if( splice[2][idx] > maxSplitReads ) {
						maxSplitReads = splice[2][idx];
						maxIdx = idx;
					}
					idx++;
				}
				if( maxIdx> 0 ) {
					idx=Arrays.binarySearch(cov, new int[]{splice[sIdx][maxIdx]+dir}, IntArrayComparator.comparator[2] );
					if( idx < 0 ) {
						idx = -(idx+1);
						idx = Math.max(0, idx-1);//has to be reduced since we did not find a match
					}
					if( cov[idx][2] <= splice[2][maxIdx] ) {
						extend=true;
						//TODO
						if( dir == -1 ) {
							firstBase = splice[sIdx][maxIdx]-dir;
						} else {
							lastBase = splice[sIdx][maxIdx]-dir;
						}
						start = splice[1-sIdx][maxIdx]; 
					}
				}
			}
			if( !extend ) {
				if( dir == -1 ) {
					firstBase -= dir;
				} else {
					lastBase -= dir;
				}
			}
			
			//add utr;
			if( firstBase <= lastBase ) {
				String[] split = {
					t.split[0],
					toolName,
					utr,
					""+ firstBase,
					""+ lastBase,
					".",
					t.split[6],
					".",
					null
				};
				
				if( d == -1 ) {
					utrs.addFirst( split );
				} else {
					utrs.addLast( split );
				}
			}
		}
		
		int u = 1;
		while( utrs.size() > 0 ) {
			String[] split;
			split = utrs.removeFirst();
			split[8] = /*"ID=" + t.id + "." + utr + "."+ u++ +";"+*/ "Parent="+t.id;
			t.addUTR(d==-1, split);
		}
		
		return utrs.size()>0;
	}
	
	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getShortName(),
					new FileParameter( "genome", "The genome file (FASTA), i.e., the target sequences in the blast run. Should be in IUPAC code", "fasta", true ),
					new FileParameter( "annotation", "The predicted genome annotation file (GFF)", "gff", true ),
					
					new SimpleParameter( DataType.STRING, "tag", "A user-specified tag for transcript predictions in the third column of the returned gff. It might be beneficial to set this to a specific value for some genome browsers.", true, GeMoMa.TAG ),					
					
					new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
							new FileParameter( "introns file", "Introns (GFF), which might be obtained from RNA-seq", "gff", false )
						), "introns", "", 1 ) ),
					new SimpleParameter( DataType.INT, "reads", "if introns are given by a GFF, only use those which have at least this number of supporting split reads", true, new NumberValidator<Integer>(1, Integer.MAX_VALUE), 1 ),

					new ParameterSetContainer( "coverage", "", new ExpandableParameterSet( new SimpleParameterSet(	
						new SelectionParameter( DataType.PARAMETERSET, 
								new String[]{"UNSTRANDED", "STRANDED"},
								new Object[]{
									//unstranded coverage
									new SimpleParameterSet(
											new FileParameter( "coverage_unstranded", "The coverage file contains the unstranded coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true )
									),
									//stranded coverage
									new SimpleParameterSet(
											new FileParameter( "coverage_forward", "The coverage file contains the forward coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true ),
											new FileParameter( "coverage_reverse", "The coverage file contains the reverse coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true )
									)
								},  "coverage file", "experimental coverage (RNA-seq)", true
						)
					), "coverage", "", 1 ) ),
					
					new SelectionParameter( DataType.PARAMETERSET, 
							new String[]{"COMPOSED", "SIMPLE", "NO"},
							new Object[]{
								//generic renaming
								new SimpleParameterSet(
										new SimpleParameter( DataType.STRING, "prefix", "the prefix of the generic name", true ),
										new SimpleParameter( DataType.STRING, "infix", "the infix of the generic name", true, "G" ),
										new SimpleParameter( DataType.STRING, "suffix", "the suffix of the generic name", true, "0" ),
										new SimpleParameter( DataType.INT, "digits", "the number of informative digits", true, new NumberValidator<Integer>(4, 10), 5 ),
										new SimpleParameter( DataType.STRING, "delete infix", "a comma-separated list of infixes that is deleted from the sequence names before building the gene/transcript name", true, "" )
								),
								//simple renaming
								new SimpleParameterSet(
										new SimpleParameter( DataType.STRING, "prefix", "the prefix of the generic name", true ),
										new SimpleParameter( DataType.INT, "digits", "the number of informative digits", true, new NumberValidator<Integer>(4, 10), 5 )
								),
								//no renaming
								new SimpleParameterSet()
							},
							new String[] {
								"composed renaming: <prefix><contig/chr*><infix><number><suffix>; recommended for whole genome (re-)annotations",
								"simple renaming: <prefix><number>; recommended for gene family annotations",
								"no renaming"
							},
							"rename", "allows to generate gernic gene and transcripts names (cf. attribute &quot;Name&quot;)", true
					)
			);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	static String toolName = "AnnotationFinalizer";
	
	public String getToolName() {
		return toolName;
	}
	
	public String getToolVersion() {
		return GeMoMa.VERSION;
	}
	
	public String getShortName() {
		return toolName;
	}

	public String getDescription() {
		return "predicts UTRs for a given set of CDS annotations";
	}

	public String getHelpText() {
		return 
			"**What it does**\n\nThis tool finalizes an annotation. It predicts for UTRs for annotated coding sequences and allows to generate generic gene and transcript names. Stranded RNA-seq libraries are recommended for UTR prediction to avoid too long predictions for genes in close proximity. Please use *ERE* to preprocess the mapped reads.\n\n"
			+ GeMoMa.REF;
	}

	private static final String defResult = "final_annotation";
	
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff3", defResult),
		};
	}
}