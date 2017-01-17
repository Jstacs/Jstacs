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
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import javax.naming.OperationNotSupportedException;

import de.jstacs.DataType;
import de.jstacs.algorithms.alignment.Alignment;
import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.PairwiseStringAlignment;
import de.jstacs.algorithms.alignment.cost.AffineCosts;
import de.jstacs.algorithms.alignment.cost.Costs;
import de.jstacs.algorithms.alignment.cost.MatrixCosts;
import de.jstacs.algorithms.graphs.UnionFind;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.tools.ui.galaxy.Galaxy;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Time;
import projects.gemoma.Tools.Ambiguity;

/**
 * Parsing tblastn hits to gene models using a dynamic programming algorithm with extensions for splice sites, start and stop codon.
 * 
 * @author Jens Keilwagen
 */
public class GeMoMa implements JstacsTool {

	private static DecimalFormat decFormat = new DecimalFormat("###.####",DecimalFormatSymbols.getInstance(Locale.US));

	/**
	 * The index of the score in the blast output.
	 */
	private static final int scoreIndex = 13;
		
	//user?
	private double contigThreshold;//threshold for initial solutions (filtering of contigs)
	private double regionThreshold;//threshold for regions
	private double hitThreshold;//threshold for hits
	
	private final int MIN_INTRON_LENGTH = 30;//the minimal intron length used in the DP to determine possible connections between hits of the same query exon
	private int MAX_INTRON_LENGTH;//the maximal intron length used in the DP to determine possible connections
	private int INTRON_GAIN_LOSS;
	private final int MAX_GAP = 5;
	private final int missingAA = 10;//the minimal number of missing AA used for a self-loop in the DP 
	private final int ignoreAAForSpliceSite = 30;// the maximal number of AA in a (blast) hit to be ignore
	
	private static boolean approx;
	private static boolean avoidStop; // avoid STOP in codons in gene models?
	private static double eValue;// the e-value for filtering tblastn hits
	private static int predictions;//how many predictions should be made

	//in
	private static HashMap<String, String> cds, protein;
	private static HashMap<String, String> seqs;
	private static HashMap<String,Character> code;
	private static HashMap<String,HashMap<String,int[][]>> transcriptInfo;
	private static HashMap<String, String> selected = null; //optional

	//out
	private BufferedWriter gff, predicted, blastLike, genomic;
	private String prefix, tag;
	
	//alignment
	private DiscreteAlphabet aaAlphabet;
	private AlphabetContainer alph;
	private MyAlignment align, align1, align2;
	private double[][] matrix;
	private int gapOpening;
	private int gapExtension;
	private Costs cost;
	private Tools.Ambiguity ambiguity = Ambiguity.AMBIGUOUS;

	//splicing
	private static HashMap<String, int[][][]> donorSites;
	private static HashMap<String, int[][][]> acceptorSites;
	private static boolean sp;
	//coverage
	private static HashMap<String, int[][]>[] coverage;	
	//canonical 
	private static final String[] DONOR = {"GT", "GC"};
	private static final String ACCEPTOR = "AG";
	private static int[] intronic = { DONOR[0].length(), ACCEPTOR.length()};
	
	//statistics
	private int numberOfLines,
		numberOfTestedStrands,
		numberOfPairwiseAlignments;
	private int[] stats = new int[5];

	private boolean writeGenomic;
	
	//for logging
	private Time time;
	private boolean verbose;
	private Protocol protocol;
	
	private int maxSize;
	private long timeOut, maxTimeOut;
	
	public GeMoMa( int maxSize, long timeOut, long maxTimeOut ) {
		this.maxSize = maxSize;
		this.timeOut = timeOut;
		this.maxTimeOut = maxTimeOut;
	}
	
	/**
	 * The main method for running the tool.
	 * 
	 * @param args the parameters for the tool.
	 * 
	 * @throws Exception forwarded from {@link CLI#run(String[])}
	 */
	public static void main(String[] args) throws Exception{
		File jarfile = new File(Galaxy.class.getProtectionDomain().getCodeSource().getLocation().toURI().getPath());
		if( jarfile != null ) {
			System.out.println("jar time stamp: " + new Date(jarfile.lastModified()) + "\n" );
		}
		
		int maxSize = -1;
		long timeOut=3600, maxTimeOut=60*60*24*7;
		
		if( args.length == 0 ) {
			System.out.println( "If you start with the tool with \"CLI\" as first parameter you can use the command line interface, otherwise you can use the Galaxy interface.");
		} else {
			if( args[0].equalsIgnoreCase("CLI") || args[0].equalsIgnoreCase("wiki") ) {
				CLI cli = new CLI( "CLI", null, new Extractor(maxSize), new ExtractRNAseqEvidence(), new GeMoMa(maxSize, timeOut, maxTimeOut), new GeMoMaAnnotationFilter() );
				if( args[0].equalsIgnoreCase("CLI") ) {
					String[] part = new String[args.length-1];
					System.arraycopy(args, 1, part, 0, part.length);
					cli.run(part);
				} else {
					System.out.println("Creating tables for the Jstacs wiki as separate files.");
					cli.wiki();
				}
			} else {
				if( "--create".equals(args[0]) ) {
					if( args.length != 1 ) {
						try {
							maxSize = Integer.parseInt(args[1]);
							timeOut = Long.parseLong(args[2]);
							maxTimeOut = Long.parseLong(args[3]);
							System.out.println("Try to parse <maxSize> <timeOut> <maxTimeOut> for the Galaxy integration");
						} catch (Exception e){
							System.out.println("Forwarding the arguments: " +  Arrays.toString(args) );
						}
					}
					System.out.println(maxSize + "\t" + timeOut + "\t" + maxTimeOut );
				}
				Galaxy galaxy = new Galaxy("", false, new Extractor(maxSize), new ExtractRNAseqEvidence(), new GeMoMa(maxSize, timeOut, maxTimeOut), new GeMoMaAnnotationFilter() );
				galaxy.run(args);
			}
		}
	}
	
	private InputStream getInputStream( Parameter parameter, String alternative ) throws FileNotFoundException {
		InputStream in;
		if( parameter.isSet() ) {
			in = new FileInputStream( parameter.getValue().toString() );
		} else {
			in = Extractor.class.getClassLoader().getResourceAsStream( alternative );
		}
		return in;
	}
	
	private static File createTempFile( String prefix ) throws IOException {
		File f = File.createTempFile(prefix, "_GeMoMa.temp", new File("."));//default temp directory?
		f.deleteOnExit();
		return f;
	}
	
	private static void readCoverage( int index, Parameter p, Protocol protocol, boolean verbose ) throws IOException {
		//if( p!= null && p.isSet() ) {
			protocol.append("read coverage file: " + p.getName() + "\n");
			coverage[index] = new HashMap<String, int[][]>();
			BufferedReader r = new BufferedReader( new FileReader( p.getValue().toString() ) );
			int threshold = 1;
			String chr = null;
			ArrayList<int[]> list = new ArrayList<int[]>();
			String line;
			int i = 0;
			while( (line=r.readLine()) != null ) {
				if( i==0 && line.startsWith("track type=bedgraph") ) {
					//ignore
				} else{
					String[] split = line.split( "\t" );
					if( chr == null || !split[0].equals(chr)) {
						if( chr != null && list.size()>0 ) {
							//Collections.sort(list, CovComparator.def);
							coverage[index].put( chr, list.toArray( new int[0][] ) );
							list.clear();
						}
						chr = split[0];
					}
					int reads;
					try {
						reads = Integer.parseInt(split[3]);
					} catch( Exception e ) {
						if( verbose ) {
							protocol.appendWarning("Could not parse number of reads. Set coverage to " +threshold + ". Line: " + line + "\n" );
						}
						reads = threshold;
					}
					if( reads >= threshold ) {
						list.add( new int[]{Integer.parseInt(split[1])+1, Integer.parseInt(split[2]), reads} );
					}
				}
				i++;
			}
			r.close();
			if( chr != null && list.size()>0 ) {
				//Collections.sort(list, CovComparator.def);
				coverage[index].put( chr, list.toArray( new int[0][] ) );
			}
		//}
	}
	
	@Override
	public ToolResult run( ParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {
		this.protocol=protocol;
		
		BufferedReader r = null;
		String line, old=null;
		
		MAX_INTRON_LENGTH = (Integer) parameters.getParameterForName( "maximum intron length" ).getValue();
		eValue = (Double) parameters.getParameterForName( "e-value" ).getValue();
		contigThreshold = (Double) parameters.getParameterForName( "contig threshold" ).getValue();
		regionThreshold = (Double) parameters.getParameterForName( "region threshold" ).getValue();
		hitThreshold = (Double) parameters.getParameterForName( "hit threshold" ).getValue();
		
		predictions = (Integer) parameters.getParameterForName( "predictions" ).getValue();
		gapOpening = (Integer) parameters.getParameterForName("gap opening").getValue();
		gapExtension = (Integer) parameters.getParameterForName("gap extension").getValue();
		INTRON_GAIN_LOSS = (Integer) parameters.getParameterForName("intron-loss-gain-penalty").getValue();
		avoidStop = (Boolean) parameters.getParameterForName("avoid stop").getValue();
		approx = (Boolean) parameters.getParameterForName("approx").getValue();		
		verbose = (Boolean) parameters.getParameterForName("verbose").getValue();
		prefix = (String) parameters.getParameterForName("prefix").getValue();
		tag = (String) parameters.getParameterForName("tag").getValue();
		timeOut = (Long) parameters.getParameterForName( "timeout" ).getValue();
		
		writeGenomic = (Boolean) parameters.getParameterForName("genomic").getValue();
		
		Parameter p = parameters.getParameterForName("selected"); 
		if( p.isSet() ) {
			selected = Tools.getSelection( p.getValue().toString(), maxSize, protocol );
			protocol.append("selected: " + selected.size() + (selected.size()<100? "\t"+ selected : "")+"\n");
		}
		
		//outputs
		File gffFile = createTempFile("gff");
		File predictedFile = createTempFile("protein");
		File alignFile = createTempFile("align");
		File genomicFile = createTempFile("genomic");
		
		gff = new BufferedWriter( new FileWriter( gffFile ) );
		gff.append("##gff-version 3");
		gff.newLine();
		predicted = new BufferedWriter( new FileWriter( predictedFile ) );
		blastLike = new BufferedWriter( new FileWriter( alignFile ) );
		genomic = new BufferedWriter( new FileWriter( genomicFile ) );
				
		p = parameters.getParameterForName("query proteins");
		
		TranscriptPredictor tp = new TranscriptPredictor(
				progress,
				(String) parameters.getParameterForName("assignment").getValue(),
				(String) parameters.getParameterForName("cds parts").getValue(),
				(String) parameters.getParameterForName("target genome").getValue(),
				
				getInputStream(parameters.getParameterForName("genetic code"), "projects/gemoma/test_data/genetic_code.txt" ),
				getInputStream(parameters.getParameterForName("substitution matrix"), "projects/gemoma/test_data/BLOSUM62.txt" ),

				(String) (p.isSet() ? p.getValue() : null)
		);
		if( selected != null ) {
			progress.setLast(selected.size());
		}
		
		p = parameters.getParameterForName("introns"); 
		if( p!= null && p.isSet() ) {
			int threshold = (Integer) parameters.getParameterForName("reads").getValue();
			r = new BufferedReader( new FileReader( p.getValue().toString() ) );
			HashMap<String, ArrayList<int[]>[]> spliceHash = new HashMap<String, ArrayList<int[]>[]>();
			ArrayList<int[]>[] h;
			
			String[] donor = null;
			String acceptor = null;
			
			while( (line=r.readLine()) != null ) {
				if( line.charAt(0) != '#' ) {
					String[] split = line.split( "\t" );
					int reads;
					try {
						reads = Integer.parseInt(split[5]);
					} catch( NumberFormatException nfe ) {
						if( verbose ) {
							protocol.appendWarning("Could not parse number of reads. Set " +threshold + ": " + line  + "\n");
						}
						reads = threshold;
					}
					if( reads >= threshold ) {
						h = spliceHash.get(split[0]);
						if( h == null ) {
							h = new ArrayList[]{ new ArrayList<int[]>(), new ArrayList<int[]>() };
							spliceHash.put( split[0], h );
						}
						/*
						int idx= split[6].charAt(0)=='+'?0:1;
						
						//donor = first element, acceptor = second
						h[idx].add( new int[]{Integer.parseInt(split[3+idx]), Integer.parseInt(split[4-idx])} );
						*/
						char c = split[6].charAt(0);
						int a=Integer.parseInt(split[3]);
						int b =Integer.parseInt(split[4]);
						
						if( c=='.' ) {
							//try to identify on which strand the intron is located
							if( donor == null ) {
								donor = new String[]{Tools.rc(DONOR[0]),Tools.rc(DONOR[1])};
								acceptor = Tools.rc(ACCEPTOR);
							}
							String s = seqs.get(split[0]);
							String x = s.substring(a-1,a+1);
							String y = s.substring(b-3,b-1);
							boolean fwd =  (x.equals(DONOR[0]) || x.equals(DONOR[1])) && y.equals(ACCEPTOR);
							boolean bwd =  (y.equals(donor[0]) || y.equals(donor[1])) && x.equals(acceptor);
							if( fwd && !bwd ) {
								c='+';
							} else if( bwd && !fwd ) {
								c='-';
							}
							//System.out.println(fwd + "\t" + bwd + "\t"+x + " .. " + y);
						}
						
						if( c =='+' || c=='.' ) {
							h[0].add( new int[]{a, b, reads} );
						}
						if( c =='-' || c=='.' ) {
							h[1].add( new int[]{b, a, reads} );
						}
					}
				}
			}
			r.close();
			
			//reformat
			donorSites = new HashMap<String, int[][][]>();
			acceptorSites = new HashMap<String, int[][][]>();
			
			Entry<String, ArrayList<int[]>[]> e;
			int[][][] vals, help;
			int[] site;
			Iterator<Entry<String, ArrayList<int[]>[]>> it = spliceHash.entrySet().iterator();
			int num = 0;
			while( it.hasNext() ) {
				e = it.next();
				h = e.getValue();
				num += h[0].size() + h[1].size();
				vals = new int[2][][];
				help = new int[2][][];
				for( int k = 0; k < vals.length; k++ ) {
					help[k] = new int[h[k].size()][3];
					for( int m = 0; m<h[k].size(); m++ ) {
						site = h[k].get(m);
						help[k][m][0] = site[0];
						help[k][m][1] = site[1];
						help[k][m][2] = site[2];
					}
					
					Arrays.sort(help[k], IntArrayComparator.comparator[1]);
					vals[k] = new int[3][h[k].size()];
					for( int m = 0; m<h[k].size(); m++ ) {
						vals[k][0][m] = help[k][m][0];
						vals[k][1][m] = help[k][m][1];
						vals[k][2][m] = help[k][m][2];
					}
				}
				acceptorSites.put(e.getKey(), vals);
				
				vals = new int[2][][];
				for( int k = 0; k < vals.length; k++ ) {
					Arrays.sort(help[k], IntArrayComparator.comparator[0]);
					vals[k] = new int[3][h[k].size()];
					for( int m = 0; m<h[k].size(); m++ ) {
						vals[k][0][m] = help[k][m][0];
						vals[k][1][m] = help[k][m][1];
						vals[k][2][m] = help[k][m][2];
					}
				}
				donorSites.put(e.getKey(), vals);
			}
			protocol.append("possible introns from RNA-seq (split reads>="+threshold+"): " + num + "\n");
			sp = (Boolean) parameters.getParameterForName("splice").getValue();
		} else {
			acceptorSites = donorSites = null;
		}
		
		coverage = new HashMap[2];
		SimpleParameterSet sps = (SimpleParameterSet) parameters.getParameterForName("coverage").getValue();
		if( sps.getNumberOfParameters() == 2 ) {
			readCoverage(0, sps.getParameterForName("coverage_forward"), protocol, verbose);
			readCoverage(1, sps.getParameterForName("coverage_reverse"), protocol, verbose);
		} else if( sps.getNumberOfParameters() == 1 ) {
			readCoverage(0, sps.getParameterForName("coverage_unstranded"), protocol, verbose);
			coverage[1] = coverage[0];
		} else {
			coverage = null;
		}
		
		//read blast output and compute result
		boolean okay = true;
		ArrayList<TextResult> res = new ArrayList<TextResult>();
		Exception e = null;
		try {
			//collect blast hits per transcript, split for chromosome, strand and cds part
			HashMap<String,HashMap<Integer,ArrayList<Hit>>[]> hash  = new HashMap<String, HashMap<Integer,ArrayList<Hit>>[]>();
			r = new BufferedReader( new FileReader( (String) parameters.getParameterForName("tblastn results").getValue() ) );
			time = Time.getTimeInstance(null);
			protocol.append(tp.getHeading()+"\n");
			while( (line=r.readLine()) != null ) {
				if( line.length() > 0 ) {
					if( old == null || !line.startsWith(old) ) {
						if( old != null ) {
							//compute
							tp.compute(old, hash);
							//clear
							hash.clear();
						}
						old = line.substring(0, line.indexOf('\t'));
						if( transcriptInfo != null ) {
							old = old.substring(0, old.lastIndexOf('_')+1);
						}
					}
					//parse blast hit
					addHit(hash, line);
				}
			}
			tp.compute(old, hash);		
		} catch ( Throwable er ) {
			protocol.appendThrowable(er);
			if( er instanceof Exception ) {
				e = (Exception) er;
			} else {
				e = new Exception(er.getMessage());
				e.setStackTrace( er.getStackTrace() );
			}
			okay=false;
		} finally {
			
			//show statistics
			protocol.append("\n" + Arrays.toString(stats) + "\n" );
					
			//close output;
			gff.close();
			predicted.close();
			blastLike.close();
			genomic.close();
			
			if( okay ) {
				
				res.add( new TextResult("predicted annotation", "Result", new FileParameter.FileRepresentation(gffFile.getAbsolutePath()), "gff", getToolName(), null, true) );
				res.add( new TextResult("predicted protein", "Result", new FileParameter.FileRepresentation(predictedFile.getAbsolutePath()), "fasta", getToolName(), null, true) );

				if( (Boolean) parameters.getParameterForName("align").getValue() ) {
					res.add( new TextResult("alignment", "Result", new FileParameter.FileRepresentation(alignFile.getAbsolutePath()), "tabular", getToolName(), null, true) );
				}
				if( writeGenomic ) {
					res.add( new TextResult("genomic", "Result", new FileParameter.FileRepresentation(genomicFile.getAbsolutePath()), "fasta", getToolName(), null, true) );
				}
			}
			
			if( r != null ) {
				r.close();
			}
			tp.close();
		}
		if( okay ) {
			return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
		} else {
			throw e;
		}
	}

	/**
	 * Class for sorting (BLAST) {@link Hit}s.
	 * 
	 * @author Jens Keilwagen
	 */
	private static class HitComparator implements Comparator<Hit> {
		/**
		 * The default instances
		 */
		static HitComparator[] comparator = {
				new HitComparator(true),
				new HitComparator(false)
		};
		
		/**
		 * Strand of the {@link Hit}s.
		 */
		private boolean forward;
		
		private HitComparator( boolean f ) {
			forward=f;
		}
		
		public int compare(Hit o1, Hit o2) {
			if( forward ) {
				return o1.targetStart - o2.targetStart;
			} else {
				return o2.targetEnd - o1.targetEnd;
			}
		}	
	}
	
	/**
	 * Class for sorting (BLAST) {@link Hit}s.
	 * 
	 * @author Jens Keilwagen
	 */
	private static class IntArrayComparator implements Comparator<int[]> {
		/**
		 * The default instances
		 */
		static IntArrayComparator[] comparator = {
				new IntArrayComparator(0),
				new IntArrayComparator(1)
		};
		
		private int idx;
		
		private IntArrayComparator( int idx ) {
			this.idx = idx;
		}
		
		public int compare( int[] o1, int[] o2) {
			return Integer.compare( o1[idx], o2[idx] );
		}	
	}
	
	/**
	 * Parses one line of the tab-separated tblastn results, filters by e-Value and creates if a {@link Hit} if the e-value is smaller than a threshold.
	 * 
	 * @param hash this {@link HashMap} contains the hits
	 * @param line a line of the blast output
	 */
	@SuppressWarnings("unchecked")
	void addHit( HashMap<String, HashMap<Integer,ArrayList<Hit>>[]> hash, String line ) {
		String[] split = line.split("\t");
		
		if( Double.parseDouble(split[10]) > eValue ) {//skip if e-Value of the hit is too large
			return;
		}
		
		//get list for insertion
		boolean forward = Integer.parseInt(split[8]) < Integer.parseInt(split[9]);
		
		HashMap<Integer,ArrayList<Hit>>[] current = hash.get(split[1]);
		if( current == null ) {
			current = new HashMap[] {
				new HashMap<Integer,ArrayList<Hit>>(),
				new HashMap<Integer,ArrayList<Hit>>()
			};
			hash.put(split[1],current);
		}
		int idx = forward ? 0 : 1;
		Integer part = transcriptInfo == null ? 0 : new Integer( split[0].substring(1+split[0].lastIndexOf("_")) );
		ArrayList<Hit> lines = current[idx].get(part);
		if( lines == null ) {
			lines = new ArrayList<Hit>();
			current[idx].put(part,lines);
		}
		
		//insert
		Hit h = new Hit( split[0], split[1], 
				Integer.parseInt(split[6]), Integer.parseInt(split[7]), Integer.parseInt(split[22]),  
				Integer.parseInt(split[8]), Integer.parseInt(split[9]), 
				Integer.parseInt(split[scoreIndex]), split[20], split[21], "tblastn;" );
		
		lines.add(h);
		//h.addSplitHits(lines);
		numberOfLines++;
	}
	
	
	
	/**
	 * Class representing tblastn hits. Hence, the query is an amino acid sequence and the target is a genomic DNA sequence.
	 * 
	 * @author Jens Keilwagen
	 * 
	 * @see HitComparator
	 */
	private class Hit implements Cloneable {		
		/**
		 * The strand of the hit, i.e., the strand of the genomic sequence which gives a good alignment
		 */
		boolean forward;
		/**
		 * The IDs for query and target sequence.
		 */
		String queryID, targetID;
		/**
		 * The start and end position for the query and the target sequence.
		 */
		int queryStart, queryEnd, targetStart, targetEnd;
		/**
		 * The length of the query, i.e. the length of the part of the coding sequence represented by one exon.
		 */
		int queryLength;
		/**
		 * The aligned query and target sequences.
		 */
		String queryAlign, targetAlign;
		/**
		 * The score of this hit.
		 */
		int score;
		/**
		 * The part (of the gene) of this hit.
		 */
		int part;
		/**
		 * The offsets for first/last exon until start/stop codon
		 */
		int firstOffset, lastOffset;
		int firstAddScore, lastAddScore;
		
		/**
		 * Field for addition information for the hit, e.g., predictor; splice sites, start, and stop codon extensions
		 */
		String info;
		
		
		//these variables reflect the neighborhood of the hit
		/**
		 * the potential up- and downstream amino acid sequence of the hit until the next stop codon occurs inframe.
		 */
		StringBuffer up, down;
		StringBuffer seqUp, seqDown;
		/**
		 * The maximal genomic position of start or end of this hit. The start and end of the contig/chromosome as well as the stop codons are used to derive these values.
		 */
		int maxUpstreamStart, maxDownstreamEnd;
		IntList[] accCand, accCandScore;
		IntList[][] donCand, donCandScore;
		
		boolean border, splice;
		boolean de, ae; //donor/acceptor evidence
		HashSet<Integer> preDonor, postAcceptor;
		char firstAA;
		
		int phase;
		
		/**
		 * The standard constructor. 
		 */
		private Hit( String queryID, String targetID, int queryStart, int queryEnd, int queryLength, int blastTargetStart, int blastTargetEnd, int score,
				String queryAlign, String targetAlign, String info ) {
			this.queryID = queryID;
			this.targetID = targetID;
			
			part = transcriptInfo == null ? 0 : Integer.parseInt(queryID.substring(queryID.lastIndexOf('_')+1));
			
			this.queryStart = queryStart;
			this.queryEnd = queryEnd;
			
			this.queryLength = queryLength;
			
			this.targetStart = blastTargetStart;
			this.targetEnd = blastTargetEnd;
			
			this.score = score;
			
			this.queryAlign = queryAlign;
			this.targetAlign = targetAlign.replaceAll("[BJZ]","X");//replace ambiguous amino acids to "X"
			
			this.info = info;
			
			forward = blastTargetStart < blastTargetEnd;
			
			if( !forward ) {
				//switch
				int h = targetEnd;
				targetEnd = targetStart;
				targetStart = h;
			}
			
			accCand=null;
			donCand=null;
			accCandScore=null;
			donCandScore=null;
			maxUpstreamStart = maxDownstreamEnd = Integer.MIN_VALUE;
			firstOffset=lastOffset=0;
			firstAddScore= lastAddScore=0;
			
			seqUp = new StringBuffer();
			seqDown = new StringBuffer();
			
			border = splice = false;
			de = ae = false;
			preDonor = postAcceptor = null;
		}
		
		public Hit clone() throws CloneNotSupportedException {
			Hit clone = (Hit) super.clone();
			clone.accCand = ArrayHandler.clone(accCand);
			clone.donCand = ArrayHandler.clone(donCand);
			clone.accCandScore = ArrayHandler.clone(accCandScore);
			clone.donCandScore = ArrayHandler.clone(donCandScore);
			clone.up = up==null? null : new StringBuffer( up );
			clone.down = down==null? null : new StringBuffer( down );
			clone.seqUp = seqUp==null? null : new StringBuffer( seqUp );
			clone.seqDown = seqDown==null? null : new StringBuffer( seqDown );
			clone.preDonor = preDonor==null? null : (HashSet<Integer>) preDonor.clone();
			clone.postAcceptor = postAcceptor==null? null : (HashSet<Integer>) postAcceptor.clone();
			return clone;
		}
		
		public String toString() {
			return queryID + "\t" + targetID + "\t" + (forward?"+":"-")
					+ "\t" + queryStart + "\t" + queryEnd +"\t"+queryLength + "\t" + targetStart + "\t" + targetEnd
					+ "\t" + score + "\t" + queryAlign + "\t" + targetAlign + "\t" + info;
		}
		
		public boolean anyEvidence( Hit downstream ) {
			boolean res = false;
			int phase=-1, i=-1, pos=-1;
			for( phase = 0; !res && phase < 3; phase++ ) {
				for( i = 0; !res && i < downstream.accCand[phase].length(); i++ ) {
					pos = ((forward?downstream.targetStart:(downstream.targetEnd+1))+(forward?-1:1)*downstream.accCand[phase].get(i));
					res |= postAcceptor.contains(pos);
				}
			}
			return res;
		}
	
		/**
		 * This method splits the current hit at stop codons (*) and adds the (refined) parts into the given list.
		 * 
		 * @param add the list of valid {@link Hit}s
		 * 
		 * @throws WrongAlphabetException if the is a problem in the score computation
		 */
		public void addSplitHits(ArrayList<Hit> add) throws WrongAlphabetException {
			int idx = targetAlign.indexOf('*');
			if( idx < 0 || queryAlign.charAt(idx)=='*' ) {
				add.add(this);
			} else {
				//split at STOP codon?
				int old=0, qStart=queryStart, tStart = forward? targetStart : targetEnd, factor = forward ? 1 : -1, f, a, b;
				do {
					a=b=0;
					if( idx < 0 ) {
						f = targetAlign.length();
					} else {
						f = idx;
						if( queryAlign.charAt(idx)=='*' ) {
							f++;
						} else {
							if( queryAlign.charAt(idx)=='-' ) {
								a=1;
								b=0;
							} else {
								a=b=1;
							}
						}
					}
						
					if( f-old > 0 ) {
						String t = targetAlign.substring(old, f);
						String q = queryAlign.substring(old, f);
						
						//shorten remaining alignment from both side if there are any gaps
						int off1= 0, oq1=0, off2 = t.length()-1, oq2=0, n = 0;
						while( off1 < t.length() && (t.charAt(off1) =='-' || q.charAt(off1)=='-') ) {
							if( q.charAt(off1)=='-') {
								oq1++;
							}
							off1++;
						}
						while( off2 >= 0 && (t.charAt(off2) =='-' || q.charAt(off2)=='-') ) {
							if( q.charAt(off2)=='-') {
								oq2++;
							}
							n++;
							off2--;
						}
						off2++;
						
						if( off1 > off2 ) {
							//weg damit
							int qg = getGapLength(q);
							int tg = getGapLength(t);
							qStart+=q.length()-qg+b;
							tStart+=factor*3*(t.length()-tg+a);
						} else {							
							q=q.substring(off1, off2);
							t=t.substring(off1, off2);
						
							int sc = getScore(q, t);
							int qg = getGapLength(q);
							int tg = getGapLength(t);
							
							int tS = tStart+factor*3*oq1;
							if( sc > 0 )  {
								Hit h =  new Hit(queryID, targetID, qStart+(off1-oq1), qStart+(off1-oq1)+q.length()-qg-1, queryLength, tS, tS +factor*(3*(t.length()-tg)-1), sc, q, t, info + "split by '*';" );
								add.add( h );
							}
							qStart+=(off1-oq1) + q.length()-qg + (n-oq2) + b;
							tStart+=factor*3*(oq1 + t.length()-tg + oq2 + a);
						}
						
						
					} else { //multiple STOP codons in a row
						qStart+=b;
						tStart+=factor*3*a;
					}
					old=idx+1;
					idx = targetAlign.indexOf('*',old); 
				} while( old > 0 );
			}
		}
		
		/**
		 * This method fills several variables mainly used for splice sites prediction.
		 * 
		 * @param chr the current chromosome
		 * 
		 * @throws CloneNotSupportedException forwarded from {@link ArrayHandler}
		 * @throws WrongAlphabetException forwarded from {@link Sequence}
		 * 
		 *  
		 * @see Hit#accCand
		 * @see Hit#accCandScore
		 * @see Hit#donCand
		 * @see Hit#donCandScore
		 * 
		 * @see Hit#up
		 * @see Hit#down
		 */
		public void prepareSpliceCandidates( String chr ) throws CloneNotSupportedException, WrongAlphabetException {
			if( !splice ) {
				if( accCand == null )  {
					accCand = ArrayHandler.createArrayOf(new IntList(), 3);
					accCandScore = ArrayHandler.createArrayOf(new IntList(), 3);
					donCand = new IntList[2][];
					donCandScore = new IntList[2][];
				}
				int l = Math.abs(targetStart-1-targetEnd), add, t = 3*ignoreAAForSpliceSite;
				int eDon = 0;
				int eAcc = 0;
				
				//add = the length which is used inside the exon to find a splice site
				if( l / 3 < t ) {
					add=l/3;
					add-=add%3; //;)
				} else {
					add=t;
				}

				String a = targetAlign.replaceAll( "-", "" );
				
				Sequence tPart, qPart;
				//acceptor
				up = new StringBuffer();
				Tools.getMaximalExtension(chr, forward, false, forward?-1:1, forward ? targetStart+add-1 : targetEnd-add, '*','*', seqUp, up, a.substring(0,add/3), eAcc, intronic[1]-1, MAX_INTRON_LENGTH, code );
				int s = seqUp.length()-add-eAcc;
				if( forward ) {
					maxUpstreamStart = targetStart - s+intronic[1];
				} else {
					maxUpstreamStart = targetEnd + s-intronic[1];
				}
				up=up.delete(up.length()-(add/3), up.length());
			
				int anz = 0, z, start=-10, end=-10, ref=-10, f=0, idx, d, m, old;
				int[][][] vals;
				//known acceptors
				if( acceptorSites != null ) {
					vals = acceptorSites.get(this.targetID);
					if( vals != null ) {
						preDonor = new HashSet<Integer>();
						
						z = forward?0:1;
						for( int i = 0; i < 3; i++) {
							accCand[i].clear();
						}

						if( forward ) {
							start = maxUpstreamStart;
							end = this.targetStart + add;
							ref = targetStart;
							f = 1;
						} else {
							start = this.targetEnd - add;
							end = maxUpstreamStart+1;
							ref = targetEnd+1;
							f = -1;
						}
	
						idx = Arrays.binarySearch(vals[z][1], start);
						if( idx < 0 ) {
							idx = -(idx+1);
						}
						old=-1;
						while( idx < vals[z][1].length && vals[z][1][idx] <= end ) {
							if( vals[z][1][idx] != old ) {
								d = f*(ref-vals[z][1][idx]);
								m = d%3;
								if( m < 0 ) {
									m+=3;
								}
								accCand[m].add( d );
								anz++;
								old=vals[z][1][idx];
							}
							preDonor.add(vals[z][0][idx]);
							idx++;
						}
					}
					ae = anz > 0;
				}

				//fall back
				if( acceptorSites == null || (sp && anz == 0) ) {
					getCandidates( seqUp, accCand, ACCEPTOR, true, true, s );
				}

				int max = Integer.MIN_VALUE;
				for( int p = 0; p < 3; p++ ) {
					for( int i = 0; i < accCand[p].length(); i++ ) {
						max = Math.max( max, accCand[p].get(i) );
					}
				}
				
				String query = cds.get(queryID);
				
				if( max > Integer.MIN_VALUE ) {
					qPart = Sequence.create( alph, query.substring(0,queryEnd) );
					if( max > 0 ) {
						max = (int) Math.floor(max/3d);
						tPart = Sequence.create( alph, up.substring(up.length() - max)+a );
					} else {
						max = -(int) Math.ceil(-max/3d);
						tPart = Sequence.create( alph, a.substring(-max) );
					}

					try {
						qPart = qPart.reverse();
						tPart = tPart.reverse();
					} catch( OperationNotSupportedException onse ) {
						//does not happen
						throw new RuntimeException( onse.getMessage() );
					}
					align.computeAlignment( AlignmentType.GLOBAL, qPart, tPart );
					for( int p = 0; p < 3; p++ ) {
						for( int i = 0; i < accCand[p].length(); i++ ) {
							int pos = accCand[p].get(i);
							if( pos > 0 ) {
								pos = (int)Math.floor(pos/3d);
							} else {
								pos = -(int)Math.ceil(-pos/3d);
							}
							pos = tPart.getLength() - (max - pos);
							accCandScore[p].add( (int) -align.getCost(qPart.getLength(), pos)-score );
						}
					}
				}
				
				//donor
				down = new StringBuffer();
				Tools.getMaximalExtension(chr, forward, true, forward?1:-1, forward ? targetEnd-add : (targetStart+add-1), '*','*', seqDown, down, a.substring(a.length() - add/3), eDon, intronic[0]-1, MAX_INTRON_LENGTH, code );
				s = add+eDon;
				if( forward ) {
					maxDownstreamEnd = targetEnd + seqDown.length()-s -intronic[0];
				} else {
					maxDownstreamEnd = targetStart - seqDown.length()+s +intronic[0];
				}
				down = down.delete(0,add/3);
				
				for( m = 0; m < donCand.length; m++ ) {
					donCand[m] = ArrayHandler.createArrayOf(new IntList(), 3);
					donCandScore[m] = ArrayHandler.createArrayOf(new IntList(), 3);
				}
				
				anz = 0;
				//known donors
				if( donorSites != null ) {
					vals = donorSites.get(this.targetID);
					if( vals != null ) {
						postAcceptor = new HashSet<Integer>();
						
						z = forward?0:1;
						for( int i = 0; i < 3; i++) {
							donCand[0][i].clear();
						}

						if( forward ) {
							start = this.targetEnd - add;
							end = maxDownstreamEnd+1;
							ref = targetEnd+1;
							f = -1;
						} else {
							start = maxDownstreamEnd;
							end = this.targetStart + add;
							ref = targetStart;
							f = 1;
						}
	
						idx = Arrays.binarySearch(vals[z][0], start);
						if( idx < 0 ) {
							idx = -(idx+1);
						}
						old = -1;
						while( idx < vals[z][0].length && vals[z][0][idx] <= end ) {
							if( vals[z][0][idx] != old ) {
								d = f*(ref-vals[z][0][idx]);
								m = d%3;
								if( m < 0 ) {
									m+=3;
								}
								donCand[0][m].add( d );
								anz++;
								old = vals[z][0][idx];
							}
							postAcceptor.add(vals[z][1][idx]);
							idx++;
						}
					}
					de = anz > 0;
				}
//for( int i = 0; i < 3; i++ ) { System.out.println(donCand[0][i]); donCand[0][i].clear(); } anz=0;
				//fall back
				if( donorSites == null || (sp && anz == 0) )
				{
					for( m = 0; m < DONOR.length; m++ ) {
						getCandidates( seqDown, donCand[m], DONOR[m], false, false, s );
					}
					
					//GC directly after the hit
					for( int i = 0; i < 3; i++ ) {
						if( add+i+2 <= seqDown.length() && seqDown.substring(add+i, add+i+2).equals(DONOR[1]) ) {
							donCand[0][i].add(i);
						}
					}/**/
				}					

//System.out.println(); for( int i = 0; i < 3; i++ ) { System.out.println(donCand[0][i]); }		
				max = Integer.MIN_VALUE;
				for( m = 0; m < donCand.length; m++ ) {
					for( int p = 0; p < 3; p++ ) {
						for( int i = 0; i < donCand[m][p].length(); i++ ) {
							max = Math.max( max, donCand[m][p].get(i) );
						}
					}
				}
				
				if( max > Integer.MIN_VALUE ) {					
					qPart = Sequence.create( alph, query.substring(queryStart-1) );
					if( max > 0 ) {
						tPart = Sequence.create( alph, a + down.substring(0, (int)Math.floor(max/3d)) );
					} else {
						tPart = Sequence.create( alph, a.substring(0,a.length() - (int) Math.ceil(-max/3d)) );
					}
					align.computeAlignment( AlignmentType.GLOBAL, qPart, tPart );
					
					for( m = 0; m < donCand.length; m++ ) {
						for( int p = 0; p < 3; p++ ) {
							for( int i = 0; i < donCand[m][p].length(); i++ ) {
								int pos = donCand[m][p].get(i);
								if( pos > 0 ) {
									pos = a.length() + (int)Math.floor(pos/3d);
								} else {
									pos = a.length() - (int) Math.ceil(-pos/3d);
								}
								donCandScore[m][p].add( (int) -align.getCost(qPart.getLength(), pos) -score );
							}
						}
					}
				}
				splice = true;
			}
		}
			
		private void getCandidates( StringBuffer region, IntList[] sites, String p, boolean add, boolean reverse, int ref ) {
			for( int i = 0; i < 3; i++) {
				sites[i].clear();
			}
			int idxSite=0;
			while( idxSite < region.length() && (idxSite=region.indexOf(p,idxSite))>= 0 ) {
				
				int pos = idxSite;
				if( add ) {
					 pos += p.length();
				}
				if( reverse ) {
					pos = ref-pos;
				} else {
					pos = pos - ref;
				}
				int i = pos%3;
				if( i < 0 ) {
					i = 3+i;
				}
				sites[i].add(pos);
				
				idxSite++;
			}
		}
		
		/**
		 * This method sets the splice site of a hit. Only one site is set by one method call.
		 * Also sets {@link Hit#info}.
		 * 
		 * @param donor donor splice site?
		 * @param add number of bp to be add
		 */
		public void setSpliceSite( boolean donor, int add ) {
//System.out.println(queryID + "\t" + donor + "\t" + add);
			if( donor ) {
				info+="donor";
				if( add != 0 ) {
					if( forward ) {
						info+=" (" + targetEnd+ ")";
						targetEnd+=add;
					} else {
						info+=" (" + targetStart+ ")";
						targetStart-=add;
					}
					if( add < 0 ) {
						int del = (int)Math.ceil(-add/3d), i  = targetAlign.length()-1, r = 0;
						while( i >= 0 && del > 0 ) {
							if( targetAlign.charAt(i) != '-' ) {
								del--;
							}
							if( queryAlign.charAt(i) != '-' ) {
								r++;
							}
							i--;
						}
						targetAlign = targetAlign.substring(0, i+1);
						queryAlign = queryAlign.substring(0, i+1);
						queryEnd -= r;
					}
				}
			} else {				
				info+="acceptor";
				if( add != 0 ) {
					phase = (add%3);
					if( phase < 0 ) {
						phase+=3;
					}
					
					if( forward ) {
						info+=" (" + targetStart + ")";
						targetStart-=add;
					} else {						
						info+=" (" + targetEnd + ")";
						targetEnd+=add;
					}
				}
				
				if( add < 0 ) {
					int del = (int)Math.ceil(-add/3d), i  = 0, r = 0;
					while( i < targetAlign.length() && del > 0 ) {
						if( targetAlign.charAt(i) != '-' ) {
							del--;
						}
						if( queryAlign.charAt(i) != '-' ) {
							r++;
						}
						i++;
					}
					targetAlign = targetAlign.substring(i);
					queryAlign = queryAlign.substring(i);
					queryStart += r;
				}
			}
			info+=";";
		}
		
		public void setBorderConstraints( boolean firstExon, boolean lastExon ) throws IllegalArgumentException, WrongAlphabetException {
			if( !border ) {
				if( firstExon ) {
					if( !(targetAlign.charAt(0) == 'M' && queryStart == 1) ) {
						firstOffset = up.length();
						if( firstOffset > 0 && up.indexOf("M")>=0 ) {//enlarge
							Sequence missing = Sequence.create(alph,cds.get(queryID).substring(0,queryStart));
							StringBuffer up = this.up.reverse();//TODO make this nice
							try {
								missing = missing.reverse();
							} catch( OperationNotSupportedException onse ) {
								//does not happen
								throw new RuntimeException( onse.getMessage() );
							}
							Sequence upSeq = Sequence.create(alph, up.substring(0,up.lastIndexOf("M")+1) );
							align.computeAlignment( AlignmentType.GLOBAL, missing, upSeq );
							
							int best, idx;
							if( targetAlign.charAt(0) == 'M' ) {
								best = -(gapOpening + missing.getLength()*gapExtension);
								idx=-1;
							} else {
								best = Integer.MIN_VALUE;
								idx= Integer.MIN_VALUE;
							}
							
							int j=-1, current;
							while( (j=up.indexOf("M", j+1)) >= 0 ) {
								current = (int) -align.getCost( missing.getLength(), j+1 );
								if( current > best ) {
									best = current;
									idx=j;
								}
							}
							if( idx == Integer.MIN_VALUE ) {
								protocol.append("HELP\n");
							}
							firstOffset = idx+1;
						} else {//shorten
							firstOffset = targetAlign.indexOf('M');
							if( firstOffset > 0 ) {
								firstOffset = firstOffset - getGapLength(targetAlign.substring(0, firstOffset));
								boolean okay = true;
								if( forward ) {
									if( targetStart + 3*firstOffset >= targetEnd ) {
										okay=false;
									}
								} else {
									if( targetEnd - 3*firstOffset <= targetStart ) {
										okay=false;
									}
								}
								if( okay ) {
									firstOffset *= -1;
								}
							} else {
								firstOffset = 0;
							}
						}
					}
					
					//XXX if( firstOffset != 0 ) 
					{
						Sequence q = Sequence.create(alph,cds.get(queryID).substring(0,queryEnd));
						Sequence t = Sequence.create(alph, Tools.translate(0, getDNA(3*firstOffset,0), code, false, ambiguity) );
						
						align.computeAlignment(AlignmentType.GLOBAL, q, t);
						int current = (int) -align.getCost( q.getLength(), t.getLength() );
						firstAddScore = current - score;
						firstAA=t.toString(0,1).charAt(0);
					}
				}
				if( lastExon ) {
					lastOffset = 0;
					if( targetAlign.charAt(targetAlign.length()-1) != '*' && down.length() > 0 && down.charAt(down.length()-1) == '*' ) {//enlarge
							lastOffset = down.length();
						Sequence q = Sequence.create(alph,cds.get(queryID).substring(queryStart-1));
						Sequence t = Sequence.create(alph, Tools.translate(0, getDNA(0,3*lastOffset), code, false, ambiguity) );
						
						PairwiseStringAlignment sa = align.getAlignment(AlignmentType.GLOBAL, q, t);
						//TODO numberOfPairwiseAlignments?
						int current = (int) -sa.getCost();
						lastAddScore = current - score;
					}
				}
				border=true;
			}
		}
		
		/**
		 * This method tries to extend a hit until start and stop codon if necessary.
		 * 
		 * @param firstExon if <code>true</code> extend until start codon
		 * @param lastExon if <code>true</code> extend until stop codon
		 */
		public void extend( boolean firstExon, boolean lastExon ) {
			if( firstExon ) {
				boolean problem = false;
				if( firstOffset!= 0 ) {
					if( forward ) {
						if( targetStart -3*firstOffset < targetEnd ) {
							targetStart -= 3*firstOffset;
						} else {
							problem=true;
						}
					} else {
						if( targetEnd + 3*firstOffset > targetStart ) {
							targetEnd += 3*firstOffset;
						} else {
							problem=true;
						}
					}
					if( !problem ) {
						info+="START (" + firstOffset + " aa);";
					} else {
						info+="START-problem";
					}
				}
			}
			if( lastExon ) {
				boolean problem = false;
				if( lastOffset > 0 ) {
					if( forward ) {
						if( targetStart < targetEnd + 3*lastOffset ) {
							targetEnd += 3*lastOffset;
						} else {
							problem=true;
						}
					} else {
						if( targetStart - 3*lastOffset < targetEnd ) {
							targetStart -= 3*lastOffset;
						} else {
							problem=true;
						}
					}
					if( !problem ) {
						info+="STOP (" + lastOffset + " aa);";
					} else {
						info+="STOP-problem";
					}
				}
			}
		}

		
		/**
		 * Returns the DNA of this {@link Hit}.
		 * 
		 * @param off1 the number of additional bp
		 * @param off2 the number of additional bp
		 * @return the DNA sequence corresponding to this {@link Hit}
		 */
		public String getDNA(int off1, int off2) {
			String chr = seqs.get(targetID);
			
			int offsetUp, offsetDown;
			if( forward ) {
				offsetUp = off1;
				offsetDown = off2;
			} else {
				offsetUp = off2;
				offsetDown = off1;
			}
			
			String r = chr.substring(targetStart-1-offsetUp, targetEnd+offsetDown);
			if( !forward ) {
				r = Tools.rc(r);
			}
			return r;
		}
		
		public int getLength() {
			return targetEnd-(targetStart-1);
		}
	}

	
	private static int getGapLength( String s ) {
		int g = 0;
		for( int i = 0; i < s.length(); i++ ) {
			if( s.charAt(i) == '-' ) {
				g++;
			}
		}
		return g;
	}
	
	//s1.length() == s2.length()
	private int getScore( String s1, String s2 ) throws WrongAlphabetException {
		int res = 0, gap=-1;
		for( int i = 0; i < s1.length(); i++ ) {
			char c1 = s1.charAt(i);
			char c2 = s2.charAt(i);
			if( c1 != '-' && c2 != '-' ) {
				res += matrix[aaAlphabet.getCode(""+c1)][aaAlphabet.getCode(""+c2)];
				gap=-1;
			} else {
				if( gap == -1 //new gap
						|| (gap == 1 && c2 == '-') //new gap in target
						|| (gap == 2 && c1 == '-') //new gap in query
				) {
					res += gapOpening;
					gap= c1 == '-' ? 1 : 2;
				}
				res+=gapExtension;
			}
			//protocol.appendln(c1 + "\t" + c2 + "\t" + res );
		}
		return -res;
	}
	
	private static int[][] getArray( String exonId, String phase ) {
		int[][] res = new int[2][]; 
		String[] split = exonId.split( "," );
		res[0] = new int[split.length];
		for( int i = 0; i < res[0].length; i++ ) {
			res[0][i] = Integer.parseInt(split[i].trim());
		}
		if( phase != null ) {
			split = phase.split( "," );
			res[1] = new int[split.length];
			for( int i = 0; i < res[1].length; i++ ) {
				res[1][i] = Integer.parseInt(split[i].trim());
			}
		}
		return res;
	}
	
	/**
	 * This class predicts gene models using a DP algorithm.
	 * 
	 * @author Jens Keilwagen
	 */
	private class TranscriptPredictor {	
		
		
		//variables for the DP
		/**
		 * Contains the sums of the scores of the stacked {@link Hit}s. Main DP field.
		 */
		private int[][] sums;
		private int[] parts, revParts, oldSize, phase;
		private int[] length, cumLength, revCumLength;
		private int[][] start, end, qStart, qEnd, qL;
		
		double[][] used;
		
		int[][][][] splice; //part 1 idx, anz, parts2 idx (relative), anz
		static final int NOT_COMPUTED = Integer.MAX_VALUE;
		static final int NO_SPLICE_VARIANT = Integer.MIN_VALUE;
				
		/**
		 * This filed is used to store the possible extensions for a splice variant
		 */
		private int[] refined = new int[2];
				
		/**
		 * Current gene name.
		 */
		private String geneName;
		
		/**
		 * Collects the potential solutions.
		 */
		private PriorityQueue<Solution> result;
		private Solution sol = new Solution();
		
		private ProgressUpdater progress;
		private ExecutorService executorService;
		
		private StringBuffer gffHelp;
						
		public TranscriptPredictor( ProgressUpdater progress, String assignment, String cdsFileName, String seqsFileName, InputStream codeIn, InputStream subsitutionMatrixIn, String proteinFileName ) throws Exception {
			BufferedReader r;
			String line ;
			String[] split;
			//read assignment
			this.progress = progress;
			if( assignment != null ) {
				transcriptInfo = new HashMap<String,HashMap<String,int[][]>>();
				HashMap<String,int[][]> c;
				r = new BufferedReader( new FileReader( assignment ) );
				while( (line=r.readLine())!= null ) {
					if( line.charAt(0) != '#' ) {
						split = line.split("\t");
						c = transcriptInfo.get( split[0] );
						if( c == null ) {
							c = new HashMap<String, int[][]>();
							transcriptInfo.put(split[0], c);
						}
						c.put(split[1].toUpperCase(), getArray(split[2], 
								null) );
								//TODO exonPhase split[3]) );
					}
				}
				r.close();
				this.progress.setLast(transcriptInfo.size());
				this.progress.setCurrent(0);
			} else {
				this.progress.setIndeterminate();
				transcriptInfo = null;
				parts = new int[1];
			}

			cds = Tools.getFasta(cdsFileName,15000,'\t');
			protocol.append("CDS: " + cds.size() + " / " + (transcriptInfo==null?cds.size():transcriptInfo.size()) + "\n" );// + "\t" + Arrays.toString(cdsInfo.keySet().toArray()) );
			
			//read genome
			seqs = Tools.getFasta(seqsFileName,20,' ');
			protocol.append("genome parts: " + seqs.size() + "\t" + Arrays.toString(seqs.keySet().toArray()) +"\n" );
						
			//read protein
			protein = Tools.getFasta(proteinFileName,15000,' ');		
			
			//read genetic code
			code = Tools.getCode(codeIn);
			code.put("NNN", 'X');//add for splitting at NNN (in align method)
			
			//read substitution matrix
			r = new BufferedReader( new InputStreamReader( subsitutionMatrixIn ) );
			while( (line=r.readLine()) != null && line.charAt(0)=='#' );
			String[] abc = line.split("\\s+");
			String[] abc2 = new String[abc.length-1];
			System.arraycopy( abc, 1, abc2, 0, abc2.length );
			aaAlphabet=new DiscreteAlphabet(true,abc2);
			alph = new AlphabetContainer(aaAlphabet);
			matrix = new double[abc2.length][abc2.length];
			int j = 0;
			while( (line=r.readLine()) != null ) {
				split = line.split("\\s+");
				for( int k = 1; k < split.length; k++ ) {
					matrix[j][k-1] = -Double.parseDouble(split[k]);
				}
				j++;
			}
			r.close();
			
			cost = new AffineCosts(gapOpening, new MatrixCosts(matrix, gapExtension));
				//new MatrixCosts(matrix, 5);
			align = new MyAlignment( cost, GeMoMa.this );
			align1 = new MyAlignment( cost, GeMoMa.this );
			align2 = new MyAlignment( cost, GeMoMa.this );
			
			result = new PriorityQueue<Solution>();
			
			executorService = Executors.newSingleThreadExecutor();
			
			gffHelp = new StringBuffer();
		}
		
		public void close() {
			executorService.shutdown();			
		}

		/**
		 * This method computes for all transcripts of one gene model the predicted gene models in the target organism 
		 * 
		 * @param name gene name ending with an underscore
		 * @param hash the collected tblastn {@link Hit}s
		 * 
		 * @throws Exception if something went wrong
		 */
		@SuppressWarnings("unchecked")
		public void compute( String name, HashMap<String,HashMap<Integer,ArrayList<Hit>>[]> hash ) throws Exception {
			geneName = name.substring(0,name.length()-1);
			
			HashMap<String,int[][]> transcript;
			int[][] current;
			String[] s;
			if( transcriptInfo != null ) {
				transcript = transcriptInfo.get(geneName);
				if( transcript == null ) {
					return;
				}
				s = new String[transcript.size()];
				transcript.keySet().toArray(s);
				Arrays.sort(s);
			} else {
				transcript = null;
				s = new String[]{ name.toUpperCase() };
			}
			HashMap<String, HashMap<Integer, ArrayList<Hit>>[]> data = new HashMap<String, HashMap<Integer,ArrayList<Hit>>[]>();
			HashMap<Integer,ArrayList<Hit>>[] v, w;
			ArrayList<Hit> hits;
			String[] h = new String[hash.size()];
			hash.keySet().toArray(h);
			
			//sort hits
			for( int j = 0; j < h.length; j++ ) {
				v = hash.get( h[j] );
				for( int f = 0; f < 2; f++ ) {
					sort(v[f], f==0);
				}
			}
			
			//compute for each transcript
			for( int i = 0; i < s.length; i++ ) {
				String info=null;
				if( selected==null || (info=selected.get(s[i]))!=null ) {
					time.reset();
					
					int max = 0;
					if( transcript != null ) {
						current = transcript.get(s[i]);
						parts = current[0];
						for( int p = 0; p < parts.length; p++ ) {
							if( parts[p] > max ) {
								max = parts[p];
							}
						}
						phase = current[1];
					}
					
					revParts = new int[max+1];
					if( transcript != null ) {
						Arrays.fill(revParts, -1);
						for( int p = 0; p < parts.length; p++ ) {
							revParts[parts[p]] = p;
						}
						
						length = new int[parts.length];
						for( int j = 0; j < parts.length; j++ ) {
							String st = geneName + (transcriptInfo==null ? "" : ("_" + parts[j]));
							st = cds.get( st );
							length[j] = st!=null ? st.length() : 0;
								
						}
					} else {
						revParts[0] = 0;
						length=new int[]{cds.get(name).length()};
					}
					
					cumLength = new int[length.length];
					for( int j = 1; j < length.length; j++ ) {
						cumLength[j] = cumLength[j-1] + length[j-1];
					}
					
					revCumLength = new int[length.length];
					for( int j = length.length-2; j >= 0; j-- ) {
						revCumLength[j] = revCumLength[j+1] + length[j+1];
					}
					
					int l = parts.length;
					sums = new int[l][];
					start = new int[l][];
					end = new int[l][];			
					qStart = new int[l][];
					qEnd = new int[l][];
					qL = new int[l][];
					
					if( transcriptInfo == null ) {
						computeTranscript( s[i], info, parts, hash );
					} else {
						boolean hasHits=false;
						data.clear();
						for( int j = 0; j < h.length; j++ ) {
							v = hash.get( h[j] );
							w = new HashMap[] {
									new HashMap<Integer,ArrayList<Hit>>(),
									new HashMap<Integer,ArrayList<Hit>>()
								};
							for( int f = 0; f < 2; f++ ) {
								for( int k = 0; k < parts.length; k++ ) {
									hits = v[f].get(parts[k]);
									if( hits != null ) {
										w[f].put( parts[k], clone( hits ) );
										hasHits=true;
									}
								}
							}
							data.put(h[j], w);
						}
						if( hasHits ) {
							computeTranscript( s[i], info, parts, data );
						}
					}
				}
			}
			numberOfLines=0;
		}
		
		public void sort( HashMap<Integer,ArrayList<Hit>> lines, boolean forward ) {
			Iterator<ArrayList<Hit>> it = lines.values().iterator();
			while( it.hasNext() ) {
				Collections.sort( it.next(), HitComparator.comparator[forward?0:1]  );
			}
		}
		
		public ArrayList<Hit> clone( ArrayList<Hit> original ) throws CloneNotSupportedException {
			ArrayList<Hit> clone = new ArrayList<Hit>( original.size() );
			for( int i = 0; i < original.size(); i++ ) {
				clone.add( original.get(i) );
			}
			return clone;
		}
		
		/**
		 * This method predicts gene model of a given transcript.
		 * 
		 * @param transcriptName the transcript name
		 * @param assign the assigned CDS parts
		 * @param hash the tblastn {@link Hit}s
		 * 
		 * @throws Exception if something went wrong
		 */
		private void computeTranscript( String transcriptName, String info, int[] assign, HashMap<String,HashMap<Integer,ArrayList<Hit>>[]> hash ) throws Exception {
			if( verbose ) {
				protocol.append("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
				protocol.append( geneName + "\t" + transcriptName + "\t" + Arrays.toString(parts) + "\t" + new Date() + (info==null?"":("\t"+info)) + "\n");
				protocol.append( Arrays.toString(length) + "\n");
				protocol.append( Arrays.toString(cumLength) + "\n");
				protocol.append( Arrays.toString(revCumLength) + "\n");
			}
			protocol.append( geneName+"\t"+transcriptName );
			protocol.flush();

			numberOfTestedStrands=numberOfPairwiseAlignments=0;
			oldSize = new int[parts.length];
			result.clear();
			
			HashMap<Integer,ArrayList<Hit>>[] lines;
			Collection<String> collection = hash.keySet();
			String chrom = null;
			int start = 0, end = -1, strand = -1;
			if( info != null && info.length()>0) {
				String[] pos = info.split("\t");
				if( hash.containsKey(pos[0]) ) {
					chrom = pos[0];
					collection = new ArrayList<String>();
					collection.add(chrom);
					if( pos.length > 1 ) {
						switch (pos[1]) {
							case "1": case "+": strand = 0; break;
							case "-1": case "-": strand = 1; break;
							default: strand=-1; break;
						}
					}
					if( pos.length > 2 ) {
						try {
							start = Integer.parseInt(pos[2]);
						} catch( Exception ex ) {
							start = 0;
						}
					}
					if( pos.length > 3 ) {
						try {
							end = Integer.parseInt(pos[3]);
						} catch( Exception ex ) {
							end=-1;
						}
					}

					lines = hash.get( chrom );
					if( lines != null ) {
						if( strand != -1 ) {
							if( lines[strand]!= null ) {
								filter( lines[strand], start, end );
							}
						} else {
							for( int i = 0; i < 2; i++ ) {
								if( lines[i]!= null ) {
									filter( lines[i], start, end );
								}
							}
						}
					}
					
					if( verbose ) protocol.append( "try to predict " + transcriptName + " on (" + chrom + " " + strand + ": " + start + " " + end + ") as user-specified\n");
				} else {
					if( verbose ) protocol.append( "\n" );
				}
			} else {
				if( verbose ) protocol.append( "\n" );
			}
			
			
			boolean out = true;
			final int STRAND = strand;
			final String CHROM = chrom;
			final HashMap<String,HashMap<Integer,ArrayList<Hit>>[]> HASH = hash;
			final Collection<String> COLLECTION = collection;
			int[] res = null;
			try {
				Future<int[]> f = executorService.submit(new Callable<int[]>(){
	                public int[] call() throws Exception {
						splice=null;
						used=null;
	                	HashMap<Integer,ArrayList<Hit>>[] lines;
	                	Iterator<String> it;
	        			int m = 0, k = 0, bestSumScore = Integer.MIN_VALUE, sumScore;
	        			if( STRAND == -1 ) {
	        				//forward DP: compute initial score for each chromosome/contig and strand
	        				it = COLLECTION.iterator();
	        				IntList score = new IntList();
	        				while( it.hasNext() ) {
	        					String c = it.next();
	        					lines = HASH.get( c );
	        	
	        					//check both strands
	        					for( int j = 0; j < lines.length; j++ ) {
	        						if( lines[j].size() > 0 ) {
	        							numberOfTestedStrands++;       							
	        							sumScore = forwardDP( j==0, lines[j], false );
	        							if( sumScore > bestSumScore ) {
	        								bestSumScore = sumScore;
	        							}
	        							score.add(sumScore);
	        							//verbose.writeln( c + "\t" + (j==0) + "\t" + sumScore + "\t" + new Date() );
	        						}
	        					}
	        				}
	        				
	        				//backward DP: compute final score and gene structure (b) on candidate chromosomes/contigs and strands
	        				it = COLLECTION.iterator();
	        				//Solution best = new Solution(), b = new Solution(), z = new Solution(), h;
	        				
	        				double threshold = bestSumScore*GeMoMa.this.contigThreshold;
	        				//verbose.writeln( threshold + "\t" + new Date() );
	        				while( it.hasNext() ) {
	        					String c = it.next();
	        					lines = HASH.get( c );
	        					
	        					//check both strands
	        					for( int j = 0; j < lines.length; j++ ) {
	        						if( lines[j].size() > 0 ) {
	        							if( score.get(m) >= threshold ) {
	        								k++;
	        								detailedAnalyse( c, j==0, lines[j] );
	        							}
	        							m++;
	        						}
	        					}
	        				}
	        			} else {
	        				lines = HASH.get( CHROM );
	        				if( lines != null && lines[STRAND].size()>0 ) {
	        					numberOfTestedStrands++;
	        					bestSumScore = forwardDP( STRAND==0, lines[STRAND], false );
	        					detailedAnalyse( CHROM, STRAND==0, lines[STRAND] );
	        					k=1;
	        				}
	        			}
	        			return new int[]{ bestSumScore, k };
	                }
	            });
				res = f.get(timeOut,TimeUnit.SECONDS);
	        } catch (TimeoutException e) {
	        	protocol.append( "\tInvocation did not return before timeout of " + timeOut + " seconds\n");
	        	out=false;
	        } finally {
	        	splice=null;
	        	used=null;
	        }
			
			if( out ) {
				String seq = protein != null ? protein.get(transcriptName) : null;
				int counts = 0;
				for( int i = 0; i < predictions && result.size()>0; i++ ) {
					Solution best = result.poll();
					if( best.hits.size() > 0 ) {
						counts++;
						best = best.refine(transcriptName);
									
						//write results
						blastLike.append( "#" + transcriptName + "\t" + best.score + "\t" + Arrays.toString(parts)+ (info==null?"":("\t"+info)) + "\n" );
						for( Hit t : best.hits ) {
							//blast like output
							if( verbose ) protocol.append(t.toString() + "\n");
							blastLike.append(t+"\n");
						}
						
						best.writeSummary( geneName, transcriptName, i);
						gff.append( ";score=" + best.score);
						
						gffHelp.delete(0, gffHelp.length());
						int anz = best.writeGFF( transcriptName, i, gffHelp );
						
						if( acceptorSites != null ) {
							gff.append( ";tae=" + (best.A==0? "?": decFormat.format(best.a/(double)best.A)) );
						}
						if( donorSites != null ) {
							gff.append( ";tde=" + (best.D==0? "?": decFormat.format(best.d/(double)best.D)) );
							gff.append( ";tie=" + (best.I==0? "?": decFormat.format(best.i/(double)best.I)) );
							gff.append( ";minSplitReads=" + (best.I==0? "?": best.minSplitReads) );
						}
						if( coverage != null && coverage[best.forward?0:1] != null ) {
							gff.append( ";tpc=" + decFormat.format(best.Cov/(double)best.Len) );
							gff.append( ";minCov=" + best.minC );
							gff.append( ";avgCov=" + decFormat.format(best.Sum/(double)best.Len) );
						}
						
						int id = 0, pos = 0, gap=-1, currentGap=0, maxGap=0;
						String s0, s1 = null, pred=null;
						PairwiseStringAlignment psa = null;
						if( seq != null ) {
							pred = best.getProtein();
							
							psa = align.getAlignment(AlignmentType.GLOBAL, Sequence.create(alph, seq), Sequence.create(alph, pred) );
							s0 = psa.getAlignedString(0);
							s1 = psa.getAlignedString(1);
							
							for( int p = 0; p < s1.length(); p++ ) {
								
								
								if( s0.charAt(p) == '-' ) {
									if( gap != 0 ) {
										if( currentGap > maxGap ) {
											maxGap = currentGap;
										}
										gap=0;
										currentGap=0;
									}
									currentGap++;
								} else if( s1.charAt(p) == '-' ) {
									if( gap != 1 ) {
										if( currentGap > maxGap ) {
											maxGap = currentGap;
										}
										gap=1;
										currentGap=0;
									}
									currentGap++;
								} else {
									//(mis)match
									if( gap != -1 ) {
										if( currentGap > maxGap ) {
											maxGap = currentGap;
										}
										gap=-1;
										currentGap=0;
									}
									if( s1.charAt(p) == s0.charAt(p) ) {
										id++;
										pos++;
									} else {
										if( matrix[aaAlphabet.getCode(s0.substring(p, p+1))][aaAlphabet.getCode(s1.substring(p, p+1))] > 0 ) {
											pos++;
										}
									}
								}
							}
							if( verbose ) {
								protocol.append(s0+"\n");
								protocol.append(s1+"\n");
							}
							
							//TODO additional GFF tags
							gff.append( ";iAA=" + decFormat.format(id/(double)s1.length()) );//+ ";maxGap=" + maxGap + ";alignF1=" + (2*aligned/(2*aligned+g1+g2)) ); 
						}
						
						//short info
						protocol.append( ((i == 0 && !verbose)? "" : (geneName+"\t"+transcriptName)) + "\t" + parts.length + "\t" + best.hits.size()
								+ "\t" + numberOfLines + "\t" + numberOfTestedStrands + "\t" + res[0]
								+ "\t" + res[1] + "\t" + (result.size()+1) + "\t" + numberOfPairwiseAlignments
								+ "\t" + best.score
								+ "\t" + best.hits.getFirst().targetID + "\t" + (best.forward?+1:-1)
								+ "\t" + (best.forward? best.hits.getFirst().targetStart + "\t" + best.hits.getLast().targetEnd : best.hits.getLast().targetStart + "\t" + best.hits.getFirst().targetEnd )
								+ "\t" + time.getElapsedTime() + "\t" + anz + "\t" + best.getInfo(gff) + "\t" + best.backup + "\t" + best.cut
								+ ( seq != null ?  
										"\t" + getGapCost(seq.length(), parts.length) /*-(seq.length()*gapExtension+gapOpening+parts.length*INTRON_GAIN_LOSS)*/
										+ "\t" + ((int)-psa.getCost()) + "\t" + getScore(seq, seq)
										+ "\t" + (pos/(double)s1.length()) + "\t" + (id/(double)s1.length()) + "\t" + maxGap + "\t" + seq.length() + "\t" + pred.length()
										: "" )
								+ (acceptorSites == null ? "" : ("\t" + (best.A==0? "?": decFormat.format(best.a/(double)best.A))) )
								+ (donorSites == null ? "" : ("\t" + (best.D==0? "?": decFormat.format(best.d/(double)best.D))) )
								+ (donorSites == null ? "" : ("\t" + (best.I==0? "?": decFormat.format(best.i/(double)best.I))) )
								+ (coverage == null ? "" : ("\t" + decFormat.format(best.Cov/(double)best.Len) + "\t" + best.minC))
								+ "\t" + best.similar(result.peek())
								+ "\n"
						);
						
						gff.newLine();
						
						gff.append(gffHelp);
						
						predicted.write( ">" + prefix+transcriptName + "_R" + i );
						predicted.newLine();
						predicted.write( best.getProtein() );
						predicted.newLine();
						
						genomic.flush();
						gff.flush();
						predicted.flush();
						blastLike.flush();
						protocol.flush();
					}
				}
				if( counts == 0 ) {
					protocol.append( "\n" );
					protocol.flush();
				}
			}
			progress.add(1);
		}
		
		private void filter( HashMap<Integer,ArrayList<Hit>> hash, int start, int end ) {
			if( start == 0 && end == -1 ) {
				return;
			}
			//filter out
			Hit h;
			if( verbose ) protocol.append("\n");
			for( int i = 0; i < parts.length; i++ ) {
				ArrayList<Hit> current = hash.get(parts[i]);
				if( verbose ) protocol.append(i + "\t" + parts[i] + "\t#=" + (current==null?0:current.size()) + "\n");
				if( current != null ) {
					int j = 0;
					while( j < current.size() ) {
						h = current.get(j);
						if( h.targetEnd < start || h.targetStart > end ) {
							if( verbose ) protocol.append("out\t"+ h+"\n");
							current.remove(j);
						} else {
							if( verbose ) protocol.append("in\t"+ h+"\n");
							j++;
						}
					}
				}
				
				if( verbose ) protocol.append(i + "\t" + parts[i] + "\t#=" + (current==null?0:current.size()) + "\n\n");
			}
		}
		
		public String getHeading() {
			return "gene\ttranscript\t#parts\t#predicted hits\t#blast hits\t#strands\tbest sum score"
					+ "\t#candidate strands\t#predictions\t#alignments"
					+ "\tbest final score"
					+ "\tchromosome\tstrand\tstart\tend\ttime\t#predicted parts\tfirst part\tfirst aa\tlast parts\tlast aa\t#*\tintron gain\tintron loss\tbackup\tcut"
					+ (protein != null?"\tminimal score\tcurrent score\toptimal score\t%positive\tpid\tmaxGap\tref_length\tpred_length":"")
					+ (acceptorSites == null ? "" :"\tacceptor evidence")
					+ (donorSites == null ? "" : "\tdonor evidence")
					+ (donorSites == null ? "" : "\tintron evidence")
					+ (coverage==null ? "" : "\tpercent covered\tmin coverage")
					+ "\tsimilar";
		}

		/**
		 * This method implements the dynamic programming (DP) algorithm for stacking the hits to gene models.
		 * It is implemented on one strand of one chromosome/contig.
		 * 
		 * If <code>{@link TranscriptPredictor#splice}!=null</code> consensus splice sites will be forced. 
		 * 
		 * @param forward the strand
		 * @param lines the {@link Hit}s
		 * 
		 * @return the best score
		 * 
		 * @throws WrongAlphabetException if there ambiguities
		 * 
		 * @see Hit#prepareSpliceCandidates(String)
		 * @see TranscriptPredictor#checkSpliceSites(String, Hit, Hit)
		 */
		private int forwardDP( boolean forward, HashMap<Integer,ArrayList<Hit>> lines, boolean gap ) throws WrongAlphabetException, IOException {
			if( lines.size() == 0 ) {
				for( int i = parts.length-1; i >= 0; i-- ) {
					sums[i] = start[i] = end[i] = qStart[i] = qEnd[i] = null;
				}
				bestIdx[0].clear();
				bestIdx[1].clear();
				return Integer.MIN_VALUE;
			}
			//forward DP to obtain the (initial) score for the region 
			Hit o;
			ArrayList<Hit> l;
			int b;
			String chr = null;
			for( int i = parts.length-1; i >= 0; i-- ) {
				l = lines.get(parts[i]);
				if( l == null || l.size() == 0 ) {
					sums[i] = start[i] = end[i] = qStart[i] = qEnd[i] = null;
				} else {				
					int anz = l.size();
					sums[i] = new int[anz];
					start[i] = new int[anz];
					end[i] = new int[anz];
					qStart[i] = new int[anz];
					qEnd[i] = new int[anz];
					qL[i] = new int[anz];
					if( splice != null && splice[i] == null ) {
//protocol.append( "create " + anz + "\n");
						splice[i] = new int[anz][parts.length-i][];
					}
					for( int j = anz-1; j >= 0; j-- ) {
						o = l.get(j);
						
						if( splice!= null && chr == null ) {
							chr = seqs.get(o.targetID);
						}
						
						start[i][j] = forward ? o.targetStart : o.targetEnd;
						end[i][j] = forward ? o.targetEnd : o.targetStart;
						qStart[i][j] = o.queryStart;
						qEnd[i][j] = o.queryEnd;
						qL[i][j] = o.queryEnd-o.queryStart+1;
						
						b = o.score;
						sums[i][j] = b + getCost( o, gap, false, i );
									
						//same & downstream exons
						int m=Math.min(i+MAX_GAP,parts.length), startIdx=j+1, max = 0;
						for( int k = i ; k < m; k++ ) {
							max = getMax(forward, lines, i, j, k, startIdx, chr, gap);
							if( sums[i][j] < b+max ) {
								sums[i][j]=b+max;
							}
							startIdx=0;
						}
					}
				}
			}
			return getBest( gap, lines );
		}
		
		int getCost( Hit o, boolean gap, boolean start, int i ) {
			if( !gap ) {
				return 0;
			} else {
				if( start ) {
					if( i == 0 ) {
						return o.firstAddScore;
					} else {
						return getGapCost(o.queryStart-1 + cumLength[i], i-1 );
					}
				} else {
					if( i == parts.length-1 ) {
						return o.lastAddScore;
					} else {
						return getGapCost(o.queryLength-o.queryEnd + revCumLength[i], parts.length-1-i);
					}
				}
			}
		}
		
		final int getGapCost( int cumLength, int exons ) {
			return -(cumLength > 0 ? gapOpening + gapExtension*cumLength : 0) - exons*INTRON_GAIN_LOSS;
		}
		
		final boolean check( boolean forward, int i , int j, int k, int m ) {
			return (forward?1:-1)*(start[k][m] - start[i][j]) > 0 //genomic start: first before second
					&& (forward?1:-1)*(end[k][m] - end[i][j]) > 0 //genomic end: first before second
					&& (k > i //different parts
							|| (qStart[i][j] < qStart[k][m] && qEnd[i][j] < qEnd[k][m] ) ); //same part but later aligned
		}
		
		final boolean check2( int i , int j, int k, int m ) {
			return i==k && (Math.max(qEnd[i][j] - qStart[k][m],0) / (double) Math.min(qL[i][j],qL[k][m])) > 0.5; 
		}		
		
		/**
		 * This method returns the score of the best matching hit for part <code>k</code> with the <code>j</code>-th {@link Hit} of part <code>i</code>.
		 * 
		 * @param forward the strand
		 * @param lines the {@link Hit}s
		 * @param i the current part
		 * @param j the index of the current {@link Hit} of part <code>i</code>
		 * @param k the part to be considered; <code>i&lt;=k</code>
		 * @param startIdx the start index of the {@link Hit} of part <code>k</code>
		 * @param chr the chromosome
		 * @return the best score
		 */
		private final int getMax( boolean forward, HashMap<Integer,ArrayList<Hit>> lines, int i, int j, int k, int startIdx, String chr, boolean gap) throws WrongAlphabetException {
			int max = getCost( lines.get(parts[i]).get(j), gap, false, i);
			
			ArrayList<Hit> ref = lines.get(parts[k]);
			if( ref == null ) {
				return max;
			}
			Hit upstreamHit = lines.get(parts[i]).get(j);
			
			int introns= Math.max(1,k-i);
//protocol.append(splice + "\t" + i + "\t" + j + "/" + ((splice!= null && splice[i]!=null)?splice[i].length:"?") + "\t" + (k-i) +"\n");//TODO
			if( splice != null && splice[i][j][k-i] == null ) {
				//protocol.appendln(o.targetID + "\t" + i + "\t" + j + "\t" + (k-i) );
				try {
					splice[i][j][k-i] = new int[ref.size()];
					Arrays.fill(splice[i][j][k-i], NOT_COMPUTED );
				} catch( OutOfMemoryError e ) {
					protocol.append(parts[i] + "\t" + length[i] + "\n" );
					protocol.append(splice[i].length+"\n");
					protocol.append(ref.size()+"\n");
					throw e;
				}
			}
			
			int len = 0;
			if( gap && sums[k] != null && i+1 != k ) {
				for( int v = i+1; v < k; v++ ) {
					len += length[v];
				}
			}
			
			for( int m = startIdx; sums[k]!= null && m < sums[k].length; m++ ) {
				int d = (forward?1:-1) * (start[k][m]-(end[i][j]+1)); //difference
				if( d < introns*MAX_INTRON_LENGTH ) {
					if( check(forward, i, j, k, m) ) {
						int adScore = 0;
						if( splice != null && splice[i][j][k-i] != null ) {
							if( splice[i][j][k-i][m] == NOT_COMPUTED ) {
								splice[i][j][k-i][m] = checkSpliceSites(chr, upstreamHit, ref.get(m), len);
							}
							adScore = splice[i][j][k-i][m];
						} else if( check2( i, j, k, m) ) {
							adScore=NO_SPLICE_VARIANT;
						}
						
						if( adScore!= NO_SPLICE_VARIANT//there is a splice variant
								&& max < sums[k][m]+adScore ) {//and it is promising
							max = sums[k][m]+adScore;
						}
					} else {
						if( splice != null && splice[i][j][k-i] != null ) {
							splice[i][j][k-i][m] = NO_SPLICE_VARIANT;
						}
					}
				} else {
					if( splice != null && splice[i][j][k-i] != null ) {
						Arrays.fill( splice[i][j][k-i], m, splice[i][j][k-i].length, NO_SPLICE_VARIANT );
					}
					break;
				}
			}
			
			return max;
		}
		
		/**
		 * This method performs a back tracking to obtain a solution (gene model) from the DP algorithm.
		 * 
		 * @param forward the strand
		 * @param lines the {@link Hit}s
		 * @param i the current part
		 * @param j the index of the current {@link Hit} of part <code>i</code>
		 * @param diff the allowed difference
		 * @param current the current {@link Solution}
		 * @param best the best {@link Solution}
		 */
		private void backTracking( boolean forward, HashMap<Integer,ArrayList<Hit>> lines, int i, int j, double diff, Solution current, Solution best, boolean gap ) {
			Hit o = lines.get(parts[i]).get(j);
			if( used != null ) {
				if( used[i] == null ) {
					used[i] = new double[lines.get(parts[i]).size()];
					Arrays.fill(used[i], Double.NEGATIVE_INFINITY);
				}
				if( diff >= 0 && diff >= used[i][j] ) {
					used[i][j]=diff;
				} else {
					return;
				}
			}

			if( current!= null ) {
				current.hits.add(o);
			}
			//protocol.append("check hit " + i+", " + j + "\t" + diff );
			diff = diff - (sums[i][j] - o.score);
			//protocol.appendln( "\t" + diff);
			
			//same & downstream exons
			int min=Math.min(i+MAX_GAP,parts.length), startIdx=j+1;
			for( int k = i; k < min; k++ ) {
				findNext( forward, lines, i, j, k, startIdx, diff, current, best, gap );
				startIdx=0;
			}
			
			if( current != null && sums[i][j] - (o.score + getCost(o, gap, false, i)) == 0 ) {
				if( best.hits.size() == 0 //no solution yet 
						|| best.compareTo(current)>0 //better than current solution
				) {
					best.set(current);
				}
			}
		}		
		
		/**
		 * This method adds the best matching between the <code>j</code>-th {@link Hit} of part <code>i</code> and one {@link Hit} of part<code>k</code> to the {@link Solution}. 
		 * 
		 * @param forward the strand
		 * @param lines the {@link Hit}s
		 * @param i the current part
		 * @param j the index of the current {@link Hit} of part <code>i</code>
		 * @param k the part to be considered; <code>i&lt;=k</code>
		 * @param startIdx the start index of the {@link Hit} of part <code>k</code>
		 * @param diff the allowed difference
		 * @param current the current {@link Solution}
		 * @param best the best {@link Solution}
		 * 
		 * @see TranscriptPredictor#backTracking(boolean, HashMap, int, int, double, Solution, Solution)
		 */
		private void findNext( boolean forward, HashMap<Integer,ArrayList<Hit>> lines, int i, int j, int k, int startIdx, double diff, Solution current, Solution best, boolean gap ) {
			if( sums[k] != null ) {
				int introns= Math.max(1,k-i);
				for( int m = startIdx; m < sums[k].length; m++ ) {			
					int d = (forward?1:-1) * (start[k][m]-(end[i][j]+1));
					if( d < introns*MAX_INTRON_LENGTH ) {
						if( check(forward, i, j, k, m) ){
							int adScore = 0;
							if( splice != null ) {
								adScore = splice[i][j][k-i][m];
							} else if( check2( i, j, k, m) ) {
								adScore=NO_SPLICE_VARIANT;
							}
							double newDiff = diff+adScore+sums[k][m];
							if( adScore != NO_SPLICE_VARIANT && newDiff>= 0 ) {
								//recursive
								backTracking( forward, lines, k, m, newDiff, current, best, gap );
								if( current!= null ) {
									current.hits.removeLast();
								}
							}
						}
					} else {
						break;
					}
				}
			}
		}

		/**
		 * This method does the refined analysis. It performs the following steps
		 * <ol>
		 * <li>DP algorithm</li>
		 * <li>reducing hits</li>
		 * <li>searching/adding missed parts</li>
		 * <li>DP algorithm with splicing</li>
		 * <li>back tracking</li>
		 * </ol>
		 * 
		 * @param chromosome the chromosome/contig name
		 * @param forward the strand
		 * @param lines the tblastn {@link Hit}s
		 */
		private void detailedAnalyse( String chromosome, boolean forward, HashMap<Integer,ArrayList<Hit>> lines ) throws CloneNotSupportedException, IllegalArgumentException, WrongAlphabetException, IOException {
			if( verbose ) {
				protocol.append( "all hits " + chromosome + " " + (forward?"+":"-") +"\n");
				show(lines);
			}
			//TASK 1: reduce blast hit set	
			int max = forwardDP(forward, lines, false);
			HashMap<Integer, ArrayList<Hit>> filtered = reduce(chromosome, forward, lines, avoidStop, false, max*hitThreshold);
			
			if( verbose ) {
				protocol.append( "filtered hits " + chromosome + " " + (forward?"+":"-")+"\n" );
				show(filtered);
			}
			
			//XXX split region?
			ArrayList<Hit> current, all = new ArrayList<Hit>();
			for( int j=0; j<parts.length; j++ ) {
				current = filtered.get(parts[j]);
				if( current != null && current.size() > 0 ) {
					all.addAll( current );
				}
			}
			Collections.sort( all, HitComparator.comparator[forward ? 0 : 1] );
			
			Hit now;
			int previous = -1, last=-1, prevPart = -1;
			boolean informative;
			IntList endIdx = new IntList(), startIdx = new IntList();
			for( int j=0; j<all.size(); j++ ) {
				now = all.get(j);
				informative = ((now.queryEnd-now.queryStart+1d) / Math.max(20d,now.queryLength)) >= 0.9;//XXX check this value
				if( informative ) 
				{
					if( previous >= 0 && prevPart>=revParts[now.part] ) {
						startIdx.add(last+1);
						last=previous;
						endIdx.add(j);
					}
					previous = j;
					prevPart = revParts[now.part];
				}
				if( verbose ) protocol.append( j + "\t" + informative + "\t" + now + "\n" );
			}
			startIdx.add(last+1);
			endIdx.add(all.size());
			
			boolean backup = false;
			if( endIdx.length() > 1 ) {
				int a = 0;
				HashMap<Integer,ArrayList<Hit>> region = new HashMap<Integer, ArrayList<Hit>>();
				for( int j=0; j<endIdx.length(); j++ ) {
					region.clear();
					region = new HashMap<Integer,ArrayList<Hit>>();
					int end = endIdx.get(j), k = startIdx.get(j);
					if( verbose ) protocol.append("check region: [" + k + ", " + end + ")\n");
					while( k < end ) {
						now = all.get(k);
						ArrayList<Hit> list = region.get( now.part );
						if( list == null ) {
							list = new ArrayList<Hit>();
							region.put( now.part, list );
						}
						list.add( now );
						k++;
					}
					if( forwardDP(forward, region, false) >= max*regionThreshold ) {
						result.add( detailedAnalyseRegion( chromosome, forward, region, backup ) );
						a++;
					} else {
						if( verbose ) protocol.append( "discard region\n" );
					}
				}
				if ( a > 0 ) {
					return;
				}
			}
			if( verbose ) protocol.append("check region: [" + startIdx.get(0) + ", " + endIdx.get(0) + ")\n");
			result.add( detailedAnalyseRegion( chromosome, forward, filtered, backup ) );
		}
		
		private void show( HashMap<Integer, ArrayList<Hit>> lines ) throws IOException {
			for( int h = 0; h <parts.length; h++ ) {
				ArrayList<Hit> list = lines.get(parts[h]);
				protocol.append(h + "\t" + parts[h]+"\t#=" + (list==null?0:list.size()) + "\n");
				for( int k = 0; list!= null && k < list.size(); k++ ) {
					protocol.append( list.get(k) +"\n" );
				}/**/
			}
			protocol.append("\n");
		}
		
		private HashMap<Integer, ArrayList<Hit>> reduce( String chromosome, boolean forward, HashMap<Integer,ArrayList<Hit>> lines, boolean avoidStop, boolean gap, double thresh ) throws IOException, WrongAlphabetException {
			used = new double[parts.length][];
			//mark unused blast hits
			ArrayList<Hit> list;
			for( int i = 0; i < sums.length; i++ ) {
				list = lines.get(parts[i]);
				for( int j = 0; sums[i] != null && j < sums[i].length; j++ ) {
					int c = getCost(list.get(j), gap, true, i);
					//System.out.println(i+"\t"+j + "\t" + (sums[i][j]+c - thresh));
					backTracking( forward, lines, i, j, sums[i][j]+c - thresh, null, null, gap );
				}
			}			
			
			if( splice != null ) {
				int[] anz = new int[parts.length];
				boolean changed = false;
				for( int i = 0; i < used.length; i++ ) {
					int a = 0;
					if( used[i] != null ) {
						for( int j = 0; j < used[i].length; j++ ) {
							if( used[i][j] >= 0 ) {
								a++;
							} else {
								changed = true;
							}
						}
					} else {
						list = lines.get(parts[i]);
						if( list != null && list.size()>0 ) {
							changed = true;
						}
					}
					anz[i] = a;
				}
				if( verbose ) protocol.append("#new hits: " + Arrays.toString(anz) + "\n");
				
				if( changed ) {
					int[][][][] splice2 = new int[parts.length][][][]; 
					int[] su, st, e, qS, qE, qL;
					for( int i = 0; i < used.length; i++ ) {
						if( anz[i] > 0 ) {
							splice2[i] = new int[anz[i]][parts.length-i][];
							su = new int[anz[i]];
							st = new int[anz[i]];
							e = new int[anz[i]];
							qS = new int[anz[i]];
							qE = new int[anz[i]];
							qL = new int[anz[i]];
							
							int idx = 0, idx2;
							for( int j = 0; used[i] != null && j < used[i].length; j++ ) {
								if( used[i][j] >= 0 ) {
									su[idx]=sums[i][j];
									st[idx]=start[i][j];
									e[idx]=end[i][j];
									qS[idx]=qStart[i][j];
									qE[idx]=qEnd[i][j];
									qL[idx]=this.qL[i][j];
									for( int k = i; k < used.length; k++ ) {
										splice2[i][idx][k-i] = new int[anz[k]];
										if( splice[i][j][k-i] != null && anz[k]>0 ) {
											idx2=0;
											for( int l = 0; l < used[k].length; l++ ) {
												if( used[k][l] >= 0 ) {
													splice2[i][idx][k-i][idx2++] = splice[i][j][k-i][l];
												}
											}
										}
									}
									idx++;
								}
							}					
							sums[i] = su;
							start[i] = st;
							end[i] = e;
							qStart[i] = qS;
							qEnd[i] = qE;
							this.qL[i] = qL;
						} else {
							sums[i] = start[i] = end[i] = qStart[i] = qEnd[i] = null;
						}
					}
					splice = splice2;
				}				
			}/**/
			
			
			//remove unused blast hits
			HashMap<Integer,ArrayList<Hit>> filtered = new HashMap<Integer,ArrayList<Hit>>();
			for( int i = 0; i < used.length; i++ ) {
				if( used[i] != null ) {
					list = lines.get(parts[i]);
					ArrayList<Hit> add = new ArrayList<Hit>();
					for( int j = 0; j < used[i].length; j++ ) {
						if( used[i][j] >= 0 ) {
							if( avoidStop ) {
								list.get(j).addSplitHits( add );
							} else {
								add.add(list.get(j));
							}
						}
					}
					if( add.size() > 0 ) {
						filtered.put( parts[i], add );
					}
				}
			}
			used=null;
			
			return filtered;
		}
		
		private Solution detailedAnalyseRegion( String chromosome, boolean forward, HashMap<Integer,ArrayList<Hit>> filtered, boolean backup ) throws CloneNotSupportedException, WrongAlphabetException, IOException {
			if( verbose ) show(filtered);
			ArrayList<Hit> current;
			
			//prepare
			for( int i = 0; i < parts.length; i++ ) {
				current = filtered.get(parts[i]);
				oldSize[i] = current == null ? 0 : current.size();
			}
			
			String chr = seqs.get(chromosome);
			if( acceptorSites!= null && donorSites != null ) {
				for( int h = 0; h <parts.length; h++ ) {
					ArrayList<Hit> list = filtered.get(parts[h]);
					for( int k = 0; list!= null && k < list.size(); k++ ) {
						Hit hit = list.get(k);
						hit.prepareSpliceCandidates(chr);
					}
				}
			}
			
			//TASK 2: add missing CDS parts
			int first = -1, last=parts.length;
			int end, start, old = -1;
			for( int j=0; j<parts.length; j++ ) {
				current = filtered.get(parts[j]);
				if( current != null && current.size() > 0 ) {
					if( first < 0 ) {
						first = j;
					}
					last=j;
					if( old >= 0 && old+1 != j ) {//missing part
						ArrayList<Hit> previous = filtered.get(parts[old]);
						int offset = current.size(), d;
						UnionFind uf = new UnionFind(offset+previous.size());
						for( int a = 0; a < current.size(); a++ ) {
							Hit x = current.get(a);
							for( int b = 0; b < previous.size(); b++ ) {
								Hit y = previous.get(b);
								if( forward ) {
									d = x.targetStart - y.targetEnd;
									if( d >= 0 && d < (j-old)*MAX_INTRON_LENGTH ) {
										uf.union(a, offset+b);
									}
								} else {
									d = y.targetStart - x.targetEnd;
									if( d >= 0 && d < (j-old)*MAX_INTRON_LENGTH ) {
										uf.union(a, offset+b);
									}
								}
							}
						}
						int[][] array = uf.getComponents();
						
						for( int a = 0; a < array.length; a++ ) {
							boolean de=false, ae=false;
							//determine region						
							if( forward ) {
								start=-1;
								end = Integer.MAX_VALUE;
								for( int i = 0; i < array[a].length; i++ ) {
									if( array[a][i] < offset ) {
										start = Math.max( start, current.get(array[a][i]).targetStart );
										de |= current.get(array[a][i]).de;
									} else {
										end = Math.min( end, previous.get(array[a][i]-offset ).targetEnd );
										ae |= previous.get(array[a][i]-offset).ae;
									}
								}
							} else {
								start= Integer.MAX_VALUE;
								end = -1;
								for( int i = 0; i < array[a].length; i++ ) {
									if( array[a][i] < offset ) {
										start = Math.min( start, current.get(array[a][i]).targetEnd );
										ae |= current.get(array[a][i]).ae;
									} else {
										end = Math.max( end, previous.get(array[a][i]-offset).targetStart);
										de |= previous.get(array[a][i]-offset).de;
									}
								}
							}
							
							//align cds parts in this region to find hits
							if( (start > end ) == forward ) {
								if( verbose ) protocol.append("INTERNAL\t" + (old+1) + ".." + j + "\t" + parts[old+1] + ".." + parts[j] + "\n");
								align(chromosome, forward, end, start, old+1, j, filtered, "internal", ae && de, ae && de );
							}
						}
					}
					old = j;
				}
			}
			
			//add upstream CDS part candidates
			if( first > 0 ) {
				if( verbose ) protocol.append("FIRST\t" + first + "\t" + parts[first] + "\n");
				extend(chromosome, forward, filtered, first, true, "upstream");
			}
			//add downstream CDS part candidates
			if( last < parts.length-1 ) {
				if( verbose ) protocol.append("LAST\t" + last + "\t" + parts[last] + "\n");
				extend(chromosome, forward, filtered, last, false, "downstream");
			}
			
			//check
			boolean cut = false;
			for( int i = 0; i < parts.length; i++ ) {
				current = filtered.get(parts[i]);
				if( current != null && current.size()-oldSize[i] > 20000 ) {
					while( current.size() > oldSize[i] ) {
						current.remove( current.size()-1 );
					}
					cut=true;
				}
			}
			
			//TASK 3: find new solution
			//prepare splicing for DP
			for( int h = 0; h <parts.length; h++ ) {
				ArrayList<Hit> list = filtered.get(parts[h]);
				String r = geneName + (transcriptInfo==null? "":("_" + parts[h]));
				String ref = cds.get(r);
				for( int k = 0; list!= null && k < list.size(); k++ ) {
					Hit hit = list.get(k);
					if( !hit.splice ) {
						hit.prepareSpliceCandidates(chr);
					}
					if( ref != null && ((h==0 && ref.charAt(0)=='M') || (h+1==parts.length && ref.charAt(ref.length()-1)=='*')) ) {
						hit.setBorderConstraints(h==0 && ref.charAt(0)=='M', h+1==parts.length && ref.charAt(ref.length()-1)=='*');//Laura Kelly: partial
					}
				}
				/*
				if( list != null ) {
					if( list.size() == 1 ) {
						Hit hit = list.get(0);
						if( hit.accMaxScore < -100 ) {
							protocol.appendln(hit);
							protocol.appendln(hit.accMaxScore);
							
							//alignPart(chromosome, forward, endLast, forward?hit.targetStart:hit.targetEnd, parts[h], 0, hit.queryStart, region, list, "extra-internal");
							
							System.exit(1);
						}
					}
				}  else {
					protocol.appendln( "MISSING\t" + parts[h]  + "\t" + Arrays.toString(parts) );
					System.exit(1);
				}*/
			}
			
			//sort hits			
			sort(filtered, forward);
			
			//filter if RNA-seq data
//TODO
/*
			if( parts.length > 1 && acceptorSites != null && donorSites != null ) {
				int[][][] acc = acceptorSites.get(chromosome);
				int[][][] don = donorSites.get(chromosome);
				int z = forward ? 0 : 1;
				if( acc != null && acc[z][1].length > 0
					&& don != null && don[z][0].length > 0 ) {
					//show(filtered);
					
					used = new double[parts.length][];			
					for( int j=parts.length-1; j>=0; j-- ) {
						current = filtered.get(parts[j]);
						if( current != null && current.size() > 0 ) {
							used[j] = new double[current.size()];
						}
					}
					ArrayList<Hit> next;

					for( int i = 0; i < parts.length; i++ ) {
						current = filtered.get(parts[i]);
						if( current!= null && current.size()>0 ) {
							Arrays.fill(used[i], 0);
							
							for( int a = 0; a < current.size(); a++ ) {
								Hit hi = current.get(a);
								
								//acceptor
								if( hi.ae ) {
									outerloop: for( int x = i; x >= Math.max(i-1,0); x-- ) {
										next = filtered.get(parts[x]);
										if( next != null ) {
											//outerloop: 
											for( int j = 0; j < 3; j++ ) {
												for( int b = 0; b < hi.accCand[j].length(); b++ ) {
													int pos = ((forward?hi.targetStart:(hi.targetEnd+1))+(forward?-1:1)*hi.accCand[j].get(b));
													
													for( int c = (x==i ? a: next.size())-1; c >= 0; c-- ) {
														Hit hi2 = next.get(c);
														if( hi2.de && hi2.postAcceptor.contains(pos) ) {
															used[i][a]++;
															break outerloop;
														}
													}
												}
											}
										}
									}
								}
								
								//donor
								if( hi.de ) {
									outerloop: for( int x = i; x <= Math.min(i+1,parts.length-1); x++ ) {
										next = filtered.get(parts[x]);
										if( next != null ) {
											//outerloop:
											for( int j = 0; j < 3; j++ ) {
												for( int k = 0; k < hi.donCand.length; k++ ) {
													for( int b = 0; b < hi.donCand[k][j].length(); b++ ) {
														int pos = ((forward?hi.targetEnd+1:hi.targetStart)-(forward?-1:1)*hi.donCand[k][j].get(b));
														
														for( int c = x==i ? a+1 : 0; c < next.size(); c++ ) {
															Hit hi2 = next.get(c);
															if( hi2.ae && hi2.preDonor.contains(pos) ) {
																used[i][a]++;
																break outerloop;
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
//System.out.println();		
					for( int i = 0; i < parts.length; i++ ) {
						current = filtered.get(parts[i]);
						if( current!= null && current.size()>0 ) {
							int t = (i==0 || i+1 == parts.length) ? 1 : 2;
							if( ToolBox.max(used[i])>=t ) {
//System.out.println("filter " + parts[i] );
								for( int aa = 0, a = 0; a < current.size(); a++ ) {
									if( used[i][a] >= t ) {
											aa++;
									} else {
										current.remove(aa);
									}
								}
							}
						}
					}				
/*
					
					int[] u = new int[parts.length];
					evidence = new boolean[parts.length][MAX_GAP+1];
					for( int j=parts.length-1; j>=0; j-- ) {
						current = filtered.get(parts[j]);
						if( current != null && current.size() > 0 ) {
							for( int a = current.size()-1; a >= 0; a-- ) {
								Hit hi = current.get(a);
								if( hi.ae ) {
									for( int i = 0; i < 3; i++ ) {
										for( int b = 0; b < hi.accCand[i].length(); b++ ) {
											int pos = ((forward?hi.targetStart:(hi.targetEnd+1))+(forward?-1:1)*hi.accCand[i].get(b));

											int m=Math.max(j-5,0);
											for( int k=j; k>=m; k-- ) {
												next = filtered.get(parts[k]);
												if( next != null ) {
													for( int c = k==j?a-1:next.size()-1; c >= 0; c-- ) {
														Hit hi2 = next.get(c);
														if( hi2.de && hi2.postAcceptor.contains(pos) ) {
															used[j][a]++;
															used[k][c]++;
															u[j]++;
															u[k]++;
															evidence[k][j-k]=true;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
					for( int j=parts.length-1; j>=0; j-- ) {
//if(verbose) System.out.println(j + "\t" + Arrays.toString(evidence[j]));
						if( u[j] != 0 && used[j] != null ) {
							current = filtered.get(parts[j]);
							for( int a = current.size()-1; a >= 0; a-- ) {
								if( used[j][a] == 0 ) {
									current.remove(a);
								}
							}
						}
					}
/**//*
					used=null;				
				}
			}

			if( verbose ) 
			{
				protocol.append("After discarding based on evidence\n");
				show(filtered);
			}
*/			
			//DP with splicing
			splice = new int[parts.length][][][];
			int bestValue = forwardDP(forward, filtered, true);
			filtered = reduce(chromosome, forward, filtered, avoidStop, true, bestValue);
			if( verbose ) {			
				for( int h = 0; h <parts.length; h++ ) {
					ArrayList<Hit> list = filtered.get(parts[h]);
					protocol.append( parts[h] + "\t" + Arrays.toString(sums[h]) + "\n" );
					for( int k = 0; list!= null && k < list.size(); k++ ) {
						Hit hi = list.get(k);
						protocol.append( hi + "\t" + hi.firstAddScore + "\t" + hi.lastAddScore + "\t" + hi.up + "\t" + hi.seqUp + "\t" + hi.down + "\t" + hi.seqDown + "\n" );
					}
				}
			}			
			
			//TASK 4: get best solution
			getBest(true, filtered);
			if( verbose ) {
				protocol.append( bestValue + "\n" );
				protocol.append( bestIdx[0] + ", " + bestIdx[1] +"\n" );
			}
			Solution best = new Solution( backup, cut );
			best.clear(forward);
			for( int h = 0; h < bestIdx[0].length(); h++ ) {
				int i = bestIdx[0].get(h);
				int idx = bestIdx[1].get(h);
				sol.clear(forward);
				
				if( verbose ) protocol.append( i + ", " + idx + "\n" );
				backTracking( forward, filtered, i, idx, 0, sol, best, true );
				if( verbose ) protocol.append( best + "\n" );
			}
			splice = null;
			best.setScore( bestValue );
			return best;
		}
		
		private void show( IntList[] il, String n ) {
			protocol.append( n +": " );
			for( int i = 0; i < il.length; i++ ) {
				protocol.append( il[i] + ", " );
			}
			protocol.append("\n");
		}
		
		private void extend( String chromosome, boolean forward, HashMap<Integer,ArrayList<Hit>> lines, int index, boolean upstream, String info ) throws IllegalArgumentException, WrongAlphabetException {
			int next, dir, end;
			if( upstream ) {
				dir = -1;
				end = -1;
			} else {
				dir=1;
				end = parts.length;
			}
			int d = MAX_INTRON_LENGTH/10, f=0;
			int startIndex, a = (upstream?0:d), stop;
			ArrayList<Hit> current;
			int[][] array = null;
			do {
				current = lines.get(parts[index]);
				if( current != null ) {
					array = new int[current.size()][2];
					for( int i = 0; i < array.length; i++ ) {
						Hit h = current.get(i);
						array[i][0] = upstream == forward ? h.targetStart : h.targetEnd;
						array[i][1] = (upstream ? h.ae : h.de) ? 1 : 0;
					}
					Arrays.sort(array,IntArrayComparator.comparator[0]);
					f = 1;
				} else {
					f++;
				}
				next = index+dir;
				if( upstream ) {
					stop = index;
				} else{
					stop = next+dir;
				}
				
				startIndex=0;
				boolean site = false;
				for( int i = 1; i <array.length; i++ ) {
					if( array[i][0]-array[i-1][0] > d ) {
						if( forward ) {
							align(chromosome, forward, array[startIndex][0]-f*(d-a), array[i-1][0]+f*a, next, stop, lines, info, upstream ? false : site, upstream ? site : false );
						} else {
							align(chromosome, forward, array[i-1][0]+f*(d-a), array[startIndex][0]-f*a, next, stop, lines, info, upstream ? false : site, upstream ? site : false );
						}
						startIndex=i;
						site = false;
					}
					site |= array[i][1] == 1;
				}
				if( forward ) {
					align(chromosome, forward, array[startIndex][0]-f*(d-a), array[array.length-1][0]+f*a, next, stop, lines, info, upstream ? false : site, upstream ? site : false );
				} else {
					align(chromosome, forward, array[array.length-1][0]+f*(d-a), array[startIndex][0]-f*a, next, stop, lines, info, upstream ? false : site, upstream ? site : false );
				}
				index=next;
			} while( index+dir != end );
		}
				
		private void align( String chromosome, boolean forward, int endLast, int startNext, int startIdx, int endIdx, HashMap<Integer,ArrayList<Hit>> lines, String info, boolean up, boolean down ) throws IllegalArgumentException, WrongAlphabetException {
			String chr = seqs.get(chromosome);
			String region;
			if( forward ) {
				if( endLast < 0 ) {
					endLast = 0;
				}
				if( startNext > chr.length() ) {
					startNext = chr.length();
				}
			} else {
				if( startNext < 0 ) {
					startNext = 0;
				}
				if( endLast > chr.length() ) {
					endLast = chr.length();
				}
			}
			//protocol.appendln( chromosome + "\t" + forward + "\t" + endLast + "\t" + startNext + "\t" + (endLast<startNext) );
			region = chr.substring(forward?endLast:startNext, forward?startNext:endLast);
			if( !forward ) {
				region = Tools.rc(region);
			}
			for( int i = startIdx; i < endIdx; i++ ) {//for each missing parts
				alignPart(chromosome, forward, endLast, startNext, i, -1, -1, region, lines, info, up, down);
				//protocol.append( "added " + parts[i] + "\t" + lines.get(parts[i]).size() );
			}
		}
		
		private boolean testForAlignment( boolean forward, 
				boolean a, int[][] acc, int accFrom, int accTo,
				boolean d, int[][] don, int donFrom, int donTo,
				int start, int end ) 
		{
			if( !a && !d ) {
				return true;
			}
			
			int help = -1;
			int middle, last;
			
			if( forward ) {
				middle = start;
				if( a && acc != null ) {
					help = Arrays.binarySearch(acc[1], accFrom, accTo, start );
					if( help < 0 ) {
						help = -(help+1);
					}
					if( help < acc[1].length ) {
						middle = acc[1][help];
					} else {
						middle = Integer.MAX_VALUE-1;
					}
				}
				last = middle;
				if( d && don != null ) {
					help = Arrays.binarySearch(don[0], donFrom, donTo, middle );
					if( help < 0 ) {
						help = -(help+1);
					}
					if( help < don[0].length ) {
						last = don[0][help];
					} else {
						last = Integer.MAX_VALUE;
					}
				}
				if( verbose ) {
					protocol.append( a + "\t" + start + "\t" + middle + ", " + last + "\t" + end+ "\t" + d + "\t\t= " + (last <= end) + "\n");
				}
				
				return last <= end;
			} else {
				/*TODO
				middle = end;
				if( a && acc != null ) {
					help = Arrays.binarySearch(acc[1], accFrom, accTo, start );
					if( help < 0 ) {
						help = -(help+1);
						help--;
					}
					if( help >= 0 ) {
						middle = acc[1][help];
					} else {
						middle = Integer.MIN_VALUE-1;
					}
				}
				last=Integer.MIN_VALUE;
				
				
				*/
				return true;//last >= start //FIX ME
			}
		}
		
		// @param startPos, endPos &lt; 0 leads to a complete search 
		private int alignPart( String chromosome, boolean forward, int endLast, int startNext, int idx, int startPos, int endPos, String region, HashMap<Integer,ArrayList<Hit>> lines, String info, boolean upstream, boolean downstream ) throws IllegalArgumentException, WrongAlphabetException {
			int[][][] sites;
			int[][] acc = null, don = null;
			int donFrom=-1, donTo=-1, accFrom=-1, accTo=-1;
			upstream = downstream = false;//XXX
/*
			if( acceptorSites != null ) {
				sites = acceptorSites.get(chromosome);
				if( sites != null )  {
					acc = sites[forward?0:1];
					
					accFrom = Arrays.binarySearch( acc[1], forward ? endLast : startNext );
					if( accFrom < 0 ) {
						accFrom = -(accFrom+1);
					}

					accTo = Arrays.binarySearch( acc[1], forward ? startNext : endLast );
					if( accTo < 0 ) {
						accTo = -(accTo+1);
					}
				}
			}

			if( donorSites != null ) {
				sites = donorSites.get(chromosome);
				if( sites != null )  {
					don = sites[forward?0:1];
					
					donFrom = Arrays.binarySearch( don[0], forward ? endLast : startNext );
					if( donFrom < 0 ) {
						donFrom = -(donFrom+1);
					}

					donTo = Arrays.binarySearch( don[0], forward ? startNext : endLast );
					if( donTo < 0 ) {
						donTo = -(donTo+1);
					}
				}
			}
/**/
			Sequence cdsSeq = null;
			String r = geneName + (transcriptInfo==null? "":("_" + parts[idx]));
			String m = cds.get(r);
				
			//collect candidate regions/alignments
			double best = 0;
			int anz = 0;
			if( m != null ) {
				int len = m.length();
				boolean findStart = idx==0 && startPos<=0 && m.charAt(0)=='M'; //Laura Kelly: partial?
				if( startPos >= 0 && endPos>= 0 ) {
					m = m.substring(startPos, endPos);
				}
				
				ArrayList<Hit> cand = lines.get(parts[idx]);
				boolean add = cand == null;
				if( add ) {
					cand = new ArrayList<Hit>();
				}
				int oldSize = cand.size();
				
				cdsSeq = Sequence.create(alph, m);
				int tested = (upstream||downstream)?0:1;
				do {
					//from all reading frames of this strand
					for( int off, splitStart, l, k = 0; k < 3; k++ ) {
						splitStart=k;
						//split the complete region at STOP codons
						String trans = Tools.translate(k, region, code, false, ambiguity);
						//protocol.appendln(trans);
						String[] split = trans.split("\\*");
						for( int a = 0; a < split.length; a++ ) {
							if( idx+1==parts.length && a+1 < split.length ) {
								split[a] +="*";
								off=0;
							} else {
								off=(a+1 < split.length) ? 1 : 0;
							}
							
							int e, x;
							if( forward ) {
								e = endLast + (forward?1:-1)*(splitStart-2);
								x = e + (forward?1:-1)*((split[a].length()+off)*3-1);
							} else {
								x = endLast + (forward?1:-1)*(splitStart-2);
								e = x + (forward?1:-1)*((split[a].length()+off)*3-1);
							}
							
							if( split[a].length()>0 //split not empty
									&& (!findStart || m.length() > 2*missingAA || split[a].indexOf('M')>=0 ) //either not first part or long enough or contains M
									//XXX && testForAlignment( forward, upstream, acc, accFrom, accTo, downstream, don, donFrom, donTo, e, x)//TODO splicing?
							) {
								//do an optimal local alignment
								Sequence intronPartSeq = Sequence.create(alph, split[a]);
								PairwiseStringAlignment sa = align.getAlignment(AlignmentType.LOCAL,cdsSeq, intronPartSeq);
								double c = sa.getCost();
								if( c<0 && c <= best*hitThreshold ) {//keep it if it is good enough
									//position within the genome/chromosome/contig
									//TODO das ist doof!!!!!
									String xx = sa.getAlignedString(1).replaceAll("-","");
									l = 3*xx.length();
									int z=split[a].indexOf('M');
									int endAlign2 = -1;
									
									while( (endAlign2 = split[a].indexOf( xx, endAlign2+1 )) >= 0 ) {
										//TODO check
										if( !findStart || sa.getStartIndexOfAlignmentForFirst() > missingAA || (z>=0 && z<endAlign2+xx.length()) )
										{										
											int pos, start, end;
											if( forward ) {
												pos = endLast + (splitStart+3*endAlign2);
												start = pos+1;
												end = pos+l;
											} else {
												pos = endLast-(splitStart+3*endAlign2);
												end = pos-l+1;
												start = pos;
											}

											if( (int)start<0 || (int)end < 0 ) {
												protocol.append(endLast+"\n");
												protocol.append(startNext+"\n");
												throw new RuntimeException("negative pos");
											}
											
											//System.out.println( r + "\t" + chromosome + "\t" + len + "\t" + start + "\t" + end + "\t" + info + "\t" + sa.getAlignedString(0) + "\t" + sa.getAlignedString(1));
											Hit h = new Hit(r, chromosome, 
													sa.getStartIndexOfAlignmentForFirst()+1, sa.getEndIndexOfAlignmentForFirst(), len,
													start, end, 
													(int)-sa.getCost(), sa.getAlignedString(0), sa.getAlignedString(1), info+";");
											cand.add( h );
											
											if( sa.getCost() < best ) {
												best = sa.getCost();
											}
										}
									}
								}
							}
							splitStart+=3*(split[a].length()+off);
						}
					}
	
					//remove 
					int j = cand.size()-1;
					while( j >= 0 && cand.size() > oldSize ) {
						Hit o = cand.get(j);
						if( o.score <= -hitThreshold*best ) {
							cand.remove(j);
						}
						j--;
					}
					if( add && cand.size() > oldSize ) {
						lines.put( parts[idx], cand );
					}
					anz = cand.size() - oldSize;
					tested++;
					
					//to ensure that next round each possible sequence is tested (in an alignment)
					upstream = downstream = false;
				} while( anz == 0 && (acc!=null && don!=null && tested < 2) );
			}
			return anz;
		}
		
		String getProteinSeqFromExons( Hit first, Hit second ) {
			String m;
			String firstCDS = cds.get( first.queryID );
			int partStart = first.part, partEnd = second.part;
			if( partStart == partEnd ) {
				m = firstCDS.substring(first.queryStart-1, second.queryEnd);				
			} else {
				m = firstCDS.substring(first.queryStart-1);
				for( int idx = revParts[partStart]+1; idx < revParts[partEnd]; idx++ ) {
					String p=cds.get(geneName + (transcriptInfo==null?"":("_" + parts[idx])) );
					if( p != null ) {
						m += p;
					}
				}
				m += cds.get( second.queryID ).substring(0, second.queryEnd);
			}
			return m;
		}
		
		boolean[] spliceType = new boolean[4];
		
		/**
		 * This method tries to find possible splice variant between the {@link Hit}s.
		 * 
		 * @param chr the chromosome
		 * @param first the first {@link Hit}
		 * @param second the second {@link Hit}
		 * 
		 * @return the score for the best splice variant between both {@link Hit}s
		 */
		private int checkSpliceSites( String chr, Hit first, Hit second, int delta)
				throws WrongAlphabetException {
			int d = revParts[second.part] - revParts[first.part], score;

			if( d <= 1 ) {
				boolean b = d!=0;
				Arrays.fill(spliceType, b);
				spliceType[3] = !b;
					
				/*
				Arrays.fill(spliceType, false);
				if( d == 0 ) {
					spliceType[3] = true;
				} else {
					if( phase != null ) {
						spliceType[phase[revParts[second.part]]] = true;
					}
				}/**/
				
				//find splice variant that has same exon phase as the reference
				score = checkSpecificSpliceSiteTypes(chr, first, second, delta);
				if( score == NO_SPLICE_VARIANT ) { /// if not found, search any splice variant
					for( int i = 0; i < spliceType.length; i++ ) {
						spliceType[i] = !spliceType[i];
					}
					score = checkSpecificSpliceSiteTypes(chr, first, second, delta);
				}
			} else {
				Arrays.fill(spliceType, true);
				score = checkSpecificSpliceSiteTypes(chr, first, second, delta);
			}
			return score;
		}
		
		
		private int checkSpecificSpliceSiteTypes( String chr, Hit first, Hit second, int delta ) throws WrongAlphabetException {	
			refined[0] = refined[1]= 0;
			boolean f = first.forward;
			int best=NO_SPLICE_VARIANT, current;
			
			Sequence targetSeq, cdsSeq = null;
			String region = null;
						
			//"intron-loss"
			int diff = first.forward ? (second.targetStart-1 - first.targetEnd) : (first.targetStart-1 - second.targetEnd);
			if( spliceType[3] && diff >= 0 
					&& diff % 3 == 0 //in-frame
					&& (
							( first.forward && first.maxDownstreamEnd > second.targetStart && first.targetEnd > second.maxUpstreamStart )
							||( !first.forward && first.maxDownstreamEnd < second.targetEnd && first.targetStart < second.maxUpstreamStart)
						) //there is no STOP codon
			) {
				if( first.forward ) {
					region = chr.substring( first.targetStart-1, second.targetEnd );
				} else {
					region = chr.substring( second.targetStart-1, first.targetEnd );
					region = Tools.rc(region);
				}
				targetSeq = Sequence.create(alph, Tools.translate(0, region, code, false, ambiguity));
				cdsSeq = Sequence.create(alph, getProteinSeqFromExons(first, second));
								
				align.computeAlignment(AlignmentType.GLOBAL, targetSeq, cdsSeq );
				current = (int) align.getCost(targetSeq.getLength(), cdsSeq.getLength());
				best = -current-(first.score+second.score);
				refined[0] = diff;
				refined[1] = 0;
				
				best-=(revParts[second.part] - revParts[first.part])*INTRON_GAIN_LOSS;
			}			

			boolean same = second.part == first.part, computed = false;			
			
			if( (!same || best == NO_SPLICE_VARIANT) && (spliceType[0] || spliceType[1] || spliceType[2]) ) {
				int gap = getGapCost( delta, revParts[second.part] - 1 - revParts[first.part] );			
				int gain_or_loss = Math.abs(revParts[second.part] - 1 - revParts[first.part]) * INTRON_GAIN_LOSS;
				
				int a, b, remainingFirst, p = 0, d; 
				int length=Integer.MAX_VALUE, len;
				int end = f ? first.targetEnd : first.targetStart;
				int start = f ? second.targetStart : second.targetEnd;
				String remaining1, remaining2=null, triplett;
				char aa;
				boolean set=false;
				do {
					for( int remainingSecond = 0; remainingSecond < 3; remainingSecond++ ) {
						if( spliceType[remainingSecond] ) {
							remainingFirst = (3-remainingSecond)%3;
							//protocol.appendln( DONOR[p] + "\t" + remainingSecond + "\t" + second.accCand[remainingSecond] + "\t" + first.donCand[p][remainingFirst]);
							for( int j = 0; j < second.accCand[remainingSecond].length(); j++ ) {
								a = second.accCand[remainingSecond].get(j);
								if( remainingSecond != 0 ) {
									current = f ? (second.targetStart-1 - a) : (second.targetEnd + a - remainingSecond);
									remaining2 = chr.substring( current, current + remainingSecond );
									if( !f ) {
										remaining2 = Tools.rc(remaining2);
									}
								}
								for( int k = 0; k < first.donCand[p][remainingFirst].length(); k++ ) {
									b = first.donCand[p][remainingFirst].get(k);
									d = f ? (start-a - (end + b)) : (end-b-(start+a));
									if( d >= MIN_INTRON_LENGTH ) { //not overlapping after extension
										if( remainingSecond != 0 ) {
											current = f ? (first.targetEnd + b- remainingFirst) : (first.targetStart-1 - b );								
											remaining1 = chr.substring( current, current + remainingFirst );
											if( !f ) {
												remaining1 = Tools.rc(remaining1);
											}
	
											triplett = remaining1+remaining2;
											aa=Tools.translate(triplett, code, ambiguity);
										} else {
											remaining1 = remaining2 = triplett=null;
											aa='!';
										}
										
										if( /*triplett == null ||*/ aa != '*' ) {// no STOP
											if( !same ) {
												current = first.donCandScore[p][remainingFirst].get(k)
														+ second.accCandScore[remainingSecond].get(j)
														+ gap;
												//TODO possible double or triple gap opening cost
											} else {
												//gap cost are computed in the alignment
												if( !approx ) {
													//slow
													if( cdsSeq == null ) {
														cdsSeq = Sequence.create(alph, getProteinSeqFromExons(first, second) );		
													}
													if( first.forward ) {
														region = chr.substring( first.targetStart-1, first.targetEnd+b )
															+ chr.substring( second.targetStart-1-a, second.targetEnd );
													} else {
														region = chr.substring( second.targetStart-1, second.targetEnd+a )
																+ chr.substring( first.targetStart-1-b, first.targetEnd );
														region = Tools.rc(region);
													}
													targetSeq = Sequence.create(alph, Tools.translate(0, region, code, false, ambiguity));
													
													align1.computeAlignment(AlignmentType.GLOBAL, targetSeq, cdsSeq );
													current = (int) align1.getCost(targetSeq.getLength(), cdsSeq.getLength()) - gain_or_loss;
												} else {
													if( !computed ) {
														if( cdsSeq == null ) {
															cdsSeq = Sequence.create(alph, getProteinSeqFromExons(first, second) );		
														}
														align1.computeAlignment(AlignmentType.GLOBAL, cdsSeq, Sequence.create(alph, first.targetAlign.replaceAll("-", "")+first.down) );
														try {
															align2.computeAlignment(AlignmentType.GLOBAL, cdsSeq.reverse(), Sequence.create(alph, second.up + second.targetAlign.replaceAll("-", "")).reverse() );
														} catch( OperationNotSupportedException onse ) {
															//does not happen
															throw new RuntimeException( onse.getMessage() );
														}
														computed = true;
													}											
													//fast										
													region = "";
													int l1 = first.targetEnd+b - (first.targetStart-1);
													int l2 = second.targetEnd - (second.targetStart-1-a);
													double v = MyAlignment.getBestValue(align1, align2, l1/3, l2/3, aa);
													if( Double.isInfinite(v) ) {
														current=NO_SPLICE_VARIANT;
													} else {
														current = (int) v;
													}
												}
												if( current!=NO_SPLICE_VARIANT ) {
													//costs for the splice variant equals the difference between the real costs and the precomputed costs of the individual parts
													current = (-current-gain_or_loss)//real costs (for the alignment)
															-first.score-second.score;//precomputed costs (of the individual parts) 
												}
											}
											
											len = Math.abs(a) + Math.abs(b);
											
											//protocol.appendln(a + "\t" + b + "\t" + current + "\t" + best + "\t" + len + "\t" + length + "\t" + triplett + "\t" + (triplett==null?"?":aa) );
											
											if( current != NO_SPLICE_VARIANT && 
													((current > best)  || (current == best && len < length)) //best solution so far
											) {
												best = current;
												length = len;											
	
												refined[0]=first.donCand[p][remainingFirst].get(k);
												refined[1]=second.accCand[remainingSecond].get(j);
												set=true;
											}
										}										
									}
								}
							}
						}
					}
					p++;
				}while( p < first.donCand.length && !set );
			}			
			return best;
		}
		
		IntList[] bestIdx = new IntList[]{new IntList(), new IntList()}; 
		
		int getBest( boolean gap, HashMap<Integer,ArrayList<Hit>> hits ) {
			bestIdx[0].clear();
			bestIdx[1].clear();
			int bestValue = Integer.MIN_VALUE;
			
			for( int i = 0; i < sums.length; i++ ) {
				ArrayList<Hit> current = hits.get(parts[i]);
				//if( (i == 0 && sums[i] != null ) || (i>0 && bestValue <= maxs[i] + cumLength[i]) ) 
				{
					for( int j = 0; sums[i] != null && j < sums[i].length; j++ ) {
						int add = getCost( current.get(j), gap, true, i );
						
						if( bestValue <= sums[i][j] + add ) {
							if( bestValue < sums[i][j] + add ) {
								bestIdx[0].clear();
								bestIdx[1].clear();
								bestValue = sums[i][j] + add;
							}
							bestIdx[0].add(i);
							bestIdx[1].add(j);
						}
					}
				}
			}
			return bestValue;
		}
		
		/**
		 * This class represents a predicted gene model.
		 * 
		 * @author Jens Keilwagen
		 */
		class Solution implements Comparable<Solution>, Cloneable {
			LinkedList<Hit> hits;
			boolean forward, backup, cut;
			int score, a, d, A, D, i, I, minSplitReads;
			int Len, Cov, minC, Sum;
			
			public Solution() {
				this( false, false );
			}
			
			public Solution( boolean b, boolean c ) {
				hits = new LinkedList<Hit>();
				score = Integer.MIN_VALUE;
				backup = b;
				cut = c;
			}

			public void setScore( int score ) {
				this.score = score;
			}

			public int getNumberOfAA() {
				int l = 0;
				Iterator<Hit> it = hits.iterator();
				while( it.hasNext() ) {
					l+=it.next().getLength();
				}
				return l / 3;
			}
			
			public String getDNA() {
				StringBuffer sb = new StringBuffer();
				Iterator<Hit> it = hits.iterator();
				while( it.hasNext() ) {
					sb.append(it.next().getDNA(0,0));
				}
				return sb.toString();
			}
			
			public String getProtein() {
				return Tools.translate(0, getDNA(), code, false, ambiguity);
			}

			public Solution clone() throws CloneNotSupportedException {
				Solution clone = (Solution) super.clone();
				clone.hits = new LinkedList<Hit>();
				Iterator<Hit> it = hits.iterator();
				while( it.hasNext() ) {
					clone.hits.add( it.next().clone() );
				}
				return clone;
			}
			
			public void clear(boolean f) {
				hits.clear();
				forward = f;
				a=d=A=D=i=I=0;
			}
			
			public int matchParts() {
				int p = 0;
				boolean[] match = new boolean[parts.length];
				Arrays.fill(match, false);
				for( int i = 0; i < hits.size(); i++ ) {
					int part = hits.get(i).part;
					if( !match[revParts[part]] ) {
						match[revParts[part]]=true;
						p++;
					}
				}
				return p;
			}
			
			public int getMin() {
				int min = Integer.MAX_VALUE;
				for( int i = 0; i < hits.size(); i++ ) {
					Hit o = hits.get(i);
					int v = forward ? o.targetStart : o.targetEnd;
					if( v < min ) {
						min=v;
					}
				}
				return min;
			}
			
			public int getMax() {
				int max = Integer.MIN_VALUE;
				for( int i = 0; i < hits.size(); i++ ) {
					Hit o = hits.get(i);
					int v = forward ? o.targetStart : o.targetEnd;
					if( v > max ) {
						max=v;
					}
				}
				return max;
			}
			
			public int getDifference() {
				return getMax()-getMin();
			}
		
			public void set( Solution s ) {
				clear( s.forward );
				hits.addAll( s.hits );
			}
			
			public int compareTo(Solution o) {//TODO check: do something better
				int d = similar( o );
				if( d == 0 ) {
					d = getDifference() - o.getDifference();
				}
				return d;
			}
			
			public int similar(Solution o) {
				//if something is undefined
				if( o == null ) {
					return Integer.MIN_VALUE;
				}
				if( score == Integer.MIN_VALUE ) {
					return 1;
				}
				if( o.score == Integer.MIN_VALUE ) {
					return -1;
				}
				//first check: higher score?
				int d= o.score - this.score;
				if( d == 0 ) {
					//second check: higher number of matched parts?
					d = o.matchParts() - matchParts();
					if( d == 0 ) {
						//third check: smaller number of parts?
						d = hits.size() - o.hits.size();
						if( d == 0 ) {
							//fourth check: correct start?
							Hit h = hits.get(0);
							char firstAA='!', oFirstAA='!';
							if( revParts[h.part]==0 ) {
								firstAA = h.firstAA;
							}
							h = o.hits.get(0);
							if( revParts[h.part]==0 ) {
								oFirstAA = h.firstAA;
							}
							if( firstAA != oFirstAA ) {
								if( firstAA == 'M' ) {
									d=1;
								} else if( oFirstAA == 'M' ) {
									d=-1;
								}
							}
						}
					}
				}
				return d;
			}
			
			public String toString() {
				Iterator<Hit> it = hits.iterator();
				String res = "";
				while( it.hasNext() ) {
					res += it.next() + "\n";
				}
				return res;
			}
			
			public Solution refine( String transcriptName ) throws WrongAlphabetException, CloneNotSupportedException {
				if( hits.size() == 0 ) {
					return this;
				}
				//refine = add start codon, splice sites, and stop codon to the annotation
				Solution res = new Solution();
				res.clear(forward);
				res.setScore(score);

				Hit l = hits.get(0), c;
				res.hits.add(l.clone());
				
				String chr = seqs.get(l.targetID);
				//splicing				
				for( int i = 1; i < hits.size(); i++ ) {//iterate over found parts
					c = hits.get(i);
					res.hits.add(c.clone());
					int delta = 0;
					for( int j = revParts[l.part]+1; j < revParts[c.part]; j++ ) {
						delta += length[j];
					}
					if( checkSpliceSites(chr, l, c, delta/*, true*/ )!= NO_SPLICE_VARIANT ) {
						res.hits.get(i-1).setSpliceSite(true, refined[0]);
						res.hits.get(i).setSpliceSite(false, refined[1]);
					}
					l=c;
				}
				
				//Laura Kelly: partial?
				String prot = protein == null ? null : protein.get(transcriptName), ref;
				if( prot == null ) {
					String r = geneName + (transcriptInfo==null? "":("_" + parts[0]));
					ref = cds.get(r);
				} else {
					ref = prot; 
				}				
				if( ref==null || ref.charAt(0)=='M' ) {
					c = res.hits.get(0);
					c.extend( c.part == parts[0], false ); //START
				}
				if( prot == null ) {
					String r = geneName + (transcriptInfo==null? "":("_" + parts[parts.length-1]));
					ref = cds.get(r);
				} else {
					ref = prot; 
				}
				if( ref==null || ref.charAt(ref.length()-1)=='*' ) {
					c = res.hits.get(res.hits.size()-1);
					c.extend( false, c.part == parts[parts.length-1]); //STOP
				}
				return res;
			}
			
			public String getInfo(BufferedWriter gff) throws IOException {		
				Iterator<Hit> it = hits.iterator();
				Hit oldH = null, h;
				String dna = "";
				boolean intronGain=false, intronLoss = false;
				while( it.hasNext() ) {
					h=it.next();				
					if( oldH != null ) {
						if( h.part == oldH.part ) {
							if( forward ) {
								if( oldH.targetEnd+1<h.targetStart ) {
									intronGain = true;
								}
							} else {
								if( oldH.targetStart-1>h.targetEnd ) {
									intronGain = true;
								}
							}
						} else { //oldH.part != h.part ) {
							if( forward ) {
								if( oldH.targetEnd+1==h.targetStart ) {
									intronLoss = true;
								}
							} else {
								if( oldH.targetStart-1==h.targetEnd ) {
									intronLoss = true;
								}
							}
						}
					}
					dna += h.getDNA(0,0);
					oldH = h;
				}
				String protein = Tools.translate(0, dna, code, false, ambiguity);
				int idx = -1, anz=0;
				while( (idx=protein.indexOf('*',idx+1))>=0 ) {
					anz++;
				}
				
				gff.append(";start=" + protein.charAt(0) + ";stop=" + protein.charAt(protein.length()-1) );
				
				String res = (hits.getFirst().part == parts[0]) + "\t" + protein.charAt(0) + "\t" 
							+ (hits.getLast().part == parts[parts.length-1]) + "\t" + protein.charAt(protein.length()-1) + "\t" 
							+ anz + "\t" + intronGain +"\t" + intronLoss;
		
				stats[0]++;
				stats[1]+=protein.charAt(0)=='M'?1:0;
				stats[2]+=protein.charAt(protein.length()-1)=='*'?1:0;
				stats[3]+=intronGain?1:0;
				stats[4]+=intronLoss?1:0;
				
				return res;
			}
			
			public void writeSummary( String geneName, String transcriptName, int i ) throws Exception {
				Hit first = hits.get(0);
				String chr = seqs.get(first.targetID);
				gff.append( first.targetID + "\tGeMoMa\t"+tag+"\t" );
				StringBuffer genomicRegion = new StringBuffer();
				int off = 300, s;
				
				if( forward ) {
					gff.append( first.targetStart + "\t" + hits.getLast().targetEnd );
					s = Math.max(first.targetStart-off,1);
					genomicRegion = new StringBuffer( chr.substring(s-1, Math.min( chr.length(), hits.getLast().targetEnd+off)).toLowerCase() );
				} else {
					gff.append( hits.getLast().targetStart + "\t" + first.targetEnd );
					s = Math.min(chr.length(), first.targetEnd+off);
					genomicRegion = new StringBuffer( Tools.rc(chr.substring(Math.max(hits.getLast().targetStart-off,1)-1, 
							s )).toLowerCase() );
				}

				gff.append( "\t.\t" + (forward?"+":"-") + "\t.\tID=" + prefix+transcriptName + "_R" + i + ";ref-gene=" + geneName + ";AA="+getNumberOfAA() );
				if( writeGenomic ) {//> v1.3.1
					genomic.append(">" + prefix+transcriptName + "_R" + i );
					genomic.newLine();
					for( Hit t : hits ) {
						int a = t.targetStart-s, b = t.targetEnd-s;
						if( forward ) {
							genomicRegion.replace(a, b+1, genomicRegion.substring(a, b+1).toUpperCase() );
						} else {
								genomicRegion.replace(-b, -a+1, genomicRegion.substring(-b, -a+1).toUpperCase() ); //TODO Bug Sven
						}
					}
					genomic.append(genomicRegion);
					genomic.newLine();
				}
			}
			
			public int writeGFF( String transcriptName, int pred, StringBuffer sb ) throws Exception {
				a=d=A=D=i=I=0;
				minSplitReads = Integer.MAX_VALUE;
				
				int start=-1, end = -1, phase = -1, parts = 0;
				String id = null;
				boolean first = true, ae = false, de = false;
				int[][] donSites = null;
				int[][] cov = null;
				
				a = A = d = D = Len = Cov = Sum = 0;
				minC= Integer.MAX_VALUE;
				
				int last = -1;
				String pref = prefix+transcriptName + "_R" + pred;
				
				for( Hit t : hits ) {
					if( id == null ) {
						id = t.targetID;
						
						int[][][] sites = donorSites != null ? donorSites.get(id): null;
						if( sites != null ) {
							donSites = sites[forward?0:1];
						}
						cov = (coverage != null && coverage[forward?0:1]!= null) ? coverage[forward?0:1].get(id) : null;
					}
					
					if( start < 0 ) {
						start = t.targetStart;
						end = t.targetEnd;
						phase=t.phase;
						ae = t.ae;
						de = t.de;
					} else if ( forward && end+1 == t.targetStart ) {
						end = t.targetEnd;
						de = t.de;
					} else if( !forward && start-1 == t.targetEnd ) {
						start = t.targetStart;
						de = t.de;
					} else {
						int l = end-start+1;
						int covered = 0, min = Integer.MAX_VALUE;
						if( cov != null ) {
							int idx = Arrays.binarySearch(cov, new int[]{start}, CovComparator.def );
							if( idx < 0 ) {
								idx = -(idx+1);
								idx = Math.max(0, idx-1);
							}
							
							int[] inter = cov[idx];
							//System.out.println("hier " + start + " .. " + end + "\t" + Arrays.toString(inter) );
							int p = start;
							outerloop: while( p <= end ) {
								while( p > inter[1] ) {
									idx++;
									if( idx < cov.length ) {
										inter = cov[idx];
									} else {
										min = 0;
										break outerloop;
									}
								}
								if( inter[0]<= p && p <=inter[1] ) {
									int h = Math.min(inter[1],end)+1;
									int a=h-p;
									covered+=a;
									Sum+=inter[2] * a;
									min = Math.min(min, inter[2]);
									
									p=h;
								} else {//p<inter[0] && p<=inter[1]
									min = 0;
									
									p = Math.min(inter[0],end+1);
								}
								//System.out.println(p + "\t" + Arrays.toString(inter) + "\t" + covered + "\t" + min);
							}
							
							Len +=l;
							Cov +=covered;
							minC = Math.min(minC, min);
						}				
						
						sb.append( id + "\tGeMoMa\tCDS\t" + start + "\t" + end + "\t.\t" + (forward?"+":"-") + "\t" +phase+ "\tID=" +pref+"_cds"+parts+ ";Parent=" + pref 
								+ (acceptorSites==null || first?"":(";ae="+ae))
								+ (donorSites==null?"":(";de="+de))
								+ (cov == null?"":(";pc=" + decFormat.format(covered/(double)l)+";minCov=" + min))
								+ "\n"
						);
						
						D++;
						d += de ? 1 : 0;
						if( !first ) {
							A++;
							a += ae ? 1 : 0;
							if( donSites != null ) {
								I++;
								int v = forward ? start : (end+1);
								int idx = Arrays.binarySearch( donSites[0], last );
								if( idx >= 0 ) {
									while( idx < donSites[0].length && donSites[0][idx] == last && donSites[1][idx] != v ) {
										idx++;
									}
									if( idx < donSites[0].length && donSites[0][idx] == last && donSites[1][idx] == v ) {
										minSplitReads = Math.min(minSplitReads, donSites[2][idx] );
										i++;
									} else {
										minSplitReads = 0;
									}
								} else {
									minSplitReads = 0;
								}
							}
						}
						parts++;
						last = forward ? end+1 : start;
						start = t.targetStart;
						end = t.targetEnd;
						phase=t.phase;
						ae = t.ae;
						de = t.de;
						first = false;
					}
				}
				
				int l = end-start+1;
				int covered = 0, min = Integer.MAX_VALUE;
				if( cov != null ) {
					int idx = Arrays.binarySearch(cov, new int[]{start}, CovComparator.def );
					if( idx < 0 ) {
						idx = -(idx+1);
						idx = Math.max(0, idx-1);
					}
					
					int[] inter = cov[idx];
					//System.out.println("hier " + start + " .. " + end + "\t" + Arrays.toString(inter) );
					int p = start;
					outerloop: while( p <= end ) {
						while( p > inter[1] ) {
							idx++;
							if( idx < cov.length ) {
								inter = cov[idx];
							} else {
								min = 0;
								break outerloop;
							}
						}
						if( inter[0]<= p && p <=inter[1] ) {
							int h = Math.min(inter[1],end)+1;
							int a=h-p;
							covered+=a;
							Sum+=inter[2] * a;
							min = Math.min(min, inter[2]);
							
							p=h;
						} else {//p<inter[0] && p<=inter[1]
							min = 0;
							
							p = Math.min(inter[0],end+1);
						}
						//System.out.println(p + "\t" + Arrays.toString(inter) + "\t" + covered + "\t" + min);
					}
					
					Len +=l;
					Cov +=covered;
					minC = Math.min(minC, min);
				}
				
				sb.append( id + "\tGeMoMa\tCDS\t" + start + "\t" + end + "\t.\t" + (forward?"+":"-") + "\t" +phase+ "\tID=" +pref+"_cds"+parts+ ";Parent=" + pref
						+ (acceptorSites==null || first?"":(";ae="+ae))
						+ (cov != null?"":(";pc=" + decFormat.format(covered/(double)l)+";minCov=" + min))
						+ "\n"
				);
				if( !first ) {
					A++;
					a += ae ? 1 : 0;
					
					I++;
					if( donSites != null ) {
						int v = forward ? start : (end+1);
						int idx = Arrays.binarySearch( donSites[0], last );
						if( idx >= 0 ) {
							while( idx < donSites[0].length && donSites[0][idx] == last && donSites[1][idx] != v ) {
								idx++;
							}
							if( idx < donSites[0].length && donSites[0][idx] == last && donSites[1][idx] == v ) {
								minSplitReads = Math.min(minSplitReads, donSites[2][idx]);
								i++;
							}  else {
								minSplitReads = 0;
							}
						} else {
							minSplitReads = 0;
						}
					}
				}
				parts++;
								
				return parts;
			}
		}
	}
	
	static class MyAlignment extends Alignment {
		
		private static IntList first = new IntList(), second = new IntList();
		private GeMoMa instance;
		public MyAlignment(Costs costs,GeMoMa instance) {
			super(costs);
			this.instance = instance;
		}

		public static double getBestValue( MyAlignment a1, MyAlignment a2, int end1, int end2, char aa ) throws WrongAlphabetException {
			if( a1.type == AlignmentType.GLOBAL && a2.type == AlignmentType.GLOBAL //global alignments
					&& a1.startS1 == 0 && a2.startS2 == 0 //both from the beginning 
			) {
				first.clear();
				for( int i = 1; i < a1.l1; i++ ) {
					if( a1.d[0][i][end1] < a1.d[1][i][end1] && a1.d[0][i][end1] < a1.d[2][i][end1] ) {
						first.add(i);
					}
				}
				
				second.clear();
				for( int i = 1; i < a1.l1; i++ ) {
					if( a2.d[0][i][end2] < a2.d[1][i][end2] && a2.d[0][i][end2] < a2.d[2][i][end2] ) {
						second.add(i);
					}
				}
				int pos1, pos2;
				double best = Double.POSITIVE_INFINITY, current;
				for( int i = 0; i < first.length(); i++ ) {
					pos1 = first.get(i);
					for( int j = 0; j < second.length(); j++ ) {
						pos2=second.get(j);
						int l = (a1.l1-pos2)-pos1;//XXX ?
						if( l >= 0 ) {
							double cost = l> 0 ? a1.aCosts.getGapCostsFor(l) : 0;//approx;
							current=a1.d[0][pos1][end1] + a2.d[0][pos2][end2] + cost;
							//protocol.appendln(pos1 + "\t" + pos2 + "\t" + (a1.l1-pos2) + "\t" + l + "\t" + current + "\t" + best);
							if( current < best ) {
								/*
								StringAlignment sa1 = a1.getAlignment(new int[]{0,pos1,end1} );
								StringAlignment sa2 = a2.getAlignment(new int[]{0,pos2,end2} );
								protocol.appendln(a1.s1);
								protocol.appendln(sa1.getAlignedString(0) + "~" + new StringBuffer(sa2.getAlignedString(0)).reverse() );
								protocol.appendln(sa1.getAlignedString(1) + "~" + new StringBuffer(sa2.getAlignedString(1)).reverse());
								/**/
								best = current;
							}
						} else { 
							break;
						}
					}
				}
				return best;
			} else {
				throw new IllegalArgumentException( "The given alignments can not be used" );
			}
		}
		
		public boolean computeAlignment( AlignmentType type, Sequence s1, int startS1, int endS1, Sequence s2, int startS2, int endS2 ) {
			if( this.type == type 
				&& this.s1.toString().equals(s1.toString()) && this.startS1 == startS1 && this.l1 == endS1-startS1
				&& this.s2.toString().equals(s2.toString()) && this.startS2 == startS2 && this.l2 == endS2-startS2
					) {
				//same as last computed alignment
				return false;
			}
			instance.numberOfPairwiseAlignments++;
			return super.computeAlignment(type, s1, startS1, endS1, s2, startS2, endS2);
		}		
	}

	public ParameterSet getToolParameters() {
		try{
			return new SimpleParameterSet(
					new FileParameter( "tblastn results", "The sorted tblastn results", "tabular", true ),
					new FileParameter( "target genome", "The target genome file (FASTA), i.e., the target sequences in the blast run. Should be in IUPAC code", "fasta", true ),
					new FileParameter( "cds parts", "The query cds parts file (FASTA), i.e., the cds parts that have been blasted", "fasta", true ),
					new FileParameter( "assignment", "The assignment file, which combines parts of the CDS to transcripts", "tabular", false ),

					new FileParameter( "introns", "Introns (GFF), which might be obtained from RNAseq", "gff", false ),
					new SimpleParameter( DataType.INT, "reads", "if introns are given by a GFF, only use those which have at least this number of supporting split reads", true, new NumberValidator<Integer>(1, Integer.MAX_VALUE), 1 ),
					new SimpleParameter( DataType.BOOLEAN, "splice", "if no intron is given by RNAseq, compute candidate splice sites or not", true, true ),
					
					//TODO
					new SelectionParameter( DataType.PARAMETERSET, 
							new String[]{"NO", "UNSTRANDED", "STRANDED"},
							new Object[]{
								//no coverage
								new SimpleParameterSet(),
								//unstranded coverage
								new SimpleParameterSet(
										new FileParameter( "coverage_unstranded", "The coverage file contains the unstranded coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true )
								),
								//stranded coverage
								new SimpleParameterSet(
										new FileParameter( "coverage_forward", "The coverage file contains the forward coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true ),
										new FileParameter( "coverage_reverse", "The coverage file contains the reverse coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true )
								)
							},  "coverage", "experimental coverage (RNAseq)", true
					),
					
					/*splice site models
					new FileParameter( "classifier donor", "The path to the donor splice site classifier (XML)", "xml", false ),
					new FileParameter( "classifier acceptor", "The path to the donor splice site classifier (XML)", "xml", false ),					
					/**/
					new FileParameter( "query proteins", "optional query protein file (FASTA) for computing the optimal alignment score against complete protein prediction", "fasta", false ),
					new FileParameter( "genetic code", "optional user-specified genetic code", "tabular", false ),
					new FileParameter( "substitution matrix", "optional user-specified substitution matrix", "tabular", false ),
					
					new SimpleParameter( DataType.INT, "gap opening", "The gap opening cost in the alignment", true, 11 ),
					new SimpleParameter( DataType.INT, "gap extension", "The gap extension cost in the alignment", true, 1 ),
					new SimpleParameter( DataType.INT, "maximum intron length", "The maximum length of an intron", true, 15000 ),
					new SimpleParameter( DataType.INT, "intron-loss-gain-penalty", "The penalty used for intron loss and gain", true, 25 ),
			
					new SimpleParameter( DataType.DOUBLE, "e-value", "The e-value for filtering blast results", true, 1E2 ),
					new SimpleParameter( DataType.DOUBLE, "contig threshold", "The threshold for evaluating contigs", true, new NumberValidator<Double>(0d, 1d), 0.9 ),
					new SimpleParameter( DataType.DOUBLE, "region threshold", "The threshold for evaluating regions", true, new NumberValidator<Double>(0d, 1d), 0.9 ),
					new SimpleParameter( DataType.DOUBLE, "hit threshold", "The threshold for adding additional hits", true, new NumberValidator<Double>(0d, 1d), 0.9 ),
					
					new SimpleParameter( DataType.INT, "predictions", "The (maximal) number of predictions per transcript", true, 1 ), 
					new FileParameter( "selected", "The path to list file, which allows to make only a predictions for the contained transcript ids. The first column should contain transcript IDs as given in the annotation. Remaining columns can be used to determine a target region that should be overlapped by the prediction, if columns 2 to 5 contain chromosome, strand, start and end of region", "tabular,txt", maxSize>-1 ), 
					new SimpleParameter( DataType.BOOLEAN, "avoid stop", "A flag which allows to avoid stop codons in a transcript (except the last AS)", true, true ),
					new SimpleParameter( DataType.BOOLEAN, "approx", "whether an approximation is used to compute the score for intron gain", true, true ),
			
					new SimpleParameter( DataType.BOOLEAN, "align", "A flag which allows to output a tab-delimited file, which contains the results in a blast-like format (deprecated)", true, false ),
					new SimpleParameter( DataType.BOOLEAN, "genomic", "A flag which allows to output a fasta file containing the genomic regions of the predictions", true, false ),
					new SimpleParameter( DataType.STRING, "prefix", "A prefix to be used for naming the predictions", true, "" ),
					new SimpleParameter( DataType.STRING, "tag", "A user-specified tag for transcript predictions in the third column of the returned gff. It might be beneficial to set this to a specific value for some genome browsers.", true, "prediction" ),
					
					new SimpleParameter( DataType.BOOLEAN, "verbose", "A flag which allows to output wealth of additional information per transcript", true, false ),
					new SimpleParameter( DataType.LONG, "timeout", "The (maximal) number of seconds to be used for the predictions of one transcript, if exceeded GeMoMa does not ouput a prediction for this transcript.", true, new NumberValidator<Long>((long) 0, maxTimeOut), timeOut )
			);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	public String getToolName() {
		return "GeneModelMapper";
	}
	
	public String getToolVersion() {
		return "1.3.2";
	}
	
	public String getShortName() {
		return "GeMoMa";
	}

	public String getDescription() {
		return "builds gene models from tblastn results";
	}

	public String getHelpText() {
		return 
			"**What it does**\n\nThis tool is the main part of GeMoMa, a homology-based gene prediction tool. GeMoMa builds gene models from tblastn results.\n\n"
				//typical usage
				+ "As first step, you should run **Extractor** obtaining *cds parts* and *assignment*. Second, you should run **tblastn** with *cds parts* as query. Finally, these results are then used in **GeMoMa**.\n"
				//protein
				+ "If you like to run GeMoMa ignoring intron position conservation, you should blast protein sequences and feed the results in *query cds parts* and leave *assignment* unselected.\n\n"
				//multiple predictions
				+ "If you like to obtain multiple predictions per gene model of the reference organism, you should set *predictions* accordingly. In addition, we suggest to decrease the value of *contig threshold* allowing GeMoMa to evaluate more candidate contigs/chromosomes.\n\n"
				//runtime
				+ "If you change the values of *contig threshold*, *region threshold* and *hit threshold*, this will influence the predictions as well as the runtime of the algorithm. The lower the values are, the slower the algorithm is.\n\n"
			+ "**References**\n\nFor more information please visit http://www.jstacs.de/index.php/GeMoMa or contact jens.keilwagen@julius-kuehn.de.\n"
				+"If you use this tool, please cite\n\n*Using intron position conservation for homology-based gene prediction.*\n Keilwagen et al., NAR, 2016, http://nar.oxfordjournals.org/content/44/9/e89";
	}
	
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", "predicted annotation"),
				new ResultEntry(TextResult.class, "fasta", "predicted protein")
		};
	}
	
	private static class CovComparator implements Comparator<int[]> {

		static CovComparator def = new CovComparator();
		
		@Override
		public int compare(int[] o1, int[] o2) {
			return Integer.compare(o1[0], o2[0]);
		}
		
	}
}