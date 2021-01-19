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
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
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
import java.util.List;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.concurrent.Future;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import javax.naming.OperationNotSupportedException;

import de.jstacs.DataType;
import de.jstacs.algorithms.alignment.Alignment;
import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.PairwiseStringAlignment;
import de.jstacs.algorithms.alignment.StringAlignment;
import de.jstacs.algorithms.alignment.cost.AffineCosts;
import de.jstacs.algorithms.alignment.cost.MatrixCosts;
import de.jstacs.algorithms.graphs.UnionFind;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
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
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.tools.ui.galaxy.Galaxy;
import de.jstacs.utils.IntList;
import projects.gemoma.Tools.Ambiguity;

//Please check //XXX improve splice site scoring

/**
 * Predicting gene models from search hits (e.g. from tblastn) using a dynamic programming algorithm with extensions for splice sites, start and stop codon.
 * 
 * @author Jens Keilwagen
 */
public class GeMoMa extends GeMoMaModule {
	
	public static final String TAG = "mRNA";
	
	public static DecimalFormat decFormat = new DecimalFormat("###.####",DecimalFormatSymbols.getInstance(Locale.US));

	/**
	 * The index of the score in the blast output.
	 */
	private static final int scoreIndex = 13;

	public static final String timeoutMsg = "Invocation did not return before timeout";
	StringBuffer timeOutWarning;

	//global variables (always same for the same target species)
	//
	static HashMap<String, String> seqs;
	private static HashMap<String,Character> code;
	static HashMap<String, String> selected = null; //optional
	//splicing
	static HashMap<String, int[][][]> donorSites;
	static HashMap<String, int[][][]> acceptorSites;
	//coverage
	static HashMap<String, int[][]>[] coverage;	

	//fill
	static synchronized void fill( Protocol protocol, boolean verbose, int maxSize, String targetGenome, String selectedFile, int reads, ExpandableParameterSet introns, ExpandableParameterSet cov ) throws Exception {
		//read genome
		seqs = Tools.getFasta(targetGenome,20,".*");
		String x = Arrays.toString(seqs.keySet().toArray());
		if( x.length() > 200 ) {
			x = x.substring(0, 200) + "...";
		}
		protocol.append("genome parts: " + seqs.size() + "\t" + x +"\n" );/**/

		//selected transcripts
		if( selectedFile != null ) {
			selected = Tools.getSelection( selectedFile, maxSize, protocol );
			protocol.append("selected: " + selected.size() + (selected.size()<100? "\t"+ selected : "")+"\n");
		}
		
		//introns
		if( introns != null ) {
			ArrayList<String> fName = new ArrayList<String>();
			for( int i = 0; i < introns.getNumberOfParameters(); i++ ) {
				Parameter y = ((ParameterSet)introns.getParameterAt(i).getValue()).getParameterAt(0);
				if( y.isSet() ) {
					fName.add(y.getValue().toString());
				}
			}
			if( fName.size()>0 ) {
				HashMap<String, int[][][]>[] res = readIntrons( reads, protocol, verbose, seqs, null, fName.toArray(new String[fName.size()]) );
				check( protocol, res[0], "introns" );
				donorSites = res[0];
				acceptorSites = res[1];
			} else {
				acceptorSites = donorSites = null;
			}
		}
		
		//coverage
		if( cov != null ) {
			coverage = readCoverage(cov, protocol, verbose,true);
		}
	}
	
	//checks whether the RNA-seq data matches the reference genome
	static void check( Protocol p, HashMap<String,?> h, String text ) {
		Iterator<String> it = seqs.keySet().iterator();
		int anz = 0, found = 0;
		while( it.hasNext() ) {
			anz++;
			String chr = it.next();
			if( h.containsKey(chr) ) {
				found++;
			}
		}
		double d = (double) found / (double) anz; 
		if( d < 0.5 ) {//Warning if less than 50% of the reference sequences have been covered 
			p.appendWarning("Check RNA-seq data ("+text+"): " +(Math.round(d*100)) + "% of the sequences in the reference genome are covered.\n" );
		}
	}
	
	//user?
	private double contigThreshold;//threshold for initial solutions (filtering of contigs)
	private double regionThreshold;//threshold for regions
	private double hitThreshold;//threshold for hits
	
	private int MAX_INTRON_LENGTH;//the maximal intron length used in the DP to determine possible connections
	private int INTRON_GAIN_LOSS;
	
	private boolean approx;
	private boolean avoidPrematureStop; // avoid premature STOP codons in gene models?
	private double eValue;// the e-value for filtering search hits
	private int predictions;//how many predictions should be made
	
	private static final int MIN_INTRON_LENGTH = 30;//the minimal intron length used in the DP to determine possible connections between hits of the same query exon
	private static final int MAX_GAP = 5;
	private static final int missingAA = 10;//the minimal number of missing AA used for a self-loop in the DP 
	private static final int ignoreAAForSpliceSite = 30;// the maximal number of AA in a (blast) hit to be ignore
		
	//in
	private boolean proteinAlignment;
	private HashMap<String, String> cds;
	private HashMap<String,HashMap<String,Info>> transcriptInfo;

	//out
	private String prefix, tag;
	
	//alignment
	private DiscreteAlphabet aaAlphabet;
	private AlphabetContainer alph;
	private double[][] matrix;
	private int gapOpening;
	private int gapExtension;
	private AffineCosts cost;
	private Tools.Ambiguity ambiguity = Ambiguity.AMBIGUOUS;

	//splicing
	private boolean sp;
	//canonical 
	public static final String[] DONOR = {"GT", "GC"};
	public static final String ACCEPTOR = "AG";
	private static int[] intronic = { DONOR[0].length(), ACCEPTOR.length()};
	
	//statistics
	private int noL;
	
	//for logging
	private boolean verbose;
		
	private int maxSize;
	private long timeOut, maxTimeOut;
	
	public String getHeading() {
		return "gene\ttranscript\t#parts\t#predicted hits\t#blast hits\t#strands\tbest sum score"
				+ "\t#candidate strands\t#predictions\t#alignments"
				+ "\tbest final score"
				+ "\tchromosome\tstrand\tstart\tend\t#predicted parts\tfirst part\tfirst aa\tlast parts\tlast aa\t#*\tintron gain\tintron loss\tbackup\tcut"
				+ (proteinAlignment ?"\tminimal score\tcurrent score\toptimal score\t%positive\tpid\tmaxGap\tref_length\tpred_length":"")
				+ (acceptorSites == null ? "" :"\tacceptor evidence")
				+ (donorSites == null ? "" : "\tdonor evidence")
				+ (donorSites == null ? "" : "\tintron evidence")
				+ (coverage==null ? "" : "\tpercent covered\tmin coverage")
				+ "\tsimilar";
	}
	
	public GeMoMa( int maxSize, long timeOut, long maxTimeOut ) {
		this.maxSize = maxSize;
		this.timeOut = timeOut;
		this.maxTimeOut = maxTimeOut;
		timeOutWarning = new StringBuffer();
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
		
		//checking for available updates
		System.setProperty("java.net.useSystemProxies", "true");
		String site="https://www.jstacs.de/index.php/GeMoMa";
		System.out.println( "Searching for the new GeMoMa updates ..." );
		boolean checked = false;
		try{
			URL url = new URL(site);
			URLConnection con = url.openConnection();
			con.setConnectTimeout(5000);
			con.setReadTimeout(5000);
			BufferedReader br = new BufferedReader( new InputStreamReader(con.getInputStream()) );
			String line = null;
			//int l=0;
			while( (line=br.readLine()) != null ) {
				if( line.startsWith("<h") && line.contains("Version history") ) {
					break;
				}
				//l++;
			}
			//System.out.println(l);
			if( line != null ) {
				line = br.readLine();
				int idx = line.lastIndexOf("GeMoMa");
				if( idx >= 0 ) {
					idx+=6;
					int idx2 = line.indexOf("<",idx);
					if( idx2 >= 0 ) {
						String ver = line.substring(idx, idx2).trim();
						if( !ver.equals(GeMoMa.VERSION) ) {
							if( ver.compareTo(GeMoMa.VERSION)>0 ) {
								System.out.println("You are using GeMoMa " + GeMoMa.VERSION + ", but the latest version is " + ver + ".\nYou can download the latest version from " + site );
							} else {
								System.out.println("You are using an unofficial GeMoMa " + GeMoMa.VERSION + "; the latest official version is " + ver + "." );
							}
						} else {
							System.out.println("You are using the latest GeMoMa version.");
						}
						checked=true;
					}
				}
			}
		} catch( Exception ex ) {
			//ex.printStackTrace();
		}
		if( !checked ) {
			System.out.println("Could not connect to the GeMoMa-Homepage. Please check manually at " + site);
		}
		System.out.println();
		
		File basic = new File(Tools.GeMoMa_TEMP);
		basic.mkdirs();
		
		//reading parameters
		int maxSize = -1;
		long timeOut=3600, maxTimeOut=60*60*24*7;
		File ini = new File( jarfile.getParentFile().getAbsolutePath() + File.separator + "GeMoMa.ini.xml" );
		if( ini.exists() ) {
			//read
			StringBuffer xml = FileManager.readFile(ini);
			maxSize = XMLParser.extractObjectForTags(xml, "maxSize", Integer.class);
			timeOut = XMLParser.extractObjectForTags(xml, "timeOut", Long.class);
			maxTimeOut = XMLParser.extractObjectForTags(xml, "maxTimeOut", Long.class);
		} else {
			//default an write
			StringBuffer xml = new StringBuffer();
			XMLParser.appendObjectWithTags( xml, maxSize, "maxSize" );
			xml.append("\n");
			XMLParser.appendObjectWithTags( xml, timeOut, "timeOut" );
			xml.append("\n");
			XMLParser.appendObjectWithTags( xml, maxTimeOut, "maxTimeOut" );
			FileManager.writeFile(ini, xml);
		}
		//System.out.println(maxSize + "\t" + timeOut + "\t" + maxTimeOut );
		
		JstacsTool[] tools = {
				new GeMoMaPipeline(),
				
				new ExtractRNAseqEvidence(),
				new CheckIntrons(),
				new DenoiseIntrons(),
				
				new NCBIReferenceRetriever(),
				
				new Extractor(maxSize),
				new GeMoMa(maxSize, timeOut, maxTimeOut),
				new GeMoMaAnnotationFilter(),
				new AnnotationFinalizer(),
				new AnnotationEvidence(),
				new CompareTranscripts(),
				new SyntenyChecker(),
				new AddAttribute(),
				new GAFComparison(),
				new BUSCORecomputer()
				/*,
				new TranscribedCluster()/**/
		};
		boolean[] configureThreads = new boolean[tools.length];
		Arrays.fill(configureThreads,false);
		configureThreads[0]=true;
		
		//running the program
		if( args.length == 0 ) {
			System.out.println( "If you start with the tool with \"CLI\" as first parameter you can use the command line interface, otherwise you can use the Galaxy interface.");
		} else {
			if( args[0].equalsIgnoreCase("CLI") || args[0].equalsIgnoreCase("wiki") ) {
				CLI cli = new CLI( "This jar allows to run all parts of GeneModelMapper (GeMoMa) except the external search algorithm (e.g. tblastn).\n" 
					+ MORE.replaceAll("\\*\\*.*\\*\\*\n", "").replaceAll("\\*", "")
					+ "\n\nIf you use this tool, please cite\n\n" + REF[0] + "\n" + REF[1]
					,
					"CLI", configureThreads, tools );
				
				if( args[0].equalsIgnoreCase("CLI") ) {
					String[] part = new String[args.length-1];
					System.arraycopy(args, 1, part, 0, part.length);
					cli.run(part);
				} else {
					System.out.println("Creating help for wiki");
					StringBuffer sb = cli.wikiPage("GeMoMa-"+GeMoMaModule.VERSION+".jar CLI");
					//TODO check "in a nutshell" if parameter keys changed
					sb.insert(0, "== In a nutshell ==\n\n"
							+ "GeMoMa is a modular, homology-based gene prediction program with huge flexibility. However, we also provide a pipeline allowing to use GeMoMa easily. If you like to start GeMoMa for the first time, we recommend to use the GeMoMaPipeline like this\n"
							+ " java -jar GeMoMa-"+GeMoMaModule.VERSION+".jar CLI GeMoMaPipeline threads=<threads> outdir=<outdir> GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=<target_genome> i=<reference_1_id> a=<reference_1_annotation> g=<reference_1_genome>\n" + 
							"there are several parameters that need to be set indicated with '''&lt;'''foo'''&gt;'''. You can specify\n" + 
							"* the number of threads\n" + 
							"* the output directory\n" + 
							"* the target genome\n" + 
							"* and the reference ID (optional), annotation and genome. If you have several references just repeat <code>s=own</code> and the parameter tags <code>i</code>, <code>a</code>, <code>g</code> with the corresponding values.\n" +
							"In addition, we recommend to set several parameters:\n" + 
							"* <code>GeMoMa.Score=ReAlign</code>: states that the score from mmseqs should be recomputed as mmseqs uses an approximation\n" + 
							"* <code>AnnotationFinalizer.r=NO</code>: do not rename genes and transcripts\n" + 
							"* <code>o=true</code>: output individual predictions for each reference as a separate file allowing to rerun the combination step ('''GAF''') very easily and quickly\n" + 
							"If you like to specify the maximum intron length please consider the parameters <code>GeMoMa.m</code> and <code>GeMoMa.sil</code>.\n" 
							+ "If you have RNA-seq data either from own experiments or publicly available data sets (cf. [https://www.ncbi.nlm.nih.gov/sra NCBI SRA], [https://www.ebi.ac.uk/ena EMBL-EBI ENA]), we recommend to use them. You need to map the data against the target genome with your favorite read mapper. In addition, we recommend to check the parameters of the section '''DenoiseIntrons'''.\n"
							+ "\n");
					FileManager.overwriteFile( "gemoma-wiki.page",
						sb.toString()
						.replaceAll(MORE, "")
						.replace(Extractor.EXAMPLE, "")
						.replace("c=<cds_parts>", "c=<cds_parts> a=<assignment>")	
					);
					//cli.wiki("wiki");
				}
			} else {
				Galaxy galaxy = new Galaxy("", configureThreads, true, tools );
				galaxy.run(args);
			}
		}
		tools[0].clear();
	}
	
	public static HashMap<String, int[][]> simplify( HashMap<String, ArrayList<int[]>> current ) {
		Iterator<String> it = current.keySet().iterator();
		HashMap<String, int[][]> res = new HashMap<String, int[][]>();
		ArrayList<int[]> list = new ArrayList<int[]>();
		while( it.hasNext() ) {
			String key = it.next();
			int[][] val = current.get(key).toArray(new int[0][]);
			
			//sort
			Arrays.sort(val,GeMoMa.IntArrayComparator.comparator[0]);

			//iterate
			int i = 0;
			list.clear();
			while( i < val.length ) {
				//find continuous stretch
				int j = i+1;
				int max = val[i][1];
				
				while( j < val.length && max>=val[j][0] ) {
					if( max < val[j][1] ) {
						max = val[j][1];
					}
					j++;
				}
				
				if( j-i == 1 ) {
					//use as given
					list.add( val[i] );
				} else {
					/*
					System.out.println(key + " =====================================");
					for( int k=i; k <j; k++ ) {
						System.out.println(k + "\t" + Arrays.toString(val[k]) );
					}
					System.out.println();
					/**/
					//combine
					int second=0;
					do {
						//get minimal value
						int min = max;
						for( int k = i; k<j; k++ ) {
							if( val[k][2] > 0 && val[k][0] < min ) {
								min=val[k][0];
							}
						}
						//get smallest value larger than min
						second = max;
						for( int k = i; k<j; k++ ) {
							if( val[k][2] > 0 && val[k][0] > min && val[k][0] < second ) {
								second=val[k][0]-1;
							}
							if( val[k][2] > 0 && val[k][1] < second ) {
								second=val[k][1];
							}
						}
						//compute combined coverage
						int s = 0;
						for( int k = i; k<j; k++ ) {
							if( min == val[k][0] && val[k][1] >= second ) {
								s += val[k][2];
							}
						}
						list.add( new int[]{min,second,s} );
						
						//System.out.println(key + "\t" + min + "\t" + second + "\t" + s );
						for( int k = i; k<j; k++ ) {
							if( val[k][0]==min && val[k][1]>=second ) {
								val[k][0]=second+1;
								if( val[k][0]>val[k][1] ) {
									val[k][2] = 0;
								}
							}
						}
					} while( second<max );
				}				
				i = j;
			}
			
			//set
			res.put(key, list.toArray(new int[list.size()][]));
		}
		return res;
	}
	
	public static HashMap<String, int[][]> combine( boolean clone, HashMap<String, ArrayList<int[]>> basis, HashMap<String, ArrayList<int[]>> add ) {
		Iterator<String> it = add.keySet().iterator();
		//System.out.println("~~~~~~~~~~");
		//System.out.println(basis.size() + "\t" + add.size());
		while( it.hasNext() ) {
			String key = it.next();
			ArrayList<int[]> a = add.get(key);
			ArrayList<int[]> b = basis.get(key);
			//System.out.println(key + "\t" +a.size() + "\t" + (b==null?null:b.size()) );
			if( b==null ) {
				b = new ArrayList<int[]>();
				basis.put(key, b);
			}
			//add
			if( clone ) {
				for( int i = 0; i < a.size(); i++ ) {
					b.add(a.get(i).clone());
				}
			} else {
				b.addAll(a);
			}
			//System.out.println("=>" + b.size() );
		}
		return simplify( basis );
	}
	
	
	/**
	 * This method reads the coverage form user-specified files.
	 * 
	 * @param eps the user-specified data
	 * @param protocol the protocol for reporting
	 * 
	 * @return an array of {@link HashMap} where the first entry is the coverage of the forward strand and the second entry is coverage of the reverse strand.
	 * 		The coverage for each strand is given as {@link HashMap} with contig/chromosome names as key and <code>int[][]</code> arrays as entries.
	 * 		Each <code>int[]</code> array represents a tuple of start position, end position, and number of reads.
	 * 		The <code>int[][]</code> array is sorted according to the first entry (genomic order).  
	 * 
	 * @throws IOException
	 */
	public static HashMap<String, int[][]>[] readCoverage( ExpandableParameterSet eps, Protocol protocol, boolean verbose, boolean check ) throws IOException {
		HashMap<String, ArrayList<int[]>>[] initialCoverage = new HashMap[3];
		for( int i = 0; i < initialCoverage.length; i++ ) {
			initialCoverage[i] = new HashMap<String, ArrayList<int[]>>();
		}
		int anz = 0, u = 0;
		for( int i = 0; i < eps.getNumberOfParameters(); i++ ) {
			anz++;
			SimpleParameterSet sps = (SimpleParameterSet) eps.getParameterAt(i).getValue();
			sps = (SimpleParameterSet) sps.getParameterAt(0).getValue();
			if( sps.getNumberOfParameters() == 2 ) {
				readCoverage(initialCoverage[0],sps.getParameterForName("coverage_forward").getValue().toString(), protocol, verbose);
				readCoverage(initialCoverage[1],sps.getParameterForName("coverage_reverse").getValue().toString(), protocol, verbose);
			} else if( sps.getNumberOfParameters() == 1 ) {
				readCoverage(initialCoverage[2], sps.getParameterForName("coverage_unstranded").getValue().toString(), protocol, verbose);
				u++;
			} else {
				anz--;
			}
		}
		HashMap<String, int[][]>[] coverage;
		if( anz > 0 ) {
			coverage = new HashMap[2];
			coverage[0] = combine(true, initialCoverage[0], initialCoverage[2] );
			if( anz == u ) {
				coverage[1] = coverage[0];
				if( check ) check( protocol, coverage[0], "coverage" );
			} else {
				coverage[1] = combine(false, initialCoverage[1], initialCoverage[2] );
				if( check ) {
					check( protocol, coverage[0], "forward coverage" );
					check( protocol, coverage[1], "reverse coverage" );
				}
			}
		} else {
			coverage = null;
		}
		return coverage;
	}
	
	public static void readCoverage( HashMap<String, ArrayList<int[]>> coverage, String fName, Protocol protocol, boolean verbose ) throws IOException {
		//protocol.append("read coverage file: " + p.getName() + " " + new Date() + "\n");
		BufferedReader r = new BufferedReader( new FileReader( fName ) );
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
					chr = split[0];
					list = coverage.get(chr);
					if( list==null ) {
						list = new ArrayList<int[]>();
						coverage.put( chr, list );
					}
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
				list.add( new int[]{Integer.parseInt(split[1]), Integer.parseInt(split[2])-1, reads} );
			}
			i++;
		}
		r.close();
	}
	
	private static int[] count = new int[3];
	
	/**
	 * 
	 * @param threshold minimal number of reads to pass the filter
	 * @param protocol the protocol for reporting
	 * @param verbose whether to output additional information
	 * @param seqs the genome sequence
	 * @param diNucl a {@link HashMap} that might be filled with the dinucleotides at the splices sites or <code>null</code> 
	 * @param intronGFF the input files
	 * 
	 * @return  an array of {@link HashMap} where the first and second entry are the introns sorted according to the donor and acceptor site, respectively.
	 * 		The three-dimensional array encodes the information:
	 * 		<ul>
	 * 			<li>first dimension: strand: 0=forward, 1=reverse complement</li>
	 *			<li>second dimension: type: 0=donor, 1=acceptor, 2=#split reads</li>
	 * 			<li>third dimension: *counter*</li>
	 * 		</ul>
	 * 		Hence, introns[0].get(chr)[s][1][i] contains the information for the i-th acceptor site of strand s from chromosome chr, when sorted according to the donor sites.
	 * 		<b>If parameter diNucl is not equal to <code>null</code>, no introns will be returned.</b>
	 * 
	 * @throws IOException if at least one {@link File} cannot be read correctly 
	 */
	public static HashMap<String, int[][][]>[] readIntrons( int threshold, Protocol protocol, boolean verbose, HashMap<String, String> seqs, HashMap<String,int[]> diNucl, String... intronGFF ) throws IOException {
		HashMap<String, ArrayList<int[]>[]> spliceHash = new HashMap<String, ArrayList<int[]>[]>();
		ArrayList<int[]>[] h;
		
		String[] donor = null;
		String acceptor = null;
		String line;
		BufferedReader r;
		int num = 0;
		Arrays.fill( count, 0);
		for( String file : intronGFF ) {
			//System.out.println(file);
			r = new BufferedReader( new FileReader( file ) );
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
					
					if( c=='.' || diNucl != null ) {
						//try to identify on which strand the intron is located
						if( donor == null ) {
							donor = new String[]{Tools.rc(DONOR[0]),Tools.rc(DONOR[1])};
							acceptor = Tools.rc(ACCEPTOR);
						}
						String s = seqs.get(split[0]);
						if( s == null ) {
							r.close();
							throw new IllegalArgumentException("Did not find sequence " + split[0] + ", which should contain an intron");
						}
						String x = s.substring(a-1,a+1).toUpperCase();
						String y = s.substring(b-3,b-1).toUpperCase();
						if( c=='.' ) {
							boolean fwd =  (x.equals(DONOR[0]) || x.equals(DONOR[1])) && y.equals(ACCEPTOR);
							boolean bwd =  (y.equals(donor[0]) || y.equals(donor[1])) && x.equals(acceptor);
							if( fwd && !bwd ) {
								c='+';
							} else if( bwd && !fwd ) {
								c='-';
							}
						}
						
						if( diNucl != null ) {
							if( reads>= threshold ) {
								if( c=='-' ) {
									String hh=Tools.rc(x);
									x=Tools.rc(y);
									y=hh;
								}
								
								int[] v = diNucl.get(x+"-"+y);
								if( v == null ) {
									v = new int[2];
									diNucl.put(x+"-"+y, v);
								}
								v[0]++;
								v[1]+=reads;
							}
							num++;
						}
						//System.out.println(fwd + "\t" + bwd + "\t"+x + " .. " + y);
					}

					if( diNucl == null ) {
						if( c =='+' || c=='.' ) {
							h[0].add( new int[]{a, b, reads} );
						}
						if( c =='-' || c=='.' ) {
							h[1].add( new int[]{b, a, reads} );
						}
					}

					//stats
					switch( c ) {
						case '+': count[0]++; break;
						case '-': count[1]++; break;
						case '.': count[2]++; break;
					}
				}
			}
			r.close();
		}
		
		//reformat
		HashMap<String, int[][][]> donorSites = new HashMap<String, int[][][]>();
		HashMap<String, int[][][]> acceptorSites = new HashMap<String, int[][][]>();
		
		Entry<String, ArrayList<int[]>[]> e;
		int[][][] help; //strand / type(+,-,.) / info(start,end,reads)
		int[] site;
		Iterator<Entry<String, ArrayList<int[]>[]>> it = spliceHash.entrySet().iterator();
		while( it.hasNext() ) {
			e = it.next();
			h = e.getValue();
			help = new int[2][][];
			for( int k = 0; k < 2; k++ ) {
				help[k] = new int[h[k].size()][3];
				for( int m = 0; m<h[k].size(); m++ ) {
					site = h[k].get(m);
					help[k][m][0] = site[0];
					help[k][m][1] = site[1];
					help[k][m][2] = site[2];
				}
			}
			//new version 1.3.3
			int[][][] x = combineIntrons(help, 1,threshold);
			num += x[0][0].length + x[1][0].length;
			acceptorSites.put(e.getKey(), x);
			donorSites.put(e.getKey(), combineIntrons(help,0,threshold));
		}
		
		protocol.append("possible introns from RNA-seq (split reads>="+threshold+"): " + num + "\n");
		protocol.append("+: " + count[0] + "\n");
		protocol.append("-: " + count[1] + "\n");
		protocol.append(".: " + count[2] + "\n");
		
		return new HashMap[]{donorSites, acceptorSites};
	}
	
	private static int[][][] combineIntrons( int[][][] help, int idx, int threshold ) {
		LinkedList<int[]>[] list = new LinkedList[2];
		int[][][] vals = new int[2][][];
		for( int k = 0; k < 2; k++ ) {
			Arrays.sort(help[k], IntArrayComparator.comparator[idx]);
			list[k] = new LinkedList<int[]>();
			int[] last=null;
			//combine
			for( int m = 0; m<help[k].length; m++ ) {
				if( last!=null && help[k][m][0] == last[0] && help[k][m][1] == last[1] ) {
					//same intron => add number of split reads
					last[2] += help[k][m][2];
				} else {
					//new (=different) intron
					
					//number of split reads of last intron sufficient?
					if( list[k].size()>0 && list[k].getLast()[2]<threshold ) {
						list[k].removeLast();
					}
					
					//add new
					last=help[k][m].clone();
					list[k].add(last);
				}
			}
			//number of split reads of last intron sufficient?
			if( list[k].size()>0 && list[k].getLast()[2]<threshold ) {
				list[k].removeLast();
			}
			
			vals[k]=new int[3][list[k].size()];
			for( int m = 0; m<list[k].size(); m++ ) {
				last = list[k].get(m);
				vals[k][0][m] = last[0];
				vals[k][1][m] = last[1];
				vals[k][2][m] = last[2];
			}
		}
		return vals;
	}
	
	@Override
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp ) throws Exception {
		timeOutWarning.delete(0, timeOutWarning.length());
		verbose = (Boolean) parameters.getParameterForName("verbose").getValue();
		if( seqs == null ) {
			fill( protocol, verbose, maxSize,
					(String) parameters.getParameterForName("target genome").getValue(), 
					(String) parameters.getParameterForName("selected").getValue(),
					(Integer) parameters.getParameterForName("reads").getValue(),
					(ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterForName("introns")).getValue(),
					(ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterForName("coverage")).getValue()
			);			
		} else {
			protocol.append("Using pre-set values.\n");
		}
		
		BufferedReader r = null;
		String line;

		String assignment = (String) parameters.getParameterForName("assignment").getValue();
		String searchFile = (String) parameters.getParameterForName("search results").getValue();
		File search;
		if( (Boolean) parameters.getParameterForName("sort").getValue() ) {
			search = Tools.externalSort(searchFile, 500000, 1, protocol,assignment!=null)[0];
		} else {
			search = new File( searchFile );
		}	
		
		MAX_INTRON_LENGTH = (Integer) parameters.getParameterForName( "maximum intron length" ).getValue();
		boolean staticIntronLength = (Boolean) parameters.getParameterForName( "static intron length" ).getValue();
		eValue = (Double) parameters.getParameterForName( "e-value" ).getValue();
		contigThreshold = (Double) parameters.getParameterForName( "contig threshold" ).getValue();
		regionThreshold = (Double) parameters.getParameterForName( "region threshold" ).getValue();
		hitThreshold = (Double) parameters.getParameterForName( "hit threshold" ).getValue();
		
		predictions = (Integer) parameters.getParameterForName( "predictions" ).getValue();
		gapOpening = (Integer) parameters.getParameterForName("gap opening").getValue();
		gapExtension = (Integer) parameters.getParameterForName("gap extension").getValue();
		INTRON_GAIN_LOSS = (Integer) parameters.getParameterForName("intron-loss-gain-penalty").getValue();
		avoidPrematureStop = (Boolean) parameters.getParameterForName("avoid stop").getValue();
		approx = (Boolean) parameters.getParameterForName("approx").getValue();		
		prefix = (String) parameters.getParameterForName("prefix").getValue();
		tag = (String) parameters.getParameterForName("tag").getValue();
		timeOut = (Long) parameters.getParameterForName( "timeout" ).getValue();
		
		proteinAlignment = (Boolean) parameters.getParameterForName("protein alignment").getValue();
		
		Score score = (Score) parameters.getParameterForName("Score").getValue();
		
		String[] split;
		//read substitution matrix
		r = new BufferedReader( new InputStreamReader( Tools.getInputStream(parameters.getParameterForName("substitution matrix"), "projects/gemoma/test_data/BLOSUM62.txt" ) ) );
		while( (line=r.readLine()) != null && line.charAt(0)=='#' );
		String[] abc = line.split("\\s+");
		String[] abc2 = new String[abc.length-1];
		System.arraycopy( abc, 1, abc2, 0, abc2.length );
		int i = 0;
		while( i < abc2.length && !abc2[i].equals("*") ) {
			i++;
		}
		boolean addStop = i == abc2.length;
		if( addStop ) {
			String[] h = new String[abc2.length+1];
			System.arraycopy( abc2, 0, h, 0, abc2.length );
			h[abc2.length]="*";
			abc2=h;
		}	

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
		if( addStop ) {
			for( i = 0; i < abc2.length-1; i++ ) {
				matrix[i][abc2.length-1] = 4;
			}
			Arrays.fill( matrix[abc2.length-1], 4 );
			matrix[abc2.length-1][abc2.length-1] = -1;
			protocol.append("Expand substitution matrix: Use default costs for *\n");
		}
			
		cost = new AffineCosts(gapOpening, new MatrixCosts(matrix, gapExtension));
			//new MatrixCosts(matrix, 5);

		//read assignment
		if( assignment != null ) {
			transcriptInfo = new HashMap<String,HashMap<String,Info>>();
			HashMap<String,Info> c;
			r = new BufferedReader( new FileReader( assignment ) );
			while( (line=r.readLine())!= null ) {
				if( line.length()>0 && line.charAt(0) != '#' ) {
					split = line.split("\t");
					c = transcriptInfo.get( split[0] );
					if( c == null ) {
						c = new HashMap<String, Info>();
						transcriptInfo.put(split[0], c);
					}
					c.put(split[1], new Info( split[2], split[3], split.length>=12 ? split[11] : "", split[9], aaAlphabet ) );
				}
			}
			r.close();
			progress.setLast(transcriptInfo.size());
			progress.setCurrent(0);
		} else {
			progress.setIndeterminate();
			transcriptInfo = null;
		}

		cds = Tools.getFasta((String) parameters.getParameterForName("cds parts").getValue(),15000);
		protocol.append("CDS: " + cds.size() + " / " + (transcriptInfo==null?cds.size():transcriptInfo.size()) + "\n" );// + "\t" + Arrays.toString(cdsInfo.keySet().toArray()) );
				
		//read genetic code
		code = Tools.getCode(Tools.getInputStream(parameters.getParameterForName("genetic code"), "projects/gemoma/test_data/genetic_code.txt" ));
		code.put("NNN", 'X');//add for splitting at NNN (in align method)
		
		if( selected != null ) {
			progress.setLast(selected.size());
		}
		
		if( donorSites != null ) {
			sp = (Boolean) parameters.getParameterForName("splice").getValue();
		}
		
		//read blast output and compute result
		boolean okay = true;
		ArrayList<TextResult> res = new ArrayList<TextResult>();
		Exception e = null;
		
		TranscriptPredictor tp = new TranscriptPredictor(true, parameters, temp);
		String problem=null, old=null;
		try {			
			//collect blast hits per transcript, split for chromosome, strand and cds part
			r = new BufferedReader( new FileReader( search ) );
			protocol.append(getHeading()+"\n");
			noL=0;
			HashMap<String, HashMap<Integer,ArrayList<Hit>>[]> hash = new HashMap<String, HashMap<Integer,ArrayList<Hit>>[]>();
			HashSet<String> used = new HashSet<String>();
			while( (line=r.readLine()) != null ) {
				if( line.length() > 0 ) {
					if( old == null || !line.startsWith(old) ) {
						if( old != null ) {
							tp.numberOfLines=noL;
							progress.add( tp.compute(old, hash, staticIntronLength) );
							protocol.append( tp.protocol.toString() );
							noL=0;
							//clear
							hash.clear();
						}
						old = line.substring(0, line.indexOf('\t'));
						if( transcriptInfo != null ) {
							old = old.substring(0, old.lastIndexOf('_')+1);
						} else {
							old += "\t";
						}
						if( used.contains(old) ) {
							throw new Exception("The search results seem to be unsorted. ID " + old + " has been seen before.");
						} else {
							used.add(old);
						}
					}
					
					//parse blast hit
					addHit(hash, line, score, tp);
				}
			}
			if( old != null ) {
				tp.numberOfLines=noL;
				progress.add( tp.compute(old, hash, staticIntronLength) );
				protocol.append( tp.protocol.toString() );
			}			
		} catch ( Throwable er ) {
			protocol.appendThrowable(er);
			if( er instanceof Exception ) {
				e = (Exception) er;
			} else {
				e = new Exception( "Forwarding " + er.getClass().getName() + ": " + er.getMessage());
				e.setStackTrace( er.getStackTrace() );
			}
			okay=false;
			problem=old;
		} finally {					
			//close output;
			tp.gff.close();
			
			if( timeOutWarning.length()>0 ) protocol.append( "\ntime-out warning: " + timeOutWarning );
			
			if( okay ) {
				String gff = tp.gffFile.getAbsolutePath();
				
				res.add( new TextResult(DEF_RES, "Result", new FileParameter.FileRepresentation(gff), "gff", getToolName(), null, true) );
			}
			
			if( r != null ) {
				r.close();
			}
			tp.close();
		}
		if( okay ) {
			return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
		} else {
			if( problem!=null && !(e instanceof InterruptedException)) System.err.println("\nProblem while gene: " + problem);
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
	public static class IntArrayComparator implements Comparator<int[]> {
		/**
		 * The default instances
		 */
		public static IntArrayComparator[] comparator = {
				new IntArrayComparator(0,1),
				new IntArrayComparator(1,0),
				new IntArrayComparator(0)
		};
		
		private int[] idx;
		
		private IntArrayComparator( int... idx ) {
			this.idx = idx;
		}
		
		public int compare( int[] o1, int[] o2) {
			int res=0, i = 0;
			while( i < idx.length && (res=Integer.compare(o1[idx[i]], o2[idx[i]])) == 0 ) {
				i++;
			}
			return res;
		}	
	}
	
	private double rescore( String a1, String a2, DiscreteAlphabet abc, double[][] matrix, double gapOpen, double gapExtend ) throws WrongAlphabetException, WrongLengthException {
		if( a1.length() != a2.length() ) {
			throw new WrongLengthException("aligned sequences need to have the same length");
		}
		int i = 0;
		double score = 0;
		do {
			if( a1.charAt(i) == '-' ) {
				int start = i;
				do {
					i++;
				} while( i < a1.length() && a1.charAt(i) == '-' );
				score += gapOpen + gapExtend*(i-start);
			} else if( a2.charAt(i) == '-' ) {
				int start = i;
				do {
					i++;
				} while( i < a2.length() && a2.charAt(i) == '-' );
				score += gapOpen + gapExtend*(i-start);
			} else {
				//(mis)match
				score += matrix[abc.getCode(a1.substring(i, i+1))][abc.getCode(a2.substring(i, i+1))];
				i++;
			}
		} while( i < a1.length() );
		return -score;
	}
	
	/**
	 * Parses one line of the tab-separated search results, filters by e-Value and creates if a {@link Hit} if the e-value is smaller than a threshold.
	 * 
	 * @param hash this {@link HashMap} contains the hits
	 * @param line a line of the blast output
	 * @throws WrongLengthException 
	 * @throws WrongAlphabetException 
	 */
	@SuppressWarnings("unchecked")
	void addHit( HashMap<String, HashMap<Integer,ArrayList<Hit>>[]> hash, String line, Score score, TranscriptPredictor tp ) throws WrongAlphabetException, WrongLengthException {
		String[] split = line.split("\t");
		
		if( !seqs.containsKey( split[1] ) ) {
			throw new IllegalArgumentException( "There is no sequence with sequence ID " + split[1] + " in the target genome." );
		}
		
		try {
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
			
			int s = Integer.parseInt(split[scoreIndex]);
			String queryAlign = split[20];
			String targetAlign = split[21];		
			switch( score ) {
				case Trust: break;
				case ReScore:
					s = (int) rescore(queryAlign, targetAlign, aaAlphabet, matrix, gapOpening, gapExtension);
					break;
				case ReAlign:
					Sequence q = Sequence.create( alph, queryAlign.replaceAll("-", "") );
					Sequence t = Sequence.create( alph, targetAlign.replaceAll("-", "") );
					
					StringAlignment sa = tp.align.getAlignment(AlignmentType.GLOBAL, q, t);
					
					queryAlign = sa.getAlignedString(0);
					targetAlign = sa.getAlignedString(1);
					s = -(int)sa.getCost();
					
					break;
				default:
					throw new IllegalArgumentException( score.name() );
			}
			

			//insert
			Hit h = new Hit( split[0], split[1], 
					Integer.parseInt(split[6]), Integer.parseInt(split[7]), Integer.parseInt(split[22]),  
					Integer.parseInt(split[8]), Integer.parseInt(split[9]), 
					s, queryAlign, targetAlign, "search algorithm;" );
			
			lines.add(h);
			//h.addSplitHits(lines);
			noL++;
		} catch( ArrayIndexOutOfBoundsException aiobe ) {
			ArrayIndexOutOfBoundsException moreDetails = new ArrayIndexOutOfBoundsException(
					"Please check the search algorithm input. It seems to have less columns than expected."
					+ "\nline: " + line
					+ "\nsplit: " + Arrays.toString(split)
					+ "\noriginal message: " + aiobe.getMessage()				
					+ "\n"
			);
			moreDetails.setStackTrace(aiobe.getStackTrace());
			throw moreDetails;
		}
	}
	
	
	
	/**
	 * Class representing search hits. Hence, the query is an amino acid sequence and the target is a genomic DNA sequence.
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
		}
		
		@SuppressWarnings("unchecked")
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
			return clone;
		}
		
		public String toString() {
			return queryID + "\t" + targetID + "\t" + (forward?"+":"-")
					+ "\t" + queryStart + "\t" + queryEnd +"\t"+queryLength + "\t" + targetStart + "\t" + targetEnd
					+ "\t" + score + "\t" + queryAlign + "\t" + targetAlign + "\t" + info;
		}
	
		/**
		 * This method splits the current hit at stop codons (*) and adds the (refined) parts into the given list.
		 * 
		 * @param add the list of valid {@link Hit}s
		 * 
		 * @throws WrongAlphabetException if the is a problem in the score computation
		 */
		public void addSplitHits(ArrayList<Hit> add) throws WrongAlphabetException {
			//int anz=0, start = add.size();
			
			int idx=-1;
			int old=0, qStart=queryStart, tStart = forward? targetStart : targetEnd, factor = forward ? 1 : -1, f, a, b;
			do {
				do {
					idx = targetAlign.indexOf('*', idx+1);
				} while( idx >=0 && queryAlign.charAt(idx)=='*' );
				
				a=b=0;
				if( idx < 0 ) {
					//there is no (further) premature stop codon that is not in the query sequence
					if( old == 0 ) {
						add.add(this);
						return;
					} else {
						f = targetAlign.length();
					}
				} else {
					//there is a premature stop codon that is not in the query sequence
					f = idx;
					if( queryAlign.charAt(idx)=='-' ) {
						a=1;
						b=0;
					} else {
						a=b=1;
					}
				}
					
				String t = targetAlign.substring(old, f);
				String q = queryAlign.substring(old, f);
				
				//shorten remaining alignment from both side if there are any gaps
				int off1=0, oq1=0;
				while( off1 < t.length() && (t.charAt(off1) =='-' || q.charAt(off1)=='-') ) {
					if( q.charAt(off1)=='-') {
						oq1++;
					}
					off1++;
				}
				
				if( off1 == t.length() ) {
					//the alignment does not contain a single match or mismatch (only gaps), hence, we don't need this part of the alignment
					int qg = getGapLength(q);
					int tg = getGapLength(t);
					qStart+=q.length()-qg+b;
					tStart+=factor*3*(t.length()-tg+a);
				} else {
					int off2=t.length()-1, oq2=0, n = 0;
					if( off1==t.length() ) {
						off2=0;
					} else {
						while( off2 >= 0 && (t.charAt(off2) =='-' || q.charAt(off2)=='-') ) {
							if( q.charAt(off2)=='-') {
								oq2++;
							}
							n++;
							off2--;
						}
						off2++;
					}
					
					q=q.substring(off1, off2);
					t=t.substring(off1, off2);
				
					int sc = getScore(q, t);
					int qg = getGapLength(q);
					int tg = getGapLength(t);
					
					int tS = tStart+factor*3*oq1;
					if( sc > 0 )  {
						Hit h =  new Hit(queryID, targetID, qStart+(off1-oq1), qStart+(off1-oq1)+q.length()-qg-1, queryLength, tS, tS +factor*(3*(t.length()-tg)-1), sc, q, t, info + "split by '*';" );
						add.add( h );
						//anz++;
					}
					qStart+=(off1-oq1) + q.length()-qg + (n-oq2) + b;
					tStart+=factor*3*(oq1 + t.length()-tg + oq2 + a);
				}
				old=f+1;
			} while( idx>=0 );
			/*
			if( anz>1 ) {
				System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
				System.out.println(this);
				while( start < add.size() ) {
					System.out.println(add.get(start++));
				}
			}*/
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
		public void prepareSpliceCandidates( String chr, MyAlignment align, int maxAllowedIntron ) throws CloneNotSupportedException, WrongAlphabetException {
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
				Tools.getMaximalExtension(chr, forward, false, forward?-1:1, forward ? targetStart+add-1 : targetEnd-add, '*','*', seqUp, up, a.substring(0,add/3), eAcc, intronic[1]-1, maxAllowedIntron, code );
				int s = seqUp.length()-add-eAcc;
				if( forward ) {
					maxUpstreamStart = targetStart - s+intronic[1];
				} else {
					maxUpstreamStart = targetEnd + s-intronic[1];
				}
				up=up.delete(up.length()-(add/3), up.length());
			
				int anz = getExperimentalCandidates(acceptorSites, 1, accCand, add, false, maxUpstreamStart);
				ae = anz > 0;

				//fall back
				if( acceptorSites == null || (sp && anz == 0) ) {
					getCanonincalCandidates( seqUp, accCand, ACCEPTOR, true, true, s );
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
							accCandScore[p].add( (int) -align.getCost(qPart.getLength(), pos)-score );//XXX improve splice site scoring
						}
					}
				}
				
				//donor
				down = new StringBuffer();
				Tools.getMaximalExtension(chr, forward, true, forward?1:-1, forward ? targetEnd-add : (targetStart+add-1), '*','*', seqDown, down, a.substring(a.length() - add/3), eDon, intronic[0]-1, maxAllowedIntron, code );
				s = add+eDon;
				if( forward ) {
					maxDownstreamEnd = targetEnd + seqDown.length()-s -intronic[0];
				} else {
					maxDownstreamEnd = targetStart - seqDown.length()+s +intronic[0];
				}
				down = down.delete(0,add/3);
				
				int m=0;
				for( m = 0; m < donCand.length; m++ ) {
					donCand[m] = ArrayHandler.createArrayOf(new IntList(), 3);
					donCandScore[m] = ArrayHandler.createArrayOf(new IntList(), 3);
				}
				
				anz = getExperimentalCandidates(donorSites, 0, donCand[0], add, true, maxDownstreamEnd);
				de = anz > 0;

				//fall back
				if( donorSites == null || (sp && anz == 0) )
				{
					for( m = 0; m < DONOR.length; m++ ) {
						getCanonincalCandidates( seqDown, donCand[m], DONOR[m], false, false, s );
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
								donCandScore[m][p].add( (int) -align.getCost(qPart.getLength(), pos) -score );//XXX improve splice site scoring
							}
						}
					}
				}
				splice = true;
			}
		}
			
		private int getExperimentalCandidates( HashMap<String, int[][][]> sites, int sIndex, IntList[] cand, int add, boolean invert, int max ) {
			int anz = 0;
			if( sites != null ) {
				int[][][] vals = sites.get(this.targetID);
				if( vals != null ) {
					int z, start=-10, end=-10, ref=-10, f=0, idx, d, m, old;
					
					z = forward?0:1;
					for( int i = 0; i < 3; i++) {
						cand[i].clear();
					}
					
					boolean b = invert ? !forward : forward;

					if( b ) {
						start = max;
						end = this.targetStart + add;
						ref = targetStart;
						f = 1;
					} else {
						start = this.targetEnd - add;
						end = max+1;
						ref = targetEnd+1;
						f = -1;
					}

					idx = Arrays.binarySearch(vals[z][sIndex], start);
					if( idx < 0 ) {
						idx = -(idx+1);
					}
					old=-1;
					while( idx < vals[z][sIndex].length && vals[z][sIndex][idx] <= end ) {
						if( vals[z][sIndex][idx] != old ) {
							d = f*(ref-vals[z][sIndex][idx]);
							m = d%3;
							if( m < 0 ) {
								m+=3;
							}
							cand[m].add( d );
							anz++;
							old=vals[z][sIndex][idx];
						}
						idx++;
					}
				}
			}
			return anz;
		}
		
		private void getCanonincalCandidates( StringBuffer region, IntList[] sites, String p, boolean add, boolean reverse, int ref ) {
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
		
		public void setBorderConstraints( boolean firstExon, boolean lastExon, MyAlignment align ) throws IllegalArgumentException, WrongAlphabetException {
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
/*if( idx == Integer.MIN_VALUE ) {
								protocol.append("HELP\n");
							}*/
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
						firstAA=alph.getSymbol(0, t.discreteVal(0)).charAt(0);//.toString(0,1).charAt(0);
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

	
	private int getGapLength( String s ) {
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
	
	private static class Info {
		int[] exonID;
		int[] phase;
		String[] splitAA;
		int maxIntron;
		String protein;
		
		Info( String exonID, String phase, String splitAA, String maxIntron, DiscreteAlphabet aaAlphabet ) {
			String[] split = exonID.split( ",[ ]*" );
			this.exonID = new int[split.length];
			for( int i = 0; i < this.exonID.length; i++ ) {
				this.exonID[i] = Integer.parseInt(split[i].trim());
			}
			
			//XXX remove
			phase=null;			
			if( phase != null ) {
				split = phase.split( "," );
				this.phase = new int[split.length];
				for( int i = 0; i < this.phase.length; i++ ) {
					this.phase[i] = Integer.parseInt(split[i].trim());
				}
			}
			
			this.splitAA = new String[this.exonID.length];
			Arrays.fill(this.splitAA, "");
			if( splitAA != null ) {
				split = splitAA.split( "," );
				for( int i = 0; i < split.length; i++ ) {
					if( split[i].length()>1 ) {
						throw new IllegalArgumentException("Expect only one additonal amino acids per exon-exon boundary.");
					}
					if( split[i].length()==1 && !aaAlphabet.isSymbol(split[i]) ) {
						throw new IllegalArgumentException("Unknown character " + split[i] + " in additional amino acid column of assignment file.");						
					}
					this.splitAA[i] = split[i];
				}
			}
			
			if( maxIntron== null || maxIntron.equals("NA") ) {
				this.maxIntron = 0;
			} else {
				this.maxIntron = Integer.parseInt(maxIntron);
			}
		}
		
		private Info() {
			exonID = new int[1];
			phase = new int[0];
			splitAA = new String[0];
			maxIntron=0;
		}
		
		static Info DUMMY = new Info();
		
		public String toString() {
			return "parts: "+ Arrays.toString(exonID)+"\n"
					+ "phase: "+ Arrays.toString(phase)+"\n"
					+ "splitAA: "+ Arrays.toString(splitAA)+"\n";
		}
		
		public void prepareProtein( HashMap<String,String> cds, String prefix ) {
			StringBuffer sb = new StringBuffer();
			for( int i = 0; i<exonID.length; i++ ) {
				String cdsPart = cds.get(prefix+exonID[i]);
				if( cdsPart!= null ) {
					sb.append( cdsPart );
				}
				if( i < splitAA.length ) {
					sb.append(splitAA[i]);
				}
			}
			protein = sb.toString();
		}
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
		//private int[] parts, phase;
		private int[] revParts, oldSize;
		private int[] length, cumLength, revCumLength;
		private int[][] start, end, qStart, qEnd, qL;
		private int maxAllowedIntron;
		
		double[][] used;
		
		int[][][][] splice; //part 1 idx, anz, parts2 idx (relative), anz
		static final int NOT_COMPUTED = Integer.MAX_VALUE;
		static final int NO_SPLICE_VARIANT = Integer.MIN_VALUE;
				
		private int numberOfLines,
			numberOfTestedStrands,
			numberOfPairwiseAlignments;
		
		/**
		 * This field is used to store the possible extensions for a splice variant
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
		
		private final ScheduledThreadPoolExecutor executorService;
		
		private StringBuffer gffHelp;
		
		private Info currentInfo;
		
		private BufferedWriter gff;
		private File gffFile;
		private StringBuffer protocol;
		private MyAlignment align, align1, align2;
						
		public TranscriptPredictor( boolean add, ToolParameterSet parameters, String temp ) throws Exception {
			protocol = new StringBuffer();
			gffFile = Tools.createTempFile("gff", temp);
			
			gff = new BufferedWriter( new FileWriter( gffFile ) );
			if( add ) {
				gff.append("##gff-version 3");
				gff.newLine();
				gff.append(GeMoMa.INFO + getShortName() + " " + getToolVersion() + "; ");
				String info = JstacsTool.getSimpleParameterInfo(parameters);
				if( info != null ) {
					gff.append("SIMPLE PARAMETERS: " + info );
				}
				gff.newLine();
			}
			
			align = new MyAlignment( cost, this );
			align1 = new MyAlignment( cost, this );
			align2 = new MyAlignment( cost, this );
			
			result = new PriorityQueue<Solution>();
			
			executorService = //Executors.newSingleThreadExecutor();
				new ScheduledThreadPoolExecutor(1, (runnable) -> {
		            Thread thread = new Thread(runnable, this.getClass().getName());
		            thread.setDaemon(true);
		            return thread;
		        });
			executorService.setRemoveOnCancelPolicy(true);
			
			gffHelp = new StringBuffer();
		}
		
		public void close() {
			List<Runnable> list = executorService.shutdownNow();
			if( list.size()> 0 ) {
				for( Runnable r: list ) {
					System.out.println( "open Runnable: " + r.toString() );
				}
			}
		}

		/**
		 * This method computes for all transcripts of one gene model the predicted gene models in the target organism 
		 * 
		 * @param name gene name ending with an underscore
		 * @param hash the collected search {@link Hit}s
		 * 
		 * @return the number of reference transcripts used for prediction
		 * 
		 * @throws Exception if something went wrong
		 */
		@SuppressWarnings("unchecked")
		public int compute( String name, HashMap<String,HashMap<Integer,ArrayList<Hit>>[]> hash, boolean staticIntronLength ) throws Exception {
			protocol.delete(0, protocol.length());
			
			geneName = name.substring(0,name.length()-1);
			
			HashMap<String,Info> transcript;
			String[] s;
			if( transcriptInfo != null ) {
				transcript = transcriptInfo.get(geneName);
				if( transcript == null ) {
					return 0;
				}
				s = new String[transcript.size()];
				transcript.keySet().toArray(s);
				Arrays.sort(s);
			} else {
				transcript = null;
				s = new String[]{ geneName }; //removed .toUpperCase()
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
					int max = 0;
					if( transcript != null ) {
						currentInfo = transcript.get(s[i]);
						for( int p = 0; p < currentInfo.exonID.length; p++ ) {
							if( currentInfo.exonID[p] > max ) {
								max = currentInfo.exonID[p];
							}
						}
					} else {
						currentInfo=Info.DUMMY;
					}
					if( staticIntronLength ) {
						maxAllowedIntron = MAX_INTRON_LENGTH;
					} else {
						maxAllowedIntron = currentInfo.maxIntron+MAX_INTRON_LENGTH;
					}
					
					revParts = new int[max+1];
					if( transcript != null ) {
						Arrays.fill(revParts, -1);
						for( int p = 0; p < currentInfo.exonID.length; p++ ) {
							revParts[currentInfo.exonID[p]] = p;
						}
						
						length = new int[currentInfo.exonID.length];
						for( int j = 0; j < currentInfo.exonID.length; j++ ) {
							String st = geneName + (transcriptInfo==null ? "" : ("_" + currentInfo.exonID[j]));
							st = cds.get( st );
							length[j] = st!=null ? st.length() : 0;
								
						}
					} else {
						revParts[0] = 0;
						length=new int[]{cds.get(geneName).length()};
					}
					
					cumLength = new int[length.length];
					for( int j = 1; j < length.length; j++ ) {
						cumLength[j] = cumLength[j-1] + length[j-1];
					}
					
					revCumLength = new int[length.length];
					for( int j = length.length-2; j >= 0; j-- ) {
						revCumLength[j] = revCumLength[j+1] + length[j+1];
					}
					
					int l = currentInfo.exonID.length;
					sums = new int[l][];
					start = new int[l][];
					end = new int[l][];			
					qStart = new int[l][];
					qEnd = new int[l][];
					qL = new int[l][];
					
					if( transcriptInfo == null ) {
						computeTranscript( s[i], info, hash );
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
								for( int k = 0; k < currentInfo.exonID.length; k++ ) {
									hits = v[f].get(currentInfo.exonID[k]);
									if( hits != null ) {
										w[f].put( currentInfo.exonID[k], clone( hits ) );
										hasHits=true;
									}
								}
							}
							data.put(h[j], w);
						}
						if( hasHits ) {
							currentInfo.prepareProtein(cds, geneName +(transcriptInfo==null? "":"_"));
							computeTranscript( s[i], info, data );
							currentInfo.protein=null;
						}
					}
				}
			}
			
			return s.length;
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
		
		int[] res = null;
		
		String getRefProtein(String transcriptName) {
			String ref;
			if( currentInfo!=Info.DUMMY ) {
				ref=currentInfo.protein;
			} else {
				ref=cds.get(transcriptName);
			}
			return ref;
		}
		
		/**
		 * This method predicts gene model of a given transcript.
		 * 
		 * @param transcriptName the transcript name
		 * @param assign the assigned CDS parts
		 * @param hash the search {@link Hit}s
		 * 
		 * @throws Exception if something went wrong
		 */
		private void computeTranscript( String transcriptName, String info, HashMap<String,HashMap<Integer,ArrayList<Hit>>[]> hash ) throws Exception {
			if( verbose ) {
				protocol.append("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
				protocol.append( geneName + "\t" + transcriptName + "\t" + Arrays.toString(currentInfo.exonID) + "\t" + new Date() + (info==null?"":("\t"+info)) + "\n");
				protocol.append( Arrays.toString(length) + "\n");
				protocol.append( Arrays.toString(cumLength) + "\n");
				protocol.append( Arrays.toString(revCumLength) + "\n");
			}
			protocol.append( geneName+"\t"+transcriptName );
			//protocol.flush();

			numberOfTestedStrands=numberOfPairwiseAlignments=0;
			oldSize = new int[currentInfo.exonID.length];
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
			res=null;
			Future f = null;
			Thread t = null;
			try {
				t = new Thread() {
					 public void run() {
						 try {
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
		        			res = new int[]{ bestSumScore, k };
						 } catch( Exception e ) {
							 RuntimeException re = new RuntimeException(e.getMessage());
							 re.setStackTrace(e.getStackTrace());
							 throw re;
						 } finally {
							 synchronized ( protocol ) {
								protocol.notify();
							}
						 }
					 }
				};
				
				f = executorService.submit(t);
				f.get(timeOut,TimeUnit.SECONDS);
	        } catch (TimeoutException e) {
	        	protocol.append( "\t" +timeoutMsg + " of " + timeOut + " seconds\n");
	        	protocol.append(new Date() + "\t" + System.currentTimeMillis() +"\n");
	        	out=false;

	        	timeOutWarning.append( (timeOutWarning.length()>0?", ":"") + transcriptName );
	        	
	        	//stop thread
	        	t.interrupt(); //will probably not help
	        	
	        	synchronized( protocol ) {
	        		//set variable to force an Exception 
		        	align1.clear();
		        	align2.clear();
		        	align.clear();

		        	int[][] s = sums; 
		        	sums=null;
		        	
		        	//wait for an Exception
		        	protocol.wait();
		        	
		        	//replace
		        	align.destroyed=false;
		        	align1.destroyed=false;
		        	align2.destroyed=false;
		        	sums=s;
	        	}
	        	protocol.append(new Date() + "\t" + System.currentTimeMillis() +"\n" );
	        } finally {
	        	splice=null;
	        	used=null;
	        }
			
			if( out ) {
				String seq = proteinAlignment ? getRefProtein(transcriptName) : null;
				int counts = 0;
				for( int i = 0; i < predictions && result.size()>0; i++ ) {
					Solution best = result.poll();
					if( best.hits.size() > 0 ) {
						counts++;
						best = best.refine(transcriptName);
						
						best.writeSummary( geneName, transcriptName, i);
						gff.append( ";score=" + best.score);
						
						gffHelp.delete(0, gffHelp.length());
						int anz = best.writeGFF( transcriptName, i, gffHelp );
						
						if( transcriptInfo != null ) {
							gff.append( ";ce=" + anz + ";rce=" + currentInfo.exonID.length );
						}						
						if( acceptorSites != null ) {
							gff.append( ";tae=" + (best.A==0? "NA": decFormat.format(best.a/(double)best.A)) );
						}
						if( donorSites != null ) {
							gff.append( ";tde=" + (best.D==0? "NA": decFormat.format(best.d/(double)best.D)) );
							gff.append( ";tie=" + (best.I==0? "NA": decFormat.format(best.i/(double)best.I)) );
							gff.append( ";minSplitReads=" + (best.I==0? "NA": best.minSplitReads) );
						}
						if( coverage != null && coverage[best.forward?0:1] != null ) {
							gff.append( ";tpc=" + decFormat.format(best.Cov/(double)best.Len) );
							gff.append( ";minCov=" + best.minC );
							gff.append( ";avgCov=" + decFormat.format(best.Sum/(double)best.Len) );
						}
						
						int id = 0, pos = 0, gap=-1, currentGap=0, maxGap=0;
						String s0, s1 = null, pred = best.getProtein();
						PairwiseStringAlignment psa = null;
						if( seq != null ) {
							
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
							
							//additional GFF tags
							gff.append( ";pAA=" + decFormat.format(pos/(double)s1.length()) + ";iAA=" + decFormat.format(id/(double)s1.length()) );//+ ";maxGap=" + maxGap + ";alignF1=" + (2*aligned/(2*aligned+g1+g2)) ); 
						}
						int preMatureStops=0, j=-1;
						String aa =pred.substring(0, pred.length()-1);
						while( (j=aa.indexOf('*',j+1)) >= 0 ) {
							preMatureStops++;
						}
						gff.append( ";nps=" + preMatureStops );
						
						//short info
						protocol.append( ((i == 0 && !verbose)? "" : (geneName+"\t"+transcriptName)) + "\t" + currentInfo.exonID.length + "\t" + best.hits.size()
								+ "\t" + numberOfLines + "\t" + numberOfTestedStrands + "\t" + res[0]
								+ "\t" + res[1] + "\t" + (result.size()+1) + "\t" + numberOfPairwiseAlignments
								+ "\t" + best.score
								+ "\t" + best.hits.getFirst().targetID + "\t" + (best.forward?+1:-1)
								+ "\t" + (best.forward? best.hits.getFirst().targetStart + "\t" + best.hits.getLast().targetEnd : best.hits.getLast().targetStart + "\t" + best.hits.getFirst().targetEnd )
								+ "\t" + anz + "\t" + best.getInfo(gff) + "\t" + best.backup + "\t" + best.cut
								+ ( seq != null ?  
										"\t" + getGapCost(seq.length(), currentInfo.exonID.length) /*-(seq.length()*gapExtension+gapOpening+parts.length*INTRON_GAIN_LOSS)*/
										+ "\t" + ((int)-psa.getCost()) + "\t" + getScore(seq, seq)
										+ "\t" + (pos/(double)s1.length()) + "\t" + (id/(double)s1.length()) + "\t" + maxGap + "\t" + seq.length() + "\t" + pred.length()
										: "" )
								+ (acceptorSites == null ? "" : ("\t" + (best.A==0? "NA": decFormat.format(best.a/(double)best.A))) )
								+ (donorSites == null ? "" : ("\t" + (best.D==0? "NA": decFormat.format(best.d/(double)best.D))) )
								+ (donorSites == null ? "" : ("\t" + (best.I==0? "NA": decFormat.format(best.i/(double)best.I))) )
								+ (coverage == null ? "" : ("\t" + decFormat.format(best.Cov/(double)best.Len) + "\t" + best.minC))
								+ "\t" + best.similar(result.peek())
								+ "\n"
						);
						
						gff.newLine();
						
						gff.append(gffHelp);
						
						
						gff.flush();
						//protocol.flush();
					}
				}
				if( counts == 0 ) {
					protocol.append( "\n" );
					//protocol.flush();
				}
			}
		}
		
		private void filter( HashMap<Integer,ArrayList<Hit>> hash, int start, int end ) {
			if( start == 0 && end == -1 ) {
				return;
			}
			//filter out
			Hit h;
			if( verbose ) protocol.append("\n");
			for( int i = 0; i < currentInfo.exonID.length; i++ ) {
				ArrayList<Hit> current = hash.get(currentInfo.exonID[i]);
				if( verbose ) protocol.append(i + "\t" + currentInfo.exonID[i] + "\t#=" + (current==null?0:current.size()) + "\n");
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
				
				if( verbose ) protocol.append(i + "\t" + currentInfo.exonID[i] + "\t#=" + (current==null?0:current.size()) + "\n\n");
			}
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
				for( int i = currentInfo.exonID.length-1; i >= 0; i-- ) {
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
			for( int i = currentInfo.exonID.length-1; i >= 0; i-- ) {
				l = lines.get(currentInfo.exonID[i]);
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
						splice[i] = new int[anz][currentInfo.exonID.length-i][];
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
						int m=Math.min(i+MAX_GAP,currentInfo.exonID.length), startIdx=j+1, max = 0;
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
					if( i == currentInfo.exonID.length-1 ) {
						return o.lastAddScore;
					} else {
						return getGapCost(o.queryLength-o.queryEnd + revCumLength[i], currentInfo.exonID.length-1-i);
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
			int max = getCost( lines.get(currentInfo.exonID[i]).get(j), gap, false, i);
			
			ArrayList<Hit> ref = lines.get(currentInfo.exonID[k]);
			if( ref == null ) {
				return max;
			}
			Hit upstreamHit = lines.get(currentInfo.exonID[i]).get(j);
			
			int introns= Math.max(1,k-i);
			if( splice != null && splice[i][j][k-i] == null ) {
				//protocol.appendln(o.targetID + "\t" + i + "\t" + j + "\t" + (k-i) );
				try {
					splice[i][j][k-i] = new int[ref.size()];
					Arrays.fill(splice[i][j][k-i], NOT_COMPUTED );
				} catch( OutOfMemoryError e ) {
					protocol.append(currentInfo.exonID[i] + "\t" + length[i] + "\n" );
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
				if( d < introns*maxAllowedIntron ) {
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
			Hit o = lines.get(currentInfo.exonID[i]).get(j);
			if( used != null ) {
				if( used[i] == null ) {
					used[i] = new double[lines.get(currentInfo.exonID[i]).size()];
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
			int min=Math.min(i+MAX_GAP,currentInfo.exonID.length), startIdx=j+1;
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
					if( d < introns*maxAllowedIntron ) {
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
		 * @param lines the search {@link Hit}s
		 */
		private void detailedAnalyse( String chromosome, boolean forward, HashMap<Integer,ArrayList<Hit>> lines ) throws CloneNotSupportedException, IllegalArgumentException, WrongAlphabetException, IOException {
			if( verbose ) {
				protocol.append( "all hits " + chromosome + " " + (forward?"+":"-") +"\n");
				show(lines);
			}
			//TASK 1: reduce blast hit set	
			int max = forwardDP(forward, lines, false);
			HashMap<Integer, ArrayList<Hit>> filtered = reduce(chromosome, forward, lines, avoidPrematureStop, false, max*hitThreshold);
			
			if( verbose ) {
				protocol.append( "filtered hits " + chromosome + " " + (forward?"+":"-")+"\n" );
				show(filtered);
			}
			
			//XXX split region?
			ArrayList<Hit> current, all = new ArrayList<Hit>();
			for( int j=0; j<currentInfo.exonID.length; j++ ) {
				current = filtered.get(currentInfo.exonID[j]);
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
			for( int h = 0; h <currentInfo.exonID.length; h++ ) {
				ArrayList<Hit> list = lines.get(currentInfo.exonID[h]);
				protocol.append(h + "\t" + currentInfo.exonID[h]+"\t#=" + (list==null?0:list.size()) + "\n");
				for( int k = 0; list!= null && k < list.size(); k++ ) {
					protocol.append( list.get(k) +"\n" );
				}/**/
			}
			protocol.append("\n");
		}
		
		private HashMap<Integer, ArrayList<Hit>> reduce( String chromosome, boolean forward, HashMap<Integer,ArrayList<Hit>> lines, boolean avoidStop, boolean gap, double thresh ) throws IOException, WrongAlphabetException {
			used = new double[currentInfo.exonID.length][];
			//mark unused blast hits
			ArrayList<Hit> list;
			for( int i = 0; i < sums.length; i++ ) {
				list = lines.get(currentInfo.exonID[i]);
				for( int j = 0; sums[i] != null && j < sums[i].length; j++ ) {
					int c = getCost(list.get(j), gap, true, i);
					//System.out.println(i+"\t"+j + "\t" + (sums[i][j]+c - thresh));
					backTracking( forward, lines, i, j, sums[i][j]+c - thresh, null, null, gap );
				}
			}			
			
			if( splice != null ) {
				int[] anz = new int[currentInfo.exonID.length];
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
						list = lines.get(currentInfo.exonID[i]);
						if( list != null && list.size()>0 ) {
							changed = true;
						}
					}
					anz[i] = a;
				}
				if( verbose ) protocol.append("#new hits: " + Arrays.toString(anz) + "\n");
				
				if( changed ) {
					int[][][][] splice2 = new int[currentInfo.exonID.length][][][]; 
					int[] su, st, e, qS, qE, qL;
					for( int i = 0; i < used.length; i++ ) {
						if( anz[i] > 0 ) {
							splice2[i] = new int[anz[i]][currentInfo.exonID.length-i][];
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
					list = lines.get(currentInfo.exonID[i]);
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
						filtered.put( currentInfo.exonID[i], add );
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
			for( int i = 0; i < currentInfo.exonID.length; i++ ) {
				current = filtered.get(currentInfo.exonID[i]);
				oldSize[i] = current == null ? 0 : current.size();
			}
			
			String chr = seqs.get(chromosome);
			if( acceptorSites!= null && donorSites != null ) {
				for( int h = 0; h <currentInfo.exonID.length; h++ ) {
					ArrayList<Hit> list = filtered.get(currentInfo.exonID[h]);
					for( int k = 0; list!= null && k < list.size(); k++ ) {
						Hit hit = list.get(k);
						hit.prepareSpliceCandidates(chr,align, maxAllowedIntron);
					}
				}
			}
			
			//TASK 2: add missing CDS currentInfo.exonID
			int first = -1, last=currentInfo.exonID.length;
			int end, start, old = -1;
			for( int j=0; j<currentInfo.exonID.length; j++ ) {
				current = filtered.get(currentInfo.exonID[j]);
				if( current != null && current.size() > 0 ) {
					if( first < 0 ) {
						first = j;
					}
					last=j;
					if( old >= 0 && old+1 != j ) {//missing part
						ArrayList<Hit> previous = filtered.get(currentInfo.exonID[old]);
						int offset = current.size(), d;
						UnionFind uf = new UnionFind(offset+previous.size());
						for( int a = 0; a < current.size(); a++ ) {
							Hit x = current.get(a);
							for( int b = 0; b < previous.size(); b++ ) {
								Hit y = previous.get(b);
								if( forward ) {
									d = x.targetStart - y.targetEnd;
									if( d >= 0 && d < (j-old)*maxAllowedIntron ) {
										uf.union(a, offset+b);
									}
								} else {
									d = y.targetStart - x.targetEnd;
									if( d >= 0 && d < (j-old)*maxAllowedIntron ) {
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
								if( verbose ) protocol.append("INTERNAL\t" + (old+1) + ".." + j + "\t" + currentInfo.exonID[old+1] + ".." + currentInfo.exonID[j] + "\n");
								align(chromosome, forward, end, start, old+1, j, filtered, "internal", ae && de, ae && de );
							}
						}
					}
					old = j;
				}
			}
			
			//add upstream CDS part candidates
			if( first > 0 ) {
				if( verbose ) protocol.append("FIRST\t" + first + "\t" + currentInfo.exonID[first] + "\n");
				extend(chromosome, forward, filtered, first, true, "upstream");
			}
			//add downstream CDS part candidates
			if( last < currentInfo.exonID.length-1 ) {
				if( verbose ) protocol.append("LAST\t" + last + "\t" + currentInfo.exonID[last] + "\n");
				extend(chromosome, forward, filtered, last, false, "downstream");
			}
			
			//check
			boolean cut = false;
			for( int i = 0; i < currentInfo.exonID.length; i++ ) {
				current = filtered.get(currentInfo.exonID[i]);
				if( current != null && current.size()-oldSize[i] > 20000 ) {
					while( current.size() > oldSize[i] ) {
						current.remove( current.size()-1 );
					}
					cut=true;
				}
			}
			
			//TASK 3: find new solution
			//prepare splicing for DP
			for( int h = 0; h <currentInfo.exonID.length; h++ ) {
				ArrayList<Hit> list = filtered.get(currentInfo.exonID[h]);
				String r = geneName + (transcriptInfo==null? "":("_" + currentInfo.exonID[h]));
				String ref = cds.get(r);
				for( int k = 0; list!= null && k < list.size(); k++ ) {
					Hit hit = list.get(k);
					if( !hit.splice ) {
						hit.prepareSpliceCandidates(chr,align,maxAllowedIntron);
					}
					if( ref != null && ((h==0 && ref.charAt(0)=='M') || (h+1==currentInfo.exonID.length && ref.charAt(ref.length()-1)=='*')) ) {
						hit.setBorderConstraints(h==0 && ref.charAt(0)=='M', h+1==currentInfo.exonID.length && ref.charAt(ref.length()-1)=='*',align);//Laura Kelly: partial
					}
				}
			}
			
			//sort hits			
			sort(filtered, forward);
					
			//DP with splicing
			splice = new int[currentInfo.exonID.length][][][];
			int bestValue = forwardDP(forward, filtered, true);
			filtered = reduce(chromosome, forward, filtered, avoidPrematureStop, true, bestValue);
			if( verbose ) {			
				for( int h = 0; h <currentInfo.exonID.length; h++ ) {
					ArrayList<Hit> list = filtered.get(currentInfo.exonID[h]);
					protocol.append( currentInfo.exonID[h] + "\t" + Arrays.toString(sums[h]) + "\n" );
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
		
		private void extend( String chromosome, boolean forward, HashMap<Integer,ArrayList<Hit>> lines, int index, boolean upstream, String info ) throws IllegalArgumentException, WrongAlphabetException {
			int next, dir, end;
			if( upstream ) {
				dir = -1;
				end = -1;
			} else {
				dir=1;
				end = currentInfo.exonID.length;
			}
			int d = maxAllowedIntron/10, f=0;
			int startIndex, a = (upstream?0:d), stop;
			ArrayList<Hit> current;
			int[][] array = null;
			do {
				current = lines.get(currentInfo.exonID[index]);
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
				//protocol.append( "added " + currentInfo.exonID[i] + "\t" + lines.get(currentInfo.exonID[i]).size() );
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
			String r = geneName + (transcriptInfo==null? "":("_" + currentInfo.exonID[idx]));
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
				
				ArrayList<Hit> cand = lines.get(currentInfo.exonID[idx]);
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
							if( idx+1==currentInfo.exonID.length && a+1 < split.length ) {
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
						lines.put( currentInfo.exonID[idx], cand );
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
					String p=cds.get(geneName + (transcriptInfo==null?"":("_" + currentInfo.exonID[idx])) );
					if( p != null ) {
						m += p;
					}
					if( idx+1 < revParts[partEnd] ) {
						m += currentInfo.splitAA[idx];
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
						
			//"intron-loss" //TODO 1.5: pc>=1 & (de=F || ae=F)
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
								//TODO 1.5: (de && ae)=F || intron
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
													double v = align1.getBestValue(align2, l1/3, l2/3, aa);
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
				ArrayList<Hit> current = hits.get(currentInfo.exonID[i]);
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
			int score, minSplitReads, a, A, d, D, i, I;
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
				boolean[] match = new boolean[currentInfo.exonID.length];
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
				String ref = getRefProtein(transcriptName);
				if( ref.charAt(0)=='M' ) {
					c = res.hits.get(0);
					c.extend( c.part == currentInfo.exonID[0], false ); //START
				}
				if( ref.charAt(ref.length()-1)=='*' ) {
					c = res.hits.get(res.hits.size()-1);
					c.extend( false, c.part == currentInfo.exonID[currentInfo.exonID.length-1]); //STOP
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
				
				String res = (hits.getFirst().part == currentInfo.exonID[0]) + "\t" + protein.charAt(0) + "\t" 
							+ (hits.getLast().part == currentInfo.exonID[currentInfo.exonID.length-1]) + "\t" + protein.charAt(protein.length()-1) + "\t" 
							+ anz + "\t" + intronGain +"\t" + intronLoss;
		
				return res;
			}
			
			public void writeSummary( String geneName, String transcriptName, int i ) throws Exception {
				Hit first = hits.get(0);
				gff.append( first.targetID + "\tGeMoMa\t"+tag+"\t" );				
				if( forward ) {
					gff.append( first.targetStart + "\t" + hits.getLast().targetEnd );
				} else {
					gff.append( hits.getLast().targetStart + "\t" + first.targetEnd );
				}
				gff.append( "\t.\t" + (forward?"+":"-") + "\t.\tID=" + prefix+transcriptName + "_R" + i + ";ref-gene=" + geneName + ";aa="+getNumberOfAA() );
			}
			
			public int writeGFF( String transcriptName, int pred, StringBuffer sb ) throws Exception {
				a=d=A=D=i=I=0;
				minSplitReads = Integer.MAX_VALUE;
				
				int start=-1, end = -1, phase = -1, parts = 0;
				String id = null;
				boolean first = true, ae = false, de = false;
				int[][] donSites = null;
				int[][] cov = null;
				
				Len = Cov = Sum = 0;
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
							int idx = Arrays.binarySearch(cov, new int[]{start}, IntArrayComparator.comparator[2] );
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
								} else {
									//p<inter[0] || p>inter[1]
									//because of while above: p<inter[0]
									min = 0;
									
									p = Math.min(inter[0],end+1);
								}
								//System.out.println(p + "\t" + Arrays.toString(inter) + "\t" + covered + "\t" + min);
							}
							
							Cov +=covered;
						}
						Len +=l;
						if( covered == 0 ) {
							min=0;
						}
						minC = Math.min(minC, min);
						
						sb.append( id + "\tGeMoMa\tCDS\t" + start + "\t" + end + "\t.\t" + (forward?"+":"-") + "\t" +phase+ "\t"
								//+ "ID=" +pref+"_cds"+parts +";" //not necessary "If a feature does not have any subparts, then it does not formally need an ID" (http://gmod.org/wiki/GFF3)
								+ "Parent=" + pref 
								+ (acceptorSites==null || first?"":(";ae="+ae))
								+ (donorSites==null?"":(";de="+de))
								+ (cov == null?"":(";pc=" + decFormat.format(covered/(double)l)+";minCov=" + min))
								+ "\n"
						);
						
						D++;
						d += de ? 1 : 0;
						combineIntronStat(first, ae, donSites, forward, start, end, last);
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
					int idx = Arrays.binarySearch(cov, new int[]{start}, IntArrayComparator.comparator[2] );
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
					
					Cov +=covered;
				}
				Len +=l;
				if( covered == 0 ) {
					min=0;
				}
				minC = Math.min(minC, min);
				
				sb.append( id + "\tGeMoMa\tCDS\t" + start + "\t" + end + "\t.\t" + (forward?"+":"-") + "\t" +phase+ "\t"
						//+"ID=" +pref+"_cds"+parts+ ";"
						+ "Parent=" + pref
						+ (acceptorSites==null || first?"":(";ae="+ae))
						+ (cov == null?"":(";pc=" + decFormat.format(covered/(double)l)+";minCov=" + min))
						+ "\n"
				);
				combineIntronStat(first, ae, donSites, forward, start, end, last);
				parts++;
								
				return parts;
			}
			
			void combineIntronStat( boolean first, boolean ae, int[][] donSites, boolean forward, int start, int end, int last ) {
				if( !first ) {
					A++;
					a += ae ? 1 : 0;
					I++;
					if( donSites != null ) {
							int v = forward ? start : (end+1);
							//TODO faster
							int idx = Arrays.binarySearch( donSites[0], last );
							if( idx > 0 ) {
								while( idx>0 &&  donSites[0][idx-1] == last ) {
									idx--;
								}
							}
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
					} else {
						minSplitReads = 0;
					}
				}
			}
			
		}
	}
	
	static class MyAlignment extends Alignment { //TODO static?
		
		private boolean destroyed=false;
		private IntList first = new IntList(), second = new IntList();
		private TranscriptPredictor instance;
		
		public MyAlignment(AffineCosts costs,TranscriptPredictor instance) {
			super(costs);
			this.instance = instance;
		}
		
		private void check() {
			if( destroyed ) throw new RuntimeException("A MyAlignment instance ist destroyed probably due to a time-out.");
		}

		public double getBestValue( MyAlignment a2, int end1, int end2, char aa ) throws WrongAlphabetException {
			check();
			if( this.type == AlignmentType.GLOBAL && a2.type == AlignmentType.GLOBAL //global alignments
					&& this.startS1 == 0 && a2.startS2 == 0 //both from the beginning 
			) {
				first.clear();
				for( int i = 1; i < this.l1; i++ ) {
					if( this.d[0][i][end1] < this.d[1][i][end1] && this.d[0][i][end1] < this.d[2][i][end1] ) {
						first.add(i);
					}
				}
				
				second.clear();
				for( int i = 1; i < this.l1; i++ ) {
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
						int l = (this.l1-pos2)-pos1;//XXX ?
						if( l >= 0 ) {
							double cost = l> 0 ? this.aCosts.getGapCostsFor(l) : 0;//approx;
							current=this.d[0][pos1][end1] + a2.d[0][pos2][end2] + cost;
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
			check();
			if( this.type == type 
				&& this.startS1 == startS1 && this.l1 == endS1-startS1
				&& this.startS2 == startS2 && this.l2 == endS2-startS2
			) {
					int i = startS1;
					while( i < endS1 && this.s1.discreteVal(i)==s1.discreteVal(i) ) {
						i++;						
					}
					if( i == endS1 ) {
						i = startS2;
						while( i < endS2 && this.s2.discreteVal(i)==s2.discreteVal(i) ) {
							i++;						
						}
						if( i == endS2 ) {
							//same as last computed alignment
							return false;
						}
					}
			}
			instance.numberOfPairwiseAlignments++;
			
			/*
			Sequence oldS1=this.s1;
			Sequence oldS2=this.s2;
			int oldStartS1=this.startS1;
			int oldStartS2=this.startS2;
			int oldL1=this.l1;
			int oldL2=this.l2;
			try {*/
				return super.computeAlignment(type, s1, startS1, endS1, s2, startS2, endS2);
			/*} catch( Exception e ) {
				System.out.println("old:");
				System.out.println(oldS1);
				System.out.println(oldS2);
				System.out.println(oldStartS1);
				System.out.println(oldStartS2);
				System.out.println(oldL1);
				System.out.println(oldL2);
				System.out.println();
				System.out.println("new:");
				System.out.println(s1);
				System.out.println(s2);
				System.out.println(startS1);
				System.out.println(startS2);
				System.out.println(l1);
				System.out.println(l2);
				System.out.println();
				
				printMatrix(s1,s2);
				
				throw e;
			}*/
		}
		
		void clear() {
			s1=s2=null;
			startS1=startS2=l1=l2=-1;
			destroyed=true;
		}
	}
	
	public static enum Score {
		/**
		 * Trust the score that is given from the search algorithm.
		 */
		Trust,
		/**
		 * Trust the alignment, but re-compute the score based on given gap opening, gap extension and substitution costs.
		 */
		ReScore,
		/**
		 * Use the sequences and compute a new global alignment based on given gap opening, gap extension and substitution costs.
		 */
		ReAlign;
	}

	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getShortName(), 
					new FileParameter( "search results", "The search results, e.g., from tblastn or mmseqs", "tabular", true, new FileExistsValidator() ),
					new FileParameter( "target genome", "The target genome file (FASTA), i.e., the target sequences in the blast run. Should be in IUPAC code", "fasta,fas,fa,fna,fasta.gz,fas.gz,fa.gz,fna.gz", true, new FileExistsValidator(), true ),
					new FileParameter( "cds parts", "The query CDS parts file (protein FASTA), i.e., the CDS parts that have been searched in the target genome using for instance BLAST or mmseqs", "fasta,fa,fas,fna", true, new FileExistsValidator(), true ),
					new FileParameter( "assignment", "The assignment file, which combines CDS parts to proteins", "tabular", false, new FileExistsValidator() ),
					
					new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
							new FileParameter( "introns", "Introns (GFF), which might be obtained from RNA-seq", "gff,gff3", true, new FileExistsValidator(), true )
						), "introns", "", 0 ) ),
					new SimpleParameter( DataType.INT, "reads", "if introns are given by a GFF, only use those which have at least this number of supporting split reads", true, new NumberValidator<Integer>(1, Integer.MAX_VALUE), 1 ),
					new SimpleParameter( DataType.BOOLEAN, "splice", "if no intron is given by RNA-seq, compute candidate splice sites or not", true, true ),
					
					new ParameterSetContainer( "coverage", "", new ExpandableParameterSet( new SimpleParameterSet(	
						new SelectionParameter( DataType.PARAMETERSET, 
								new String[]{"UNSTRANDED", "STRANDED"},
								new Object[]{
									//unstranded coverage
									new SimpleParameterSet(
											new FileParameter( "coverage_unstranded", "The coverage file contains the unstranded coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true, new FileExistsValidator(), true )
									),
									//stranded coverage
									new SimpleParameterSet(
											new FileParameter( "coverage_forward", "The coverage file contains the forward coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true, new FileExistsValidator(), true ),
											new FileParameter( "coverage_reverse", "The coverage file contains the reverse coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true, new FileExistsValidator(), true )
									)
								},  "coverage", "experimental coverage (RNA-seq)", true
						)
					), "coverage", "", 0 ) ),
					
					new FileParameter( "genetic code", "optional user-specified genetic code", "tabular", false, new FileExistsValidator() ),
					new FileParameter( "substitution matrix", "optional user-specified substitution matrix", "tabular", false, new FileExistsValidator() ),
					
					new SimpleParameter( DataType.INT, "gap opening", "The gap opening cost in the alignment", true, 11 ),
					new SimpleParameter( DataType.INT, "gap extension", "The gap extension cost in the alignment", true, 1 ),
					new SimpleParameter( DataType.INT, "maximum intron length", "The maximum length of an intron", true, 15000 ),
					new SimpleParameter( DataType.BOOLEAN, "static intron length", "A flag which allows to switch between "
							+ "static intron length, which can be specified by the user and is identical for all genes, and "
							+ "dynamic intron length, which is based on the gene-specific maximum intron length in the reference organism plus the user given maximum intron length"
							, true, true ),
					new SimpleParameter( DataType.INT, "intron-loss-gain-penalty", "The penalty used for intron loss and gain", true, 25 ),
			
					new SimpleParameter( DataType.DOUBLE, "e-value", "The e-value for filtering blast results", true, 1E2 ),
					new SimpleParameter( DataType.DOUBLE, "contig threshold", "The threshold for evaluating contigs", true, new NumberValidator<Double>(0d, 1d), 0.4 ),
					new SimpleParameter( DataType.DOUBLE, "region threshold", "The threshold for evaluating regions", true, new NumberValidator<Double>(0d, 1d), 0.9 ),
					new SimpleParameter( DataType.DOUBLE, "hit threshold", "The threshold for adding additional hits", true, new NumberValidator<Double>(0d, 1d), 0.9 ),
					
					new SimpleParameter( DataType.INT, "predictions", "The (maximal) number of predictions per transcript", true, 10 ), 
					new FileParameter( "selected", "The path to list file, which allows to make only a predictions for the contained transcript ids. The first column should contain transcript IDs as given in the annotation. Remaining columns can be used to determine a target region that should be overlapped by the prediction, if columns 2 to 5 contain chromosome, strand, start and end of region", "tabular,txt", maxSize>-1, new FileExistsValidator() ), 
					new SimpleParameter( DataType.BOOLEAN, "avoid stop", "A flag which allows to avoid (additional) pre-mature stop codons in a transcript", true, true ),
					new SimpleParameter( DataType.BOOLEAN, "approx", "whether an approximation is used to compute the score for intron gain", true, true ),
					new SimpleParameter( DataType.BOOLEAN, "protein alignment", "whether a protein alignment between the prediction and the reference transcript should be computed. If so two additional attributes (iAA, pAA) will be added to predictions in the gff output. These might be used in GAF. However, since some transcripts are very long this can increase the needed runtime and memory (RAM).", true, true ),
			
					new SimpleParameter( DataType.STRING, "prefix", "A prefix to be used for naming the predictions", true, "" ),
					new SimpleParameter( DataType.STRING, "tag", "A user-specified tag for transcript predictions in the third column of the returned gff. It might be beneficial to set this to a specific value for some genome browsers.", true, TAG ),
					
					new SimpleParameter( DataType.BOOLEAN, "verbose", "A flag which allows to output a wealth of additional information per transcript", true, false ),
					new SimpleParameter( DataType.LONG, "timeout", "The (maximal) number of seconds to be used for the predictions of one transcript, if exceeded GeMoMa does not output a prediction for this transcript.", true, new NumberValidator<Long>((long) 0, maxTimeOut), timeOut ),
					
					new SimpleParameter( DataType.BOOLEAN, "sort", "A flag which allows to sort the search results", true, false ),
					new EnumParameter( Score.class, "A flag which allows to do nothing, re-score or re-align the search results", true, "Trust" )
			);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	public String getToolName() {
		return "GeneModelMapper";
	}
		
	public String getShortName() {
		return "GeMoMa";
	}

	public String getDescription() {
		return "builds gene models from search results";
	}

	public String getHelpText() {
		return 
			"This tool is the main part of, a homology-based gene prediction tool. GeMoMa builds gene models from search results (e.g. tblastn or mmseqs).\n\n"
				//typical usage
				+ "As first step, you should run **Extractor** obtaining *cds parts* and *assignment*. Second, you should run a search algorithm, e.g. **tblastn** or **mmseqs**, with *cds parts* as query. Finally, these search results are then used in **GeMoMa**. Search results should be clustered according to the reference genes. The most easiest way is to sort the search results accoring to the first column. If the search results are not sorted by default (e.g. mmseqs), you should the parameter *sort*.\n"
				//protein
				+ "If you like to run GeMoMa ignoring intron position conservation, you should blast protein sequences and feed the results in *query cds parts* and leave *assignment* unselected.\n\n"
				//RNA-seq
				+ "If you like to run GeMoMa using RNA-seq evidence, you should map your RNA-seq reads to the genome and run **ERE** on the mapped reads. For several reasons, spurious introns can be extracted from RNA-seq data. Hence, we recommend to run **DenoiseIntrons** to remove such spurious introns. Finally, you can use the obtained *introns* (and *coverage*) in GeMoMa.\n\n"
				//multiple predictions
				+ "If you like to obtain multiple predictions per gene model of the reference organism, you should set *predictions* accordingly. In addition, we suggest to decrease the value of *contig threshold* allowing GeMoMa to evaluate more candidate contigs/chromosomes.\n\n"
				//runtime
				+ "If you change the values of *contig threshold*, *region threshold* and *hit threshold*, this will influence the predictions as well as the runtime of the algorithm. The lower the values are, the slower the algorithm is.\n\n"
				//filter
				+ "You can filter your predictions using **GAF**, which also allows for combining predictions from different reference organismns.\n\n"	
				//UTR + rename
				+ "Finally, you can predict UTRs and rename predictions using **AnnotationFinalizer**.\n\n"
				//pipeline
				+ "If you like to run the complete GeMoMa pipeline and not only specific module, you can run the multi-threaded module **GeMoMaPipeline**."
			+ MORE;
	}
	
	static final String DEF_RES = "predicted annotation";
	
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", DEF_RES)
		};
	}
	
	@Override
	public ToolResult[] getTestCases( String path ) {
		try {
			return new ToolResult[]{new ToolResult(FileManager.readFile(path+File.separator+"tests/gemoma/xml/gemoma-test.xml"))};
		} catch( Exception e ) {
			e.printStackTrace();
			return null;
		}
	}
}