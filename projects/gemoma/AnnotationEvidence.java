package projects.gemoma;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import de.jstacs.DataType;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
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
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.IntList;
import projects.gemoma.Extractor.Gene;
import projects.gemoma.GeMoMa.IntArrayComparator;

/**
 * This class computes the (RNAseq) evidence for a given set of annotations.
 * 
 * @author Jens Keilwagen
 */
public class AnnotationEvidence implements JstacsTool {
	
	private HashMap<String, int[][][]> donorSites;
	private static HashMap<String, int[][]>[] coverage;
	private HashMap<String, String> seqs;
	
	public ToolResult run( ParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {
		//sequence
		seqs = Tools.getFasta(parameters.getParameterForName("genome").getValue().toString(),20,' ');
		
		//annotation
		HashMap<String, HashMap<String,Gene>> annotation = Extractor.readGFF( parameters.getParameterForName("annotation").getValue().toString(), null, protocol);

		//introns
		ExpandableParameterSet eps = (ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(2)).getValue();
		ArrayList<String> fName = new ArrayList<String>();
		for( int i = 0; i < eps.getNumberOfParameters(); i++ ) {
			Parameter y = ((ParameterSet)eps.getParameterAt(i).getValue()).getParameterAt(0);
			if( y.isSet() ) {
				fName.add(y.getValue().toString());
			}
		}
		if( fName.size()>0 ) {
			HashMap<String, int[][][]>[] res = GeMoMa.readIntrons( (Integer) parameters.getParameterForName("reads").getValue(), protocol, false, seqs, fName.toArray(new String[fName.size()]) );
			donorSites = res[0];
		} else {
			donorSites = null;
		}

		//coverage
		coverage = GeMoMa.readCoverage( (ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(4)).getValue(), protocol );
		
		//compute
		String[] chr = seqs.keySet().toArray(new String[seqs.size()]);
		Arrays.sort(chr);
		
		File file = GeMoMa.createTempFile("AnnotationEvidence");
		BufferedWriter w = new BufferedWriter( new FileWriter(file) );
		for( String c: chr ) {
			HashMap<String,Gene> current = annotation.get(c);
			int[][][] sites = donorSites.get(c);
			if( current!= null && current.size() > 0 ) {
				Iterator<Gene> it = current.values().iterator();
				while( it.hasNext() ) {
					Gene g = it.next();
					g.sortExons();
		
					int[][] cov = (coverage != null && coverage[g.strand==1?0:1]!= null) ? coverage[g.strand==1?0:1].get(c) : null;
					
					Iterator<Entry<String,IntList>> cds = g.transcript.entrySet().iterator();
					while( cds.hasNext() ) {
						Entry<String,IntList> e = cds.next();
						IntList parts = e.getValue();
						w.append( g.id + "\t" + g.strand + "\t" + e.getKey() + "\t" + parts.length() + "\t" );
						
						double tie=0;
						int last=-10;
						int covered=0, l=0;
						int[][] donSites = sites[g.strand==1?0:1];
						for( int j = 0;j < parts.length(); j++ ) {
							int[] part = g.exon.get(parts.get(j));
							
							//tie
							int v = g.strand==1 ? part[1] : (part[2]+1);
							int idx = Arrays.binarySearch( donSites[0], last );
							if( idx >= 0 ) {
								while( idx < donSites[0].length && donSites[0][idx] == last && donSites[1][idx] != v ) {
									idx++;
								}
								if( idx < donSites[0].length && donSites[0][idx] == last && donSites[1][idx] == v ) {
									tie++;
								}
							}
							
							//cov
							if( cov != null ) {
								//TODO?
								int start = part[1];//t.targetStart;
								int end = part[2];//t.targetEnd;
								l += end-start+1;
								
								idx = Arrays.binarySearch(cov, new int[]{start}, IntArrayComparator.comparator[2] );
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
											break outerloop;
										}
									}
									if( inter[0]<= p && p <=inter[1] ) {
										int h = Math.min(inter[1],end)+1;
										int a=h-p;
										covered+=a;
										
										p=h;
									} else {//p<inter[0] && p<=inter[1]
										p = Math.min(inter[0],end+1);
									}
									//System.out.println(p + "\t" + Arrays.toString(inter) + "\t" + covered + "\t" + min);
								}
							}
							
							last = g.strand==1 ? part[2]+1 : part[1];
						}
						if( parts.length()==1 ) {
							w.append( "NA" );
						} else {
							w.append( GeMoMa.decFormat.format( tie/(parts.length()-1d)) );
						}
						w.append( "\t" );
						if( coverage ==null ) {
							w.append( "NA" );
						} else {
							w.append( GeMoMa.decFormat.format(covered/(double)l) );
						}
						w.newLine();
					}
				}
			}
		}
		w.close();
		
		TextResult t = new TextResult("evidence", "Result", new FileParameter.FileRepresentation(file.getAbsolutePath()), "tabular", getToolName(), null, true);
		
		return new ToolResult("", "", null, new ResultSet(t), parameters, getToolName(), new Date());
	}
	
	public ParameterSet getToolParameters() {
		try{
			return new SimpleParameterSet(
					new FileParameter( "annotation", "The genome annotaion file (GFF)", "gff", true ),
					new FileParameter( "genome", "The genome file (FASTA), i.e., the target sequences in the blast run. Should be in IUPAC code", "fasta", true ),
					new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
							new FileParameter( "introns file", "Introns (GFF), which might be obtained from RNAseq", "gff", false )
						), "introns", "", 1 ) ),
					new SimpleParameter( DataType.INT, "reads", "if introns are given by a GFF, only use those which have at least this number of supporting split reads", true, new NumberValidator<Integer>(1, Integer.MAX_VALUE), 1 ),

					new ParameterSetContainer( "coverage", "", new ExpandableParameterSet( new SimpleParameterSet(	
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
								},  "coverage file", "experimental coverage (RNAseq)", true
						)
					), "coverage", "", 1 ) )
			);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	public String getToolName() {
		return "Annotation evidence";
	}
	
	public String getToolVersion() {
		return "1.3.2";
	}
	
	public String getShortName() {
		return "Annotation evidence";
	}

	public String getDescription() {
		return "computes the evidence for annotated coding sequences";
	}

	public String getHelpText() {
		return 
			"**What it does**\n\nThis tool is the main part of GeMoMa, a homology-based gene prediction tool. AnnotationEvidence computes for each annotated coding sequence the transcription intron evidence (tie) and the transcript percentage coverage (tpc) given mapped RNAseq evidence. Please use ERE to preprocess the mapped reads."
			+ "**References**\n\nFor more information please visit http://www.jstacs.de/index.php/GeMoMa or contact jens.keilwagen@julius-kuehn.de.\n"
				+"If you use this tool, please cite\n\n*Using intron position conservation for homology-based gene prediction.*\n Keilwagen et al., NAR, 2016, http://nar.oxfordjournals.org/content/44/9/e89";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "tabular", "tie"),
		};
	}
}