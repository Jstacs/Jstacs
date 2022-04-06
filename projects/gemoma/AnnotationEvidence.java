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
import de.jstacs.utils.IntList;
import projects.gemoma.Extractor.Gene;
import projects.gemoma.Extractor.Part;
import projects.gemoma.Extractor.Transcript;
import projects.gemoma.Tools.Ambiguity;

/**
 * This class computes the (RNA-seq) evidence for a given set of annotations.
 * 
 * @author Jens Keilwagen
 */
public class AnnotationEvidence extends GeMoMaModule {
	
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp ) throws Exception {		
		if( GeMoMa.seqs == null  ) {
			GeMoMa.fill(protocol, false, -1, 
				parameters.getParameterForName("genome").getValue().toString(),
				null, 
				(Integer) parameters.getParameterForName("reads").getValue(), 
				(ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(3)).getValue(), 
				(ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(5)).getValue()
			);
		}
		String tag = parameters.getParameterForName("tag").getValue().toString();
		
		HashMap<String,Character> code = Tools.getCode(Tools.getInputStream(parameters.getParameterForName("genetic code"), "projects/gemoma/test_data/genetic_code.txt" ));
		code.put("NNN", 'X');//add for splitting at NNN (in align method)
				
		//annotation
		HashMap<String, HashMap<String,Gene>> annotation = Extractor.read( false, parameters.getParameterForName("annotation").getValue().toString(), null, protocol);

		//compute
		String[] chr = GeMoMa.seqs.keySet().toArray(new String[GeMoMa.seqs.size()]);
		Arrays.sort(chr);
		
		File file = Tools.createTempFile("AnnotationEvidence",temp);
		BufferedWriter w = new BufferedWriter( new FileWriter(file) );
		w.append( "#gene id\tchr\tstart\tend\tstrand\ttranscript id\t#exons\ttie\ttpc\tminCov\tavgCov\tminSplitReads\tnps" );
		w.newLine();
		File aFile = Tools.createTempFile("AnnotationEvidence",temp);
		BufferedWriter annot = new BufferedWriter( new FileWriter(aFile) );
		annot.append("##gff-version 3");
		annot.newLine();
		annot.append(GeMoMa.INFO + getShortName() + " " + getToolVersion() + "; ");
		String info = JstacsTool.getSimpleParameterInfo(parameters);
		if( info != null ) {
			annot.append("SIMPLE PARAMETERS: " + info );
		}
		annot.newLine();
		StringBuffer nuc = new StringBuffer();
		
		boolean coverage = GeMoMa.coverage!=null;
		boolean introns = GeMoMa.donorSites!=null;
		
		HashMap<String, String> attr = new HashMap<String,String>();
		StringBuffer old = new StringBuffer();
		int perfect=0;
		for( String c: chr ) {
			String seq = GeMoMa.seqs.get(c);
			HashMap<String,Gene> current = annotation.get(c);
			int[][][] sites = introns ? GeMoMa.donorSites.get(c) : null;
			if( current!= null && current.size() > 0 ) {
				Gene[] array = current.values().toArray(new Gene[0]);
				if( array.length==1 ) {
					array[0].precompute();
				} else {
					Arrays.sort(array);
				}
				for( int a = 0; a < array.length; a++ ) {
					Gene g = array[a];
					
					boolean hasTranscript = false;
					Iterator<Transcript> ti = g.transcript.values().iterator();
					while( ti.hasNext() ) {
						if( ti.next().b.length()>0 ) {
							hasTranscript=true;
							break;
						}
					}

					if( hasTranscript ) {
						g.sortExons();
						
						int[][] cov = (coverage && GeMoMa.coverage[g.strand==1?0:1]!= null) ? GeMoMa.coverage[g.strand==1?0:1].get(c) : null;
						
						Iterator<Entry<String,Transcript>> cds = g.transcript.entrySet().iterator();
						
						annot.append( c + "\t" + g.evidence + "\tgene\t" + g.start + "\t" + g.end + "\t.\t" + (g.strand==1?"+":"-") + "\t.\tID=" + g.id );
						annot.newLine();
						
						while( cds.hasNext() ) {
							Entry<String,Transcript> e = cds.next();
							Transcript t = e.getValue();
							IntList parts = t.b;
							if( parts.length()>0 ) {
								int tStart = g.strand>0 ? g.exon.get(parts.get(0))[1] : g.exon.get(parts.get(parts.length()-1))[1];
								int tEnd = g.strand>0 ? g.exon.get(parts.get(parts.length()-1))[2] : g.exon.get(parts.get(0))[2];
								if( t.add != null ) {
									for( int j = 0; j<t.add.size(); j++ ) {
										String[] s = t.add.get(j);
										tStart = Math.min(tStart, Integer.parseInt(s[3]) );
										tEnd = Math.max(tEnd, Integer.parseInt(s[4]));
									}
								}
								
								w.append( g.id + "\t" + c + "\t" + tStart + "\t" + tEnd + "\t"+ g.strand + "\t" + e.getKey() + "\t" + parts.length() + "\t" );
								
								double tie=0;
								int minSplitReads=Integer.MAX_VALUE;
								int last=-10;
								int covered=0, l=0, idx;
								int[][] donSites = sites==null? null : sites[g.strand==1?0:1];
								nuc.delete( 0, nuc.length() );
								double sum=0;
								int min=Integer.MAX_VALUE;
								for( int j = 0;j < parts.length(); j++ ) {
									int[] part = g.exon.get(parts.get(j));
									
									if( g.strand> 0 ) {
										nuc.append( seq.substring( part[1]-1, part[2] ) );
									} else {
										nuc.insert(0, seq.substring( part[1]-1, part[2] ) );
									}
									
									
									//tie
									if( j > 0 ) {
										int currentSplitReads=0;
										if( donSites != null ) {
											int v = g.strand==1 ? part[1] : (part[2]+1);
											
											idx = Arrays.binarySearch( donSites[0], last );
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
													tie++;
													currentSplitReads=donSites[2][idx];
												}
											}
										}
										minSplitReads=Math.min(minSplitReads,currentSplitReads);
									}
									
									//cov
									int start = part[1];//t.targetStart;
									int end = part[2];//t.targetEnd;
									l += end-start+1;
									if( cov != null ) {
										idx = Arrays.binarySearch(cov, new int[]{start}, GeMoMa.IntArrayComparator.comparator[2] );
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
													min=0;
													break outerloop;
												}
											}
											if( inter[0]<= p && p <=inter[1] ) {
												int h = Math.min(inter[1],end)+1;
												int z=h-p;
												covered+=z;
												
												sum+=inter[2] * z;
												min = Math.min(min, inter[2]);
												
												p=h;
											} else {//p<inter[0] && p<=inter[1]
												min=0;
												p = Math.min(inter[0],end+1);
											}
											//System.out.println(p + "\t" + Arrays.toString(inter) + "\t" + covered + "\t" + min);
										}
									}
									
									last = g.strand==1 ? part[2]+1 : part[1];
								}
								
								//check
								String cod;
								if( g.strand < 0 ) {
									cod = Tools.rc(nuc.toString());
								} else {
									cod = nuc.toString();
								}
								String aa = Tools.translate(0, cod, code, false, Ambiguity.AMBIGUOUS);
								
								int preMatureStops=0;
								for( int j=0; j < aa.length()-1; j++ ) {
									if( aa.charAt(j)=='*') preMatureStops++;
								}
								
								String tieString = ( parts.length() == 1 ) ? "NA" : GeMoMa.decFormat.format( tie/(parts.length()-1d));
								String tpcString = ( GeMoMa.coverage == null ) ? "NA" : GeMoMa.decFormat.format(covered/(double)l);
								String avgCovString = ( GeMoMa.coverage == null ) ? "NA" : GeMoMa.decFormat.format(sum/(double)l);
								String minCovString = ( GeMoMa.coverage == null ) ? "NA" : (""+min);
								String minSplitReadsString = ( parts.length() == 1 ) ? "NA" : (""+minSplitReads);
								
								if( (tie==parts.length()-1) && covered==l ) {
									perfect++;
								}
								
								w.append( tieString + "\t" + tpcString + "\t" + minCovString + "\t" + avgCovString + "\t" + minSplitReadsString + "\t" + preMatureStops );
								w.newLine();
								
								//annotation
								annot.append( c + "\t" + g.evidence + "+AnnotationEvidence\t"+tag+"\t" + tStart + "\t" + tEnd + "\t.\t" + (g.strand==1?"+":"-") + "\t.\t");
								//attributes
								String[] att = t.attributes==null ? EMPTY : t.attributes.split(";");
								for( int j = 0; j < att.length; j++ ) {
									att[j]=att[j].trim();
								}
								boolean gtf = att.length>=2 && (att[0].startsWith("gene_id") || att[1].startsWith("transcript_id"));
								if( gtf ) {
									//GTF -> GFF
									for( int j = 0; j < att.length; j++ ) {
										att[j]=att[j].replace(" \"", "=\"");
									}
									if( att[0].startsWith("gene_id") ) {
										att[0]=att[0].replace("gene_id", "Parent" );
										att[1]=att[1].replace("transcript_id", "ID" );
									} else {
										att[1]=att[1].replace("gene_id", "Parent" );
										att[0]=att[0].replace("transcript_id", "ID" );
									}
								}
								attr.clear();
								for( int j = 0; j < att.length; j++ ) {
									int pos = att[j].indexOf('=');
									attr.put(att[j].substring(0,pos), att[j].substring(pos+1));
								}
								
								old.delete(0, old.length());
								
								add(old,attr,"ID",e.getKey());
								add(old,attr,"Parent",g.id);
								add(old,attr,"aa",""+(l/3));
								if( introns ) {
									add(old,attr,"tie",tieString);
									add(old,attr,"minSplitReads",minSplitReadsString);
								}
								if( coverage ) {
									add(old,attr,"tpc",tpcString);
									add(old,attr,"minCov",""+min);
									add(old,attr,"avgCov",avgCovString);
								}
								add(old,attr,"nps",""+preMatureStops);
								add(old,attr,"start",""+aa.charAt(0));
								add(old,attr,"stop",""+aa.charAt(aa.length()-1));
								add(old,attr,"ce",""+parts.length());
								
								//write attributes in the same order as in the input file
								for( int j = 0; j < att.length; j++ ) {
									int pos = att[j].indexOf('=');
									String key = att[j].substring(0,pos);
									annot.append(key+"=" + attr.remove(key) + ";");
								}
								//add additional attributes
								String[] k = attr.keySet().toArray(EMPTY);
								Arrays.sort(k);
								for( String key: k ) {
									annot.append(key+"=" + attr.remove(key) + ";");
								}
								if( old.length()>0 ) annot.append("#old: " + old);
								
								/*
								annot.append( "ID=" + e.getKey() + ";Parent="+g.id + ";"
											+ "aa=" + (l/3) 
											+ (introns?";tie=" + tieString:"") 
											+ (coverage?(";tpc="+ tpcString + ";minCov=" + min + ";avgCov=" + avgCovString):"")
											+ ";nps=" + preMatureStops										
											+ ";start=" + aa.charAt(0) + ";stop=" + aa.charAt(aa.length()-1) );
								*/
								annot.newLine();
								if( t.add != null ) {
									for( int j = 0; j < t.add.size(); j++ ) {
										String[] split = t.add.get(j);
										for( int m = 0; m < split.length; m++ ) {
											annot.append( (m==0?"":"\t") + split[m] );
										}
										annot.newLine();
									}
								}
								for( int j = 0;j < parts.length(); j++ ) {
									int[] part = g.exon.get(parts.get(j));
									annot.append( c + "\t" + g.evidence + "\tCDS\t" + part[1] + "\t" + part[2] + "\t.\t" + (g.strand==1?"+":"-") + "\t"+(part[3]==Part.NO_PHASE?".":part[3])+"\tParent=" + e.getKey() );
									annot.newLine();
								}
							}
						}
					}
				}
			}
		}
		w.close();
		annot.close();
		
		protocol.append("number of detected transcripts with very good RNA-seq evidence (tpc==1 && (tie==1 || tie==NA)): " + perfect +"\n");
		
		ArrayList<TextResult> res = new ArrayList<TextResult>();
		res.add( new TextResult(defResult, "Result", new FileParameter.FileRepresentation(file.getAbsolutePath()), "tabular", getToolName(), null, true) );
		if( ((Boolean) parameters.getParameterForName("annotation output").getValue()) ) {
			res.add( new TextResult("annotation with attributes", "Result", new FileParameter.FileRepresentation(aFile.getAbsolutePath()), "gff", getToolName(), null, true) );
		}
		return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
	}
	
	private static final String[] EMPTY = new String[0];
	
	private static void add( StringBuffer old, HashMap<String,String> attr, String key, String newValue ) {
		String oldValue = attr.get(key);
		if( oldValue != null ) {
			if( oldValue.equals(newValue) ) return;
			old.append(key + "="  + oldValue+";");
		}
		attr.put(key, newValue);
	}
	
	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getShortName(),
					new FileParameter( "annotation", "The genome annotation file (GFF,GTF)", "gff,gff3,gtf,gff.gz,gff3.gz,gtf.gz", true, new FileExistsValidator(), true ),
					new SimpleParameter( DataType.STRING, "tag", "A user-specified tag for transcript predictions in the third column of the returned gff. It might be beneficial to set this to a specific value for some genome browsers.", true, GeMoMa.TAG ),
					new FileParameter( "genome", "The target genome file (FASTA). Should be in IUPAC code", "fasta,fas,fa,fna,fasta.gz,fas.gz,fa.gz,fna.gz", true, new FileExistsValidator(), true ),
					new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
							new FileParameter( "introns file", "Introns (GFF), which might be obtained from RNA-seq", "gff,gff3", false, new FileExistsValidator(), true )
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
											new FileParameter( "coverage_unstranded", "The coverage file contains the unstranded coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true, new FileExistsValidator() )
									),
									//stranded coverage
									new SimpleParameterSet(
											new FileParameter( "coverage_forward", "The coverage file contains the forward coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true, new FileExistsValidator() ),
											new FileParameter( "coverage_reverse", "The coverage file contains the reverse coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true, new FileExistsValidator() )
									)
								},  "coverage file", "experimental coverage (RNA-seq)", true
						)
					), "coverage", "", 1 ) ),
					
					new SimpleParameter( DataType.BOOLEAN, "annotation output", "if the annotation should be returned with attributes tie, tpc, and aa", true, true ),
					new FileParameter( "genetic code", "optional user-specified genetic code", "tabular", false )
			);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	public String getToolName() {
		return "Annotation evidence";
	}
	
	public String getShortName() {
		return "AnnotationEvidence";
	}

	public String getDescription() {
		return "computes the evidence for annotated coding sequences";
	}

	public String getHelpText() {
		return 
			"This tool adds attributes to the annotation, e.g., tie, tpc, aa, start, stop."
			+ " These attributes can be used, for instance, if the annotation is used in **GAF**."	
			+ " All predictions of the annotation are used. The predictions are not filtered for internal stop codons, missing start or stop codons, frame-shifts, ... ."
			+ " Please use **ERE** to preprocess the mapped reads."
			+ MORE;
	}

	private static final String defResult = "evidence";
	
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "tabular", defResult)
		};
	}

	@Override
	public ToolResult[] getTestCases( String path ) {
		// TODO missing test cases
		return null;
	}
}