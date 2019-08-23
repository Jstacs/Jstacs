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
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.IntList;
import projects.gemoma.Extractor.Gene;
import projects.gemoma.GeMoMa.IntArrayComparator;
import projects.gemoma.Tools.Ambiguity;

/**
 * This class computes the (RNA-seq) evidence for a given set of annotations.
 * 
 * @author Jens Keilwagen
 */
public class AnnotationEvidence extends GeMoMaModule {
	
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {		
		GeMoMa.fill(protocol, false, -1, 
				parameters.getParameterForName("genome").getValue().toString(),
				null, 
				(Integer) parameters.getParameterForName("reads").getValue(), (ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(2)).getValue(), 
				(ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(4)).getValue()
		);
		
		HashMap<String,Character> code = Tools.getCode(Tools.getInputStream(parameters.getParameterForName("genetic code"), "projects/gemoma/test_data/genetic_code.txt" ));
		code.put("NNN", 'X');//add for splitting at NNN (in align method)
				
		//annotation
		HashMap<String, HashMap<String,Gene>> annotation = Extractor.read( parameters.getParameterForName("annotation").getValue().toString(), null, protocol);

		//compute
		String[] chr = GeMoMa.seqs.keySet().toArray(new String[GeMoMa.seqs.size()]);
		Arrays.sort(chr);
		
		File file = Tools.createTempFile("AnnotationEvidence");
		BufferedWriter w = new BufferedWriter( new FileWriter(file) );
		w.append( "#gene id\tchr\tstart\tend\tstrand\ttranscript id\t#exons\ttie\ttpc" );
		w.newLine();
		File aFile = Tools.createTempFile("AnnotationEvidence");
		BufferedWriter annot = new BufferedWriter( new FileWriter(aFile) );
		annot.append("##gff-version 3");
		annot.newLine();
		StringBuffer nuc = new StringBuffer();
		for( String c: chr ) {
			String seq = GeMoMa.seqs.get(c);
			HashMap<String,Gene> current = annotation.get(c);
			int[][][] sites = GeMoMa.donorSites.get(c);
			if( current!= null && current.size() > 0 ) {
				Gene[] array = current.values().toArray(new Gene[0]);
				Arrays.sort(array);
				for( int a = 0; a < array.length; a++ ) {
					Gene g = array[a];
					g.sortExons();
					g.precompute();
		
					int[][] cov = (GeMoMa.coverage != null && GeMoMa.coverage[g.strand==1?0:1]!= null) ? GeMoMa.coverage[g.strand==1?0:1].get(c) : null;
					
					Iterator<Entry<String,IntList>> cds = g.transcript.entrySet().iterator();
					
					annot.append( c + "\t" + g.evidence + "\tgene\t" + g.start + "\t" + g.end + "\t.\t" + (g.strand==1?"+":"-") + "\t.\tID=" + g.id );
					annot.newLine();
					
					while( cds.hasNext() ) {
						Entry<String,IntList> e = cds.next();
						IntList parts = e.getValue();
						if( parts.length()>0 ) {
							w.append( g.id + "\t" + c + "\t" + g.start + "\t" + g.end + "\t"+ g.strand + "\t" + e.getKey() + "\t" + parts.length() + "\t" );
							
							double tie=0;
							int last=-10;
							int covered=0, l=0, idx;
							int[][] donSites = sites==null? null : sites[g.strand==1?0:1];
							nuc.delete( 0, nuc.length() );
							for( int j = 0;j < parts.length(); j++ ) {
								int[] part = g.exon.get(parts.get(j));
								
								if( g.strand> 0 ) {
									nuc.append( seq.substring( part[1]-1, part[2] ) );
								} else {
									nuc.insert(0, seq.substring( part[1]-1, part[2] ) );
								}
								
								
								//tie
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
										}
									}
								}
								
								//cov
								int start = part[1];//t.targetStart;
								int end = part[2];//t.targetEnd;
								l += end-start+1;
								if( cov != null ) {
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
											int z=h-p;
											covered+=z;
											
											p=h;
										} else {//p<inter[0] && p<=inter[1]
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
							
							String tieString = ( parts.length() == 1 ) ? "NA" : GeMoMa.decFormat.format( tie/(parts.length()-1d));
							String tpcString = ( GeMoMa.coverage == null ) ? "NA" : GeMoMa.decFormat.format(covered/(double)l);
							
							w.append( tieString + "\t" + tpcString );
							w.newLine();
							
							//TODO annotation
							annot.append( c + "\t" + g.evidence + "+AnnotationEvidence\t"+"prediction"/*TODO*/+"\t" + g.start + "\t" + g.end + "\t.\t" + (g.strand==1?"+":"-") + "\t.\tID=" + e.getKey() + ";Parent="+g.id +";tie=" + tieString + ";tpc="+ tpcString + ";AA=" + (l/3) + ";start=" + aa.charAt(0) + ";stop=" + aa.charAt(aa.length()-1) );
							annot.newLine();
							for( int j = 0;j < parts.length(); j++ ) {
								int[] part = g.exon.get(parts.get(j));
								annot.append( c + "\t" + g.evidence + "\tCDS\t" + part[1] + "\t" + part[2] + "\t.\t" + (g.strand==1?"+":"-") + "\t.\tParent=" + e.getKey() );
								annot.newLine();
							}
						}
					}
				}
			}
		}
		w.close();
		annot.close();
		
		ArrayList<TextResult> res = new ArrayList<TextResult>();
		res.add( new TextResult(defResult, "Result", new FileParameter.FileRepresentation(file.getAbsolutePath()), "tabular", getToolName(), null, true) );
		if( ((Boolean) parameters.getParameterForName("annotation output").getValue()) ) {
			res.add( new TextResult("annotation with attributes", "Result", new FileParameter.FileRepresentation(aFile.getAbsolutePath()), "gff", getToolName(), null, true) );
		}
		return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
	}
	
	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getShortName(),
					new FileParameter( "annotation", "The genome annotation file (GFF)", "gff", true ),
					new FileParameter( "genome", "The genome file (FASTA), i.e., the target sequences in the blast run. Should be in IUPAC code", "fasta", true ),
					new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
							new FileParameter( "introns file", "Introns (GFF), which might be obtained from RNA-seq", "gff", false )
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
								},  "coverage file", "experimental coverage (RNA-seq)", true
						)
					), "coverage", "", 1 ) ),
					
					new SimpleParameter( DataType.BOOLEAN, "annotation output", "if the annotation should be returned with attributes tie, tpc, and AA", true, false ),
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
			"**What it does**\n\nThis tool computes for each annotated coding sequence the transcription intron evidence (tie) and the transcript percentage coverage (tpc) given mapped RNA-seq evidence. "
			+ "All predictions of the annotation are used. The predictions are not filtested for internal stop codons, missing start or stop codons, frame-shifts, ...\n\n"
			+ "In addition, it allows to add attributes to the annotation, e.g. tie, tpc, AA, start, stop, that can be used if the annotation should be used in *GAF*.\n\n"
			+ "Please use *ERE* to preprocess the mapped reads.\n\n"
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
	public ToolResult[] getTestCases() {
		// TODO missing test cases
		return null;
	}
}