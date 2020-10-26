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
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;

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
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;

/**
 * This class computes the (RNA-seq) evidence for a given set of annotations.
 * 
 * @author Jens Keilwagen
 */
public class TranscribedCluster extends GeMoMaModule {
	
	private static HashMap<String, int[][][]>[] spliceSites;
	private static HashMap<String, int[][]>[] coverage;
	private static int tc = 0;
	
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp ) throws Exception {
		//sequence
		HashMap<String, String> seqs = Tools.getFasta(parameters.getParameterForName("genome").getValue().toString(),20);
				
		//introns
		ExpandableParameterSet eps = (ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(1)).getValue();
		ArrayList<String> fName = new ArrayList<String>();
		for( int i = 0; i < eps.getNumberOfParameters(); i++ ) {
			Parameter y = ((ParameterSet)eps.getParameterAt(i).getValue()).getParameterAt(0);
			if( y.isSet() ) {
				fName.add(y.getValue().toString());
			}
		}
		if( fName.size()>0 ) {
			spliceSites = GeMoMa.readIntrons( (Integer) parameters.getParameterForName("reads").getValue(), protocol, false, seqs, null, fName.toArray(new String[fName.size()]) );
		} else {
			spliceSites = null;
		}

		//coverage
		coverage = GeMoMa.readCoverage( (ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(3)).getValue(), protocol, false, true );
	
		int minGap = (Integer) parameters.getParameterForName("minimal gap").getValue();

		//compute
		String[] chr = seqs.keySet().toArray(new String[seqs.size()]);
		Arrays.sort(chr);

		File file = Tools.createTempFile("TranscribedCluster",temp);
		BufferedWriter w = new BufferedWriter( new FileWriter(file) );
		w.append("##gff-version 3");
		w.newLine();
		int[][] cov;
		for( String c : chr ) {
			if( coverage[0] == coverage[1] ) {
				//unstranded
				cov = (coverage != null && coverage[0]!= null) ? coverage[0].get(c) : null;
				identify( protocol, c, -1, cov, minGap, w );
			} else {
				//stranded
				cov = (coverage != null && coverage[0]!= null) ? coverage[0].get(c) : null;
				identify( protocol, c, 0, cov, minGap, w );
				cov = (coverage != null && coverage[1]!= null) ? coverage[1].get(c) : null;
				identify( protocol, c, 1, cov, minGap, w );
			}
		}
		w.close();
		TextResult t = new TextResult(defResult, "Result", new FileParameter.FileRepresentation(file.getAbsolutePath()), "gff", getToolName(), null, true);
		
		return new ToolResult("", "", null, new ResultSet(t), parameters, getToolName(), new Date());
	}
	
	static char[] strand = {'?','+','-'};
	
	private static void identify( Protocol p, String chr, int str, int[][] cov, int minGap, BufferedWriter w ) throws IOException {
		//compute
		int start = -100, end = -100, max = -100, s = str==-1 ? -2 : str;
		int[][][] donor = null, acceptor = null;
		if( spliceSites!= null ) {
			donor = spliceSites[0].get(chr);
			acceptor = spliceSites[1].get(chr);
		}
		int f = 0, r = 0;
		if( cov != null ) {
			for( int c = 0; c < cov.length; c++ ) {
				if( end+1+minGap >= cov[c][0] ) {
					if( end < cov[c][1] ) {
						end = cov[c][1];
						
						//splicing
						if( str <= 0) {
							//forward
							int old = f;
							if( donor != null && donor[0] != null ) {
								while( f < donor[0][0].length && donor[0][0][f] <= end ) {
		p.append( chr + ": " + donor[0][0][f] + " -> " + donor[0][1][f] + "\n");
									if( end < donor[0][1][f] ) {
										end = donor[0][1][f];
									}
									f++;
								}
							}
							if( old < f & str == -1 ) {
								if( s==0 || s == -2 ) {
									s=0;
								} else {
									s=-1;
								}
							}
						}
						if( str != 0 ) {
							//reverse
							int old = r;
							if( acceptor != null && acceptor[0] != null ) {
								while( r < acceptor[1][1].length && acceptor[1][1][r] <= end ) {
		p.append( chr + ": " + acceptor[1][0][r] + " <- " + acceptor[1][1][r] + "\n" );
									if( end < acceptor[1][0][r] ) {
										end = acceptor[1][0][r];
									}
									r++;
								}
							}
							if( old < r & str == -1) {
								if( s==1 || s == -2 ) {
									s=1;
								} else {
									s=-1;
								}
							}
						}
					}
					if( max < cov[c][2] ) {
						max = cov[c][2];
					}
				} else {
					if( end > 0 ) {
						p.append("=> " + chr + "\t" + start + "\t" + end + "\t" + strand[Math.max(0, s+1)] + "\t" + max+ "\n");
						w.append(chr+"\tRNAseq\ttranscribed_cluster\t" + start + "\t"+ end+ "\t" + max + "\t" + strand[Math.max(0, s+1)] + "\t.\tID=tc"+tc++ );
						w.newLine();
					}
					start = cov[c][0];
					end = cov[c][1];
					max = cov[c][2];
					s=str==-1 ? -2 : str;
				}
			}
			p.append("=> " + chr + "\t" + start + "\t" + end + "\t" + strand[Math.max(0, s+1)] + "\t" + max+ "\n");
			w.append(chr+"\tRNAseq\ttranscribed_cluster\t" + start + "\t"+ end+ "\t" + max + "\t" + strand[Math.max(0, s+1)] + "\t.\tID=tc"+tc++ );
			w.newLine();
		}
	}
	
	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getShortName(),
					new FileParameter( "genome", "The genome file (FASTA), i.e., the target sequences in the blast run. Should be in IUPAC code", "fasta", true ),
					
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
					new SimpleParameter( DataType.INT, "minimal gap", "the minimal gap between two transcribed clusters, otherwise these will be merged", true, new NumberValidator<Integer>(0, Integer.MAX_VALUE), 50 )
			);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	public String getToolName() {
		return "Transcribed Cluster";
	}
	
	public String getShortName() {
		return "TranscribedCluster";
	}

	public String getDescription() {
		return "computes the transcribed clusters from RNA-seq evidence";
	}

	public String getHelpText() {
		return 
			"**What it does**\n\nThis tool computes ... .\n\n" //TODO
			+ MORE;
	}

	private static final String defResult = "transcribedCluster";
	
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", defResult),
		};
	}
	
	@Override
	public ToolResult[] getTestCases( String path ) {
		// TODO missing test cases
		return null;
	}
}