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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;

import de.jstacs.DataType;
import de.jstacs.io.FileManager;
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
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import projects.gemoma.GeMoMa.IntArrayComparator;

/**
 * Tries to denoise the intron file by removing spurious introns.
 * 
 * @author Jens Keilwagen
 *
 */
public class DenoiseIntrons extends GeMoMaModule {

	@Override
	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getShortName(), 
					new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
							new FileParameter( "introns", "Introns (GFF), which might be obtained from RNA-seq", "gff", true )
						), "introns", "", 1 ) ),
					
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
								},  "coverage", "experimental coverage (RNA-seq)", true
						)
					), "coverage", "", 1 ) ),
					
					new SimpleParameter( DataType.INT, "maximum intron length", "The maximum length of an intron", true, 15000 ),
					new SimpleParameter( DataType.DOUBLE, "minimum expression", "The threshold for removing introns", true, new NumberValidator<Double>(0d, 1d), 0.01 ),
					new SimpleParameter( DataType.INT, "context", "The context upstream a donor and donwstream an acceptor site that is used to determine the expression of the region", true, new NumberValidator<Integer>(0, 100), 10 ) 
			);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		int context = (Integer) parameters.getParameterForName("context").getValue();
		double min = (Double) parameters.getParameterForName("minimum expression").getValue();
		int maxIntron = (Integer) parameters.getParameterForName("maximum intron length").getValue();
		
		//read coverage
		HashMap<String, int[][]>[] coverage = GeMoMa.readCoverage((ExpandableParameterSet) parameters.getParameterForName("coverage").getValue(), protocol, false, false);
		
		//old
		/*
		GeMoMa.readCoverage(initialCoverage[2], args[0], protocol, false);
		int anz =1,  u = 1;
		HashMap<String, int[][]>[] coverage;
		if( anz > 0 ) {
			coverage = new HashMap[2];
			coverage[0] = GeMoMa.combine(true, initialCoverage[0], initialCoverage[2] );
			if( anz == u ) {
				coverage[1] = coverage[0];
				//GeMoMa.check( protocol, coverage[0], "coverage" );
			} else {
				coverage[1] = GeMoMa.combine(false, initialCoverage[1], initialCoverage[2] );
				//GeMoMa.check( protocol, coverage[0], "forward coverage" );
				//GeMoMa.check( protocol, coverage[1], "reverse coverage" );
			}
		} else {
			coverage = null;
		}*/
		
		ExpandableParameterSet introns = (ExpandableParameterSet) parameters.getParameterForName("introns").getValue();
		ArrayList<String> fName = new ArrayList<String>();
		for( int i = 0; i < introns.getNumberOfParameters(); i++ ) {
			Parameter y = ((ParameterSet)introns.getParameterAt(i).getValue()).getParameterAt(0);
			if( y.isSet() ) {
				fName.add(y.getValue().toString());
			}
		}
		String combined;
		if( fName.size()==1 ) {
			combined = fName.get(0);
		} else {
			combined = File.createTempFile("combined-introns",".gff").getAbsolutePath();
			String[] in = fName.toArray(new String[0]);
			CombineIntronFiles.combine(protocol, combined, in);
		}
		
		protocol.append("\nstart denoising\n");		
		BufferedReader r = new BufferedReader( new FileReader(combined) );
		File denoise = Tools.createTempFile("denoise");
		BufferedWriter w = new BufferedWriter( new FileWriter(denoise) );
		String line;
		int del=0, all=0;
		int oldStart=-1, oldEnd=-1, oldIdx=-9999;
		String oldChr=null;
		int[][][] cov=null, single = new int[1][][], both = new int[2][][];
		double e1=-100, i1=-100, i2=-100, e2=-100;
		boolean first=true;
		while( (line=r.readLine()) != null ) {
			if( line.equalsIgnoreCase("##FASTA") ) break; //http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
			if( line.length() == 0 || line.startsWith("#") ) {
				w.append(line);
				w.newLine();
				continue; 
			} else if( first ) {
				w.append(INFO + getShortName() + " " + getToolVersion() + "; ");
				String info = JstacsTool.getSimpleParameterInfo(parameters);
				if( info != null ) {
					w.append("SIMPLE PARAMETERS: " + info );
				}
				w.newLine();
				first = false;
			}
		
			Intron i = new Intron(line);
			int idx = i.getIndex();
			boolean same= i.chr.equals(oldChr) && idx == oldIdx;
			if( !same ) {
				if( idx>= 0 ) {
					single[0] = coverage[idx].get(i.chr);
					cov = single;
				} else {
					both[0] = coverage[0].get(i.chr);
					both[1] = coverage[1].get(i.chr);
					cov = both;
				}
			}
			
			if( !(same && i.start == oldStart) ) {
				e1 = getAvg(i.start-context,i.start,cov);
				//i1 = getAvg(i.start,i.start+context,cov);
			}
			if( !(same && i.end == oldEnd) ) {
				//i2 = getAvg(i.end-context,i.end,cov);
				e2 = getAvg(i.end,i.end+context,cov);
			}
			
			//filter
			if( i.end-i.start<maxIntron && i.reads/Math.max(e1, e2) >= min ) {
				w.append(i.toString());
				w.newLine();
			} else {
				//System.out.println(i.chr + "\t" + i.start + "\t" + i.end + "\t" + (i.end-i.start) + "\t" + e1 + "\t" + i1 + "\t" + i2 + "\t" + e2  + "\t" + i.reads + "\t" + (i.reads/Math.max(e1, e2)) );
				del++;
			}
			all++;
			
			oldChr = i.chr;
			oldStart = i.start;
			oldEnd = i.end;
			oldIdx = idx;
		}
		r.close();
		w.close();
		protocol.append( "remove introns: " + del + "/" + all + " = " + (del/(double)all) + "\n");
		
		return new ToolResult("", "", null, new ResultSet(new TextResult("denoised introns", "Result", new FileParameter.FileRepresentation(denoise.getAbsolutePath()), "gff", getToolName(), null, true)), parameters, getToolName(), new Date());
	}
	
	static double getAvg( int s, int e, int[][][] cov ) {
		//System.out.println("interval: " + s + "\t" + e);
		double sum = 0;
		for( int x = 0; x < cov.length; x++ ) {
			double c = 0;
			int p = s;
			int idx = Arrays.binarySearch(cov[x], new int[]{p}, IntArrayComparator.comparator[2] );
			if( idx < 0 ) {
				idx = -(idx+1);
				idx = Math.max(0, idx-1);
			}
			int[] inter = cov[x][idx];
			
			int stop=0;
			outerloop: while( p < e && stop<20) {
				while( p > inter[1] ) {
					idx++;
					if( idx < cov[x].length ) {
						inter = cov[x][idx];
					} else {
						break outerloop;
					}
				}
				//System.out.println(p + "\t" + Arrays.toString(inter) + "\t" + c );
				if( inter[0]<= p && p <=inter[1] ) {
					int h = Math.min(inter[1],e)+1;
					int a=h-p;
					c+=inter[2] * a;
					p=h;
				} else {
					//p<inter[0] || p>inter[1]
					//because of while above: p<inter[0]
					p = Math.min(inter[0],e+1);
				}
				
				stop++;
			}
			sum += (c / (double) (e-s));
		}
		return sum / cov.length;
	}
	
	static class Intron{
		String[] split;
		String chr;
		char strand;
		int start, end, reads;
		
		public Intron( String line ) {
			split = line.split("\t");
			chr=split[0];
			start = Integer.parseInt(split[3]);
			end = Integer.parseInt(split[4]);
			reads = Integer.parseInt(split[5]);
			strand = split[6].charAt(0);
		}
		
		public int getIndex() {
			switch( strand ) {
				case '+': return 0;
				case '-': return 1;
				default: return -1;
			}
		}

		@Override
		public String toString() {
			String res = "";
			for( int i = 0; i < split.length; i++ ) {
				res += (i==0?"":"\t") + split[i];
			}
			return res;
		}
	}

	@Override
	public String getToolName() {
		return "Denoise";
	}

	@Override
	public String getShortName() {
		return getToolName();
	}

	@Override
	public String getDescription() {
		return "removes spurious introns";
	}

	@Override
	public String getHelpText() {
		return "This module allows to analyze introns extracted by ERE. Introns with a large intron size or a low relative expression are possibly artefacts and will be removed."
				+ "The result of this module can be used in the module GeMoMa and AnnotationEvidence.";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", "denoised introns")
		};
	}

	@Override
	public ToolResult[] getTestCases( String path ) {
		try {
			return new ToolResult[]{new ToolResult(FileManager.readFile(path+File.separator+"tests/gemoma/xml/denoise-test.xml"))};
		} catch( Exception e ) {
			e.printStackTrace();
			return null;
		}
	}
}