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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.IntList;


/**
 * This class allows for comparing GFFs.
 * It determines for each predicted gene model the best matching annotated gene model by evaluating the F1 measure.
 * 
 * @author Jens Keilwagen
 */
public class CompareTranscripts implements JstacsTool {
	
	@Override
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		String fName = (String) parameters.getParameterForName("assignment").getValue();
		HashMap<String, String[]> gene = Tools.getAlias(fName, 1, 0, 2);
		
		HashMap<String,Annotation> truth = readGFF( parameters.getParameterForName("annotation").getValue().toString(), false );
		HashMap<String,Annotation> prediction = readGFF( parameters.getParameterForName("prediction").getValue().toString(), true );//standard
		protocol.append( "truth: " + truth.size()  +"\n");
		protocol.append( "prediction: " + prediction.size() +"\n" );
		
		File f = GeMoMa.createTempFile("CompareTranscript");
		bestHit(gene, /*TODO*/"_R", null, truth, prediction, f,protocol);
		protocol.append("problem: " + Arrays.toString(problem));

		return new ToolResult("", "", null, new ResultSet(
				new TextResult("comparison", "Result", new FileParameter.FileRepresentation(f.getAbsolutePath()), "tabular", getToolName(), null, true)
			), parameters, getToolName(), new Date());
	}
	
	private static int problem[] = new int[2];
	
	private static String par = "Parent=";
	
	private static HashMap<String, String> info = new HashMap<String, String>();
		
	/**
	 * 
	 * @param input the input file name (GFF)
	 * @param type the type to be searched, e.g. &quot;CDS&quot;
	 * @param skip
	 * @param selected
	 * @param sep the separator for the info field in the GFF
	 * @return a {@link HashMap} containing all relevant {@link Annotation}s from the GFF
	 * 
	 * @throws Exception
	 */
	private static HashMap<String,Annotation> readGFF( String input, boolean information ) throws Exception {
		String line;
				
		HashMap<String,Annotation> annot = new HashMap<String,Annotation>();
		Annotation current;
		BufferedReader r = new BufferedReader( new FileReader(input) );
		while( (line=r.readLine()) != null ) {
			if( line.equalsIgnoreCase("##FASTA") ) break; //http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
			if( line.length() == 0 || line.startsWith("#") ) continue; 
			
			String[] split = line.split("\t");
			if( split[2].equalsIgnoreCase("CDS") ) {
				
				int idx = split[8].indexOf(par);
				if( idx >= 0 ) {
					idx+=par.length();
				} else {
					idx = split[8].indexOf("ID=");
					if( idx<0 ) {
						r.close();
						throw new IllegalArgumentException("corrupt line: " + line);
					}
					idx+=3;
				}
				int h = split[8].indexOf(';',idx);
				String[] t = split[8].substring(idx, h<0?split[8].length():h).toUpperCase().split(",");
				
				for( String transcriptID: t ) {
					current = annot.get(transcriptID);
					if( current == null ) {
						current = new Annotation( transcriptID, split[0], split[6].trim().charAt(0) == '+');
						annot.put(transcriptID, current);
					}
					current.add( split[3], split[4] );
				}
			} else if( information &&  (split[2].equalsIgnoreCase("mRNA") || split[2].equalsIgnoreCase("prediction") || split[2].equalsIgnoreCase("transcript")) ) {
				int idx = split[8].indexOf("ID=")+3;
				int h=split[8].indexOf(";",idx);
				if( h < 0 ) {
					h = split[8].length();
				}
				String id = split[8].substring(idx, h).toUpperCase();
				info.put( id, split[8] );
			}
		}
		r.close();
			
		return annot;
	}
	
	/**
	 * This class represents a gene model.
	 * 
	 * @author Jens Keilwagen
	 *
	 */
	static class Annotation implements Comparable<Annotation> {
		String id, chr;
		boolean forward;
		ArrayList<int[]> cdsParts;
		
		public Annotation( String id, String chr, boolean forward ) {
			this.id = id;
			this.chr = chr;
			this.forward = forward;
			cdsParts = new ArrayList<int[]>();
		}
		
		public void add( String start, String end ) {
			cdsParts.add( new int[]{ Integer.parseInt( start ), Integer.parseInt( end )	} );
		}
		
		static int getLength( int[] interval ) {
			return interval[1]-interval[0]+1;
		}
		
		public int getLength() {
			int l = 0;
			for( int i = 0; i < cdsParts.size(); i++ ) {
				int[] interval = cdsParts.get(i);
				l += getLength(interval);
			}
			return l;
		}
		
		public int getMin() {
			int min = Integer.MAX_VALUE;
			for( int i = 0; i < cdsParts.size(); i++ ) {
				int[] interval = cdsParts.get(i);
				if( interval[0] < min ) {
					min=interval[0];
				}
			}
			return min;
		}
		
		public int getMax() {
			int max = Integer.MIN_VALUE;
			for( int i = 0; i < cdsParts.size(); i++ ) {
				int[] interval = cdsParts.get(i);
				if( interval[1] > max ) {
					max=interval[1];
				}
			}
			return max;
		}
				
		public void set( BitSet region, int start ) {
			for( int i = 0; i < cdsParts.size(); i++ ) {
				int[] interval = cdsParts.get(i);
				region.set(interval[0]-start, interval[1]-start+1);
			}
		}
		
		public int compareTo( Annotation a ) {
			int d = chr.compareTo(a.chr);
			if( d == 0 ) {
				if( forward == a.forward ) {
					d = getMin() - a.getMin();
				} else {
					d = forward==true?-1:1;
				}
			}
			return d;
		}
		
		public String toString( boolean full ) {
			StringBuffer sb = new StringBuffer( id + "\t" + chr + "\t" + forward + "\t" + getMin() + "\t" + getMax() );
			if( full ) {
				for( int i = 0; i < cdsParts.size(); i++ ) {
					sb.append("\n" + i + "\t" + Arrays.toString(cdsParts.get(i)) );
				}
			}
			sb.append("\n"+getLength());
			return sb.toString();
		}
		
		public String toString() {
			return toString(false);
		}
	}
	
	/**
	 * This {@link PseudoAnnotation} allows for filtering the {@link Annotation}.
	 * 
	 * @author Jens Keilwagen
	 */
	static class PseudoAnnotation extends Annotation {
		
		int v;
		
		public PseudoAnnotation() {
			super(null, null, true);
		}
		
		void set(String chr, boolean forward, int v) {
			this.v = v;
			this.chr = chr;
			this.forward = forward;
		}
		
		public int getMin() {
			return v;
		}
	}
	
	/**
	 * This method determines the best matching hit for each gene model from the prediction.
	 * 
	 * @param gene a {@link HashMap} containing the information for the genes/transcripts
	 * @param sep the separator for <code>prediction</code> to cut the suffix
	 * @param alias the {@link HashMap} for the alias names or <code>null</code>
	 * @param truth the true annotation
	 * @param prediction the predicted annotation
	 * @param file the output file
	 * @param protocol for writing messages
	 * 
	 * @throws IOException if there is any problem with the files
	 */
	private static void bestHit( HashMap<String,String[]> gene, String sep, HashMap<String,String[]> alias, HashMap<String,Annotation> truth, HashMap<String,Annotation> prediction, File file, Protocol protocol ) throws IOException {
		Annotation[] a = new Annotation[truth.size()];
		truth.values().toArray(a);
		int[] end = getEnds(a);
		
		int[] counts = new int[9];
		int[] bestCounts = new int[9];
		
		HashSet<Integer> check = new HashSet<Integer>();
		IntList best = new IntList();

		PseudoAnnotation ps = new PseudoAnnotation();
		
		BufferedWriter w = new BufferedWriter( new FileWriter( file ) );
		String[] id = { "?", "?"};

		String transcript, predictionID;
		w.append("#gene\ttranscript\t#exons");//
		w.append("\tprediction\t#predicted exons\tchr\tstrand\tstart\tstop\tnumber of best hits\tf1\tinfo: id,annotated exons,tp,fn,fp,perfect exons,missed exons,superfluous exons,max exon splice error,perfect start,perfect end\tinfo");
		w.newLine();
		String[] array = prediction.keySet().toArray(new String[0]);
		Arrays.sort(array);
		
		HashSet<String> used = new HashSet<String>();
		for( int i = 0; i < array.length; i++ ) {
			predictionID = array[i];
			int index = predictionID.lastIndexOf(sep);
			transcript = index>0 ? predictionID.substring(0,index) : predictionID;
			if( alias != null ) {
				//System.out.print(id);
				transcript = alias.get(transcript)[0];
				//System.out.println(" -> " + id);
			}
			if( gene != null ) {
				id = gene.get(transcript);
			}
			if( id != null ) {
				w.append( id[0] + "\t" + transcript + "\t" + id[1] + "\t" + predictionID);
				//find overlapping transcripts
				Annotation test = prediction.get(predictionID);
				check.clear();
				if( test == null ) {
					w.append( "\tNA\tNA\tNA\tNA\tNA" );
				}
				else {
					w.append( "\t" + test.cdsParts.size() + "\t" + test.chr + "\t" + (test.forward?"+":"-") + "\t" + test.getMin() + "\t" + test.getMax() );
					
					getCandidates(test, ps, a, end, check);
					//getCandidates(test, ps, b, endDis, checkDis);
				}
				
				if( check.size()>0 ) {
					Iterator<Integer> it = check.iterator();
					while( it.hasNext() ) {
						used.add(a[it.next()].id);
					}
				}
				
				//select best
				double x = getBest(test, a, check, counts, bestCounts, best, false );
				//double y = getBest(test, b, checkDis, counts, bestCountsDis, bestDis, false );
				
				if( best.length() > 0 ) {

					double f1;
					IntList idx;
					Annotation[] c;
					f1 = x;
					idx = best;
					c = a;
					w.append( "\t" + idx.length() + "\t" + f1);
					
					for( int j = 0; j < idx.length(); j++ ) {
						int k = idx.get(j);
						w.append( (j==0?"\t":";") );
						compare( test, c[k], counts, true );
						/*if( counts[6] < Integer.MAX_VALUE && counts[1]+counts[2] < counts[6] ) {
							System.out.println(Arrays.toString(counts));
							System.out.println( (counts[1]+counts[2]) + " < " + counts[6] );
							System.out.println();
							System.out.println(test.toString(true));
							System.out.println();
							System.out.println(c[k].toString(true));
							System.out.println();
							System.exit(1);
							//TODO
						}/**/
						String s = Arrays.toString(counts).replaceAll(" ", "");
						w.append( c[k].id + "," + c[k].cdsParts.size() + "," + s.substring(1,s.length()-1) );
					}
				} else {
					w.append( "\t0\tNA\tNA" );
				}
				w.append("\t"+info.get(predictionID));
				w.newLine();
			} else {
				protocol.append(transcript+"\n");
			}
		}
		w.close();
		
		protocol.append("\n");
		int count = 0;
		for( int i = 0; i < a.length; i++ ) {
			if( !used.contains( a[i].id ) ) {
				count++;
				protocol.append(a[i].id+ "\n");
			}
		}
		protocol.append("unused: " + count +"\n\n");
	}
	
	private static int[] getEnds( Annotation[] a ) {
		Arrays.sort(a);
		int[] end = new int[a.length];
		for( int j = 0; j < a.length; j++ ) {
			if( j == 0 || !(a[j].chr.equals(a[j-1].chr) && (a[j].forward==a[j-1].forward)) ) {
				end[j] = a[j].getMax();
			} else {
				end[j] = Math.max(end[j-1],a[j].getMax());
			}
		}
		return end;
	}
	
	private static void getCandidates( Annotation test, PseudoAnnotation ps, Annotation[] a, int[] end, HashSet<Integer> check ) {
		for( int k = 0; k < test.cdsParts.size(); k++ ) {
			int[] x = test.cdsParts.get(k);
			ps.set( test.chr, test.forward, x[1] );
			int stop = getIndex(Arrays.binarySearch(a, ps));
			int start = stop-1, min = x[0];
			while( start >= 0 && end[start] >= min && (a[start].chr.equals(test.chr) && (a[start].forward==test.forward)) ) {
				check.add(start);
				start--;
			}
		}
	}
	
	private static double getBest( Annotation test, Annotation[] a, HashSet<Integer> check, int[] counts, int[] bestCounts, IntList best, boolean exon ) {
		best.clear();
		Arrays.fill(bestCounts, 0);
		double d = 0;
		Iterator<Integer> it = check.iterator();
		while( it.hasNext() ) {
			int k = it.next();
			compare(test, a[k], counts, exon );
			double h = getF1(counts);
			if( h > d ) {
				d = h;
				best.clear();
				best.add(k);
				System.arraycopy(counts, 0, bestCounts, 0, 3);
			} else if( d>0 && h == d ) {
				best.add(k);
			}
		}
		return d;
	}
	
	private static double getF1( int[] counts ) {
		return 2d*counts[0]/(double)(2d*counts[0]+counts[1]+counts[2]);
	}
	
	static int getIndex( int insert ) {
		int end;
		if( insert < 0 ) {
			end = -(insert+1);
		} else {
			end = insert+1;
		}
		return end;
	}
	
	private static boolean compare( Annotation predictedAnnot, Annotation trueAnnot, int[] counts, boolean exon ) {
		Arrays.fill(counts, 0);
		boolean chrAndStrandCorrect = (predictedAnnot.chr.equals(trueAnnot.chr) && (predictedAnnot.forward==trueAnnot.forward));
		if( chrAndStrandCorrect ) {
			int start = Math.min( predictedAnnot.getMin(), trueAnnot.getMin() )-1;
			int end = Math.max( predictedAnnot.getMax(), trueAnnot.getMax() )+1;
			BitSet predictedRegion = new BitSet(end-start+1);
			BitSet trueRegion = new BitSet(end-start+1);
			predictedAnnot.set( predictedRegion, start );
			trueAnnot.set( trueRegion, start );
			
			counts[0] = count( predictedRegion, trueRegion, false); //tp
			counts[1] = count( predictedRegion, trueRegion, true); //fn
			counts[2] = count( trueRegion, predictedRegion, true); //fp
			counts[7] = (trueAnnot.forward ? trueRegion.nextSetBit(0) == predictedRegion.nextSetBit(0) : trueRegion.previousSetBit(end-start) == predictedRegion.previousSetBit(end-start)  ) ? 1 : 0;
			counts[8] = (!trueAnnot.forward ? trueRegion.nextSetBit(0) == predictedRegion.nextSetBit(0) : trueRegion.previousSetBit(end-start) == predictedRegion.previousSetBit(end-start)  ) ? 1 : 0;
			if( exon ) {
				//perfect exon
				for( int i = 0; i < trueAnnot.cdsParts.size(); i++ ) {
					int[] border = trueAnnot.cdsParts.get(i);
					int s = border[0]-1, e = border[1]+2;
					while( s <= e && predictedRegion.get(s-start) == trueRegion.get(s-start) ) {
						s++;
					}
					if( s > e ) {
						counts[3]++;
					}
				}
				
				//missed exon
				for( int i = 0; i < trueAnnot.cdsParts.size(); i++ ) {
					int[] border = trueAnnot.cdsParts.get(i);
					int p = predictedRegion.nextSetBit(border[0]-start);
					if( p==-1 || p > border[1]-start ) {
						counts[4]++;
					}
				}
				
				for( int i = 0; i < predictedAnnot.cdsParts.size(); i++ ) {
					int[] border = predictedAnnot.cdsParts.get(i);
					int p = trueRegion.nextSetBit(border[0]-start);
					if( p==-1 || p > border[1]-start ) {
						counts[5]++;//superfluous exon
					} else {
						int dist;
						
						//left difference
						if( p > border[0]-start ) {
							dist= p - (border[0]-start);
						} else {
							dist= (border[0]-start) - (trueRegion.previousClearBit(border[0]-start)+1);
						}
						if( dist > counts[6] ) {
							counts[6] = dist;
						}
						
						//right difference
						p = trueRegion.previousSetBit(border[1]-start);
						if( p < border[1]-start ) {
							dist = (border[1]-start)-p;
						} else {
							dist= trueRegion.nextClearBit(border[1]-start)-1 - (border[1]-start);
						}
						if( dist > counts[6] ) {
							counts[6] = dist;
						}
						/*
						if( counts[4]==0 ) {
							int b0 = trueRegion.get(border[0]-start) ? trueRegion.previousClearBit(border[0]-start)+1 : trueRegion.nextSetBit(border[0]-start);
							int b1 = trueRegion.get(border[1]-start) ? trueRegion.nextClearBit(border[1]-start)-1 : trueRegion.previousSetBit(border[1]-start);
							int h1, h2;
							//System.out.println(Arrays.toString(border) + "\t" + (border[0]-start) + ", " + (border[1]-start) + "\t" + b0 + ", " + b1 +"\t" + trueRegion.nextClearBit(b0) );
							if( trueRegion.nextClearBit(b0) > b1 //continuous stretch b0 ... b1
								&& (b1 < border[1]-start || ((h1=predictedRegion.nextSetBit(border[1]-start+1)) == -1 || h1 > b1))
								&& (b0 > border[0]-start || ((h2=predictedRegion.nextSetBit(b0)) == -1 || h2 == border[0]-start))
							) {
								int v = Math.max( Math.abs(border[0]-start-b0), Math.abs(border[1]-start-b1) );
								//System.out.println( v +"\t" + counts[6] );
								if( v > counts[6] ) {
									counts[6] = v;
								}
							} else {
								counts[6] = Integer.MAX_VALUE;
							}
						} else {
							counts[6] = Integer.MAX_VALUE;
						}*/
					}
				}
				//System.out.println();
			}
		} else {
			counts[0] = 0; //tp
			counts[1] = trueAnnot.getLength(); //fn
			counts[2] = predictedAnnot.getLength();  //fp
			if( exon ) {
				counts[3] = 0;//perfectly predicted exons
				counts[4] = predictedAnnot.cdsParts.size();//predicted exons with no overlap
				counts[5] = trueAnnot.cdsParts.size();//annotated exons with no overlap
				counts[6] = Integer.MAX_VALUE;
			}
			counts[7] = counts[8] = 0;
		}
		return chrAndStrandCorrect;
	}
	
	private static int count( BitSet b1, BitSet b2, boolean flipFirst ) {
		BitSet h = (BitSet) b1.clone();
		if( flipFirst ) {
			h.flip(0,h.size());
		}
		h.and(b2);
		return h.cardinality();
	}

	@Override
	public ParameterSet getToolParameters() {
		try{
			return new SimpleParameterSet(
					new FileParameter( "prediction", "The predicted annotation", "gff", true ),
					new FileParameter( "annotation", "The true annotation", "gff", true ),
					new FileParameter( "assignment", "the transcript info for the reference of the prediction", "tabular", false )
			);
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public String getToolName() {
		return "Compare transcripts";
	}

	@Override
	public String getToolVersion() {
		return GeMoMa.version;
	}

	@Override
	public String getShortName() {
		return "CompareTranscripts";
	}

	@Override
	public String getDescription() {
		return "This tool helps to compare transcripts from different annotations (e.g. prediction vs. given annotation)";
	}

	@Override
	public String getHelpText() {
		return "This tool compares a predicted annotation with a given annotation in terms of F1 measure. If the F1 measure is 1 both annotations are in perfect agreement for this transcript. The smaller the value is the low is the agreement. If it is NA then there is no overlapping annotation.";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "tabular", "comparison")
		};
	}
}
