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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import de.jstacs.utils.IntList;


/**
 * This class allows for comparing GFFs.
 * It determines for each predicted gene model the best matching annotated gene model by evaluating the F1 measure.
 * 
 * @author Jens Keilwagen
 */
public class CompareTranscripts {

	/**
	 * Just the main.
	 * 
	 * @param args
	 * 0 ... transcript info (query organism)<br/>
	 * 1 ... transcript info (target organism)<br/>
	 * 2 ... true annotation (gff)<br/>
	 * 3 ... predicted annotation (gff)<br/>
	 * 4 ... alias file for prediction (if not found <code>null</code>)<br/>
	 * 5 ... output file<br/>
	 * 
	 * @throws Exception if something went wrong
	 */
	public static void main(String[] args) throws Exception {
		HashMap<String, String[]> gene = Tools.getAlias(args[0], 1, 0, 2);
		//System.out.println(gene.toString().substring(0,1000) + "...");
		HashMap<String, String[]> alias = Tools.getAlias(args[4], 0, 1, -1);
		System.out.println(alias==null);
		/*for( int i = 0; i < args.length; i++ ) {
			System.out.println(i+"\t" + args[i]);
		}*/
		HashMap<String,Annotation> discarded = new HashMap<String, Annotation>();
		HashMap<String,Annotation> truth = readGFF( args[2], "CDS", false, args[1], discarded, ";" );
		HashMap<String,Annotation> prediction;
		if( alias == null ) {
			prediction = readGFF( args[3], "CDS", false, null, null, ";" );//_R" );//standard
		} else {
			prediction = readGFF( args[3], "coding_exon", true, null, null, "-R" );//genBlast
		}
		System.out.println(truth.size() + " vs. " + prediction.size() );
		
		bestHit(gene, alias==null?"_R":"-R", alias, truth, discarded, prediction, args[5]);
		System.out.println("problem: " + Arrays.toString(problem));
	}
	
	private static int problem[] = new int[2];
	
	private static String par = "Parent=";
	
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
	private static HashMap<String,Annotation> readGFF( String input, String type, boolean skip, String selected, HashMap<String,Annotation> discarded, String... sep ) throws Exception {
		BufferedReader r;
		String line;
		
		HashSet<String> sel, del = new HashSet<String>();
		if( selected == null || !(new File( selected )).exists() ) {
			sel = null;
		} else {
			sel = new HashSet<String>();
			r = new BufferedReader( new FileReader(selected) );
			while( (line=r.readLine()) != null ) {
				sel.add( line.split("\t")[1] );
			}
			r.close();
			//System.out.println(sel);
		}
		System.out.println("read " + (selected==null) );
		
		HashMap<String,Annotation> hAnnot, annot = new HashMap<String,Annotation>();
		Annotation current;
		r = new BufferedReader( new FileReader(input) );
		while( (line=r.readLine()) != null ) {
			if( line.length() != 0 && !line.startsWith("#") ) {
				int idx = line.indexOf('\t')+1;
				idx = line.indexOf('\t',idx)+1;
				int end = line.indexOf('\t',idx); 
				String t = line.substring(idx,end);
				if( t.equalsIgnoreCase(type) ) {
					String[] split = line.split("\t");
					idx = split[8].indexOf(par)+par.length();
					int h=Integer.MAX_VALUE, f;
					for( int i = 0; i < sep.length; i++ ) {
						f = split[8].indexOf(sep[i],idx);
						if( f >= 0 ) {
							h = Math.min(h, f);
						}
					}
					//System.out.println(idx + "\t" + h + "\t" + split[8]);
					if( !skip || h> 0 ) {
						String transcriptID = split[8].substring(idx, h<Integer.MAX_VALUE?h:split[8].length() );

						if( sel == null || sel.contains(transcriptID) ) {
							hAnnot = annot;
						} else if( discarded != null ) {
							hAnnot = discarded;
						} else {
							hAnnot = null;
						}
						
						if( hAnnot != null ) {
							current = hAnnot.get(transcriptID);
							if( current == null ) {
								current = new Annotation( transcriptID, split[0], split[6].trim().charAt(0) == '+');
								hAnnot.put(transcriptID, current);
								if( hAnnot == discarded ) {
									del.add(transcriptID);
								}
							}
							current.add( split[3], split[4] );
						}
					}
				}
			}
		}
		r.close();
		
		if( del.size() > 0 ) {
			System.out.println(del.size());
			//System.out.println(del);
		}
		
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
	 * @param fileName the output file name
	 * 
	 * @throws IOException if there is any problem with the files
	 */
	private static void bestHit( HashMap<String,String[]> gene, String sep, HashMap<String,String[]> alias, HashMap<String,Annotation> truth, HashMap<String,Annotation> discarded, HashMap<String,Annotation> prediction, String fileName ) throws IOException {
		Annotation[] a = new Annotation[truth.size()];
		truth.values().toArray(a);
		int[] end = getEnds(a);
		
		Annotation[] b = new Annotation[discarded.size()];
		discarded.values().toArray(b);
		int[] endDis = getEnds(b);
				
		System.out.println(a.length + "\t" + b.length);

		int[] counts = new int[9];
		int[] bestCounts = new int[9];
		
		HashSet<Integer> check = new HashSet<Integer>();
		IntList best = new IntList();

		HashSet<Integer> checkDis = new HashSet<Integer>();
		int[] bestCountsDis = new int[6];
		IntList bestDis = new IntList();

		PseudoAnnotation ps = new PseudoAnnotation();
		
		BufferedWriter w = new BufferedWriter( new FileWriter( fileName ) );
		String[] id;
		String transcript, predictionID;
		w.append("#gene\ttranscript\t#exons");//
		w.append("\tprediction\t# predicted exons\tchr\tstrand\tstart\tstop\tnumber of best hits\tf1\tinfo: id,annotated exons,tp,fn,fp,perfected exons,missed exon,superfluous exons,perfect start,perfect end");
		w.newLine();
		String[] array = prediction.keySet().toArray(new String[0]);
		Arrays.sort(array);
		for( int i = 0; i < array.length; i++ ) {
			predictionID = array[i];
			int index = predictionID.lastIndexOf(sep);
			transcript = index>0 ? predictionID.substring(0,index) : predictionID;
			if( alias != null ) {
				//System.out.print(id);
				transcript = alias.get(transcript)[0];
				//System.out.println(" -> " + id);
			}
			id = gene.get(transcript);
			if( id != null ) {
				w.append( id[0] + "\t" + transcript + "\t" + id[1] + "\t" + predictionID);
				//find overlapping transcripts
				Annotation test = prediction.get(predictionID);
				check.clear();
				checkDis.clear();
				if( test == null ) {
					w.append( "\tNA\tNA\tNA\tNA\tNA" );
				}
				else {
					w.append( "\t" + test.cdsParts.size() + "\t" + test.chr + "\t" + (test.forward?"+":"-") + "\t" + test.getMin() + "\t" + test.getMax() );
					
					getCandidates(test, ps, a, end, check);
					getCandidates(test, ps, b, endDis, checkDis);
				}
				
				//select best
				double x = getBest(test, a, check, counts, bestCounts, best, false );
				double y = getBest(test, b, checkDis, counts, bestCountsDis, bestDis, false );
				
				if( best.length() > 0 || bestDis.length() > 0 ) {
					/*
					if( test.id.startsWith("AT1G01020")
							|| (bestCounts[0]+bestCounts[1]) % 3 > 0
							||  (bestCounts[0]+bestCounts[2]) % 3 > 0) {
						if( (bestCounts[0]+bestCounts[1]) % 3 > 0 ) {
							problem[0]++;
						}
						if( (bestCounts[0]+bestCounts[2]) % 3 > 0 ) {
							problem[1]++;
						}
						
						System.out.println(Arrays.toString(bestCounts));
						System.out.println();
						System.out.println( (bestCounts[0]+bestCounts[2]) % 3 );
						System.out.println("prediction:\n"+test.toString(true));
						System.out.println();
						System.out.println( (bestCounts[0]+bestCounts[1]) % 3 );
						System.out.println("truth:\n"+a[best.get(0)].toString(true));
						System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
						System.exit(1);
					}/**/
					
					double f1;
					IntList idx;
					Annotation[] c;
					if( best.length() > 0 ) {
						f1 = x;
						idx = best;
						c = a;
					} else {
						f1 = -100+y;
						idx = bestDis;
						c = b;
					}
					w.append( "\t" + idx.length() + "\t" + f1);
					
					for( int j = 0; j < idx.length(); j++ ) {
						int k = idx.get(j);
						w.append( (j==0?"\t":";") );
						compare( test, c[k], counts, true );
						if( counts[6] < Integer.MAX_VALUE && counts[1]+counts[2] < counts[6] ) {
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
				w.newLine();
			}
		}
		w.close();
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
						counts[5]++;
					} else {
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
						}
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
}
