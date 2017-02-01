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
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;

import de.jstacs.DataType;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;

/**
 * This class allows to filter GeMoMa annotation files (GFFs) to obtain the most relevant predictions.
 * 
 * @author Jens Keilwagen
 * 
 * @see GeMoMa
 */
public class GeMoMaAnnotationFilter implements JstacsTool {
	
	private static int MAX;
	
	@Override
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		String tag = parameters.getParameterForName("tag").getValue().toString();
		double relScoTh = (Double) parameters.getParameterForName("relative score filter").getValue();
		boolean complete = (Boolean) parameters.getParameterForName("complete").getValue();
		boolean noTie = (Boolean) parameters.getParameterForName("missing intron evidence filter").getValue();
		double tieTh = (Double) parameters.getParameterForName("intron evidence filter").getValue();
		double cbTh = (Double) parameters.getParameterForName("common border filter").getValue();
		double epTh = (Double) parameters.getParameterForName("evidence percentage filter").getValue();
				
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterAt(6).getValue();
		
		String line, t;
		String[] split;
				
		ArrayList<Prediction> pred = new ArrayList<Prediction>();
		Prediction current = null;

		//read genes
		BufferedReader r;
		MAX = eps.getNumberOfParameters();
		for( int k = 0; k < MAX; k++ ) {
			String fName = ((ParameterSet)eps.getParameterAt(k).getValue()).getParameterAt(0).getValue().toString();
			//System.out.println(fName);

			r = new BufferedReader(new FileReader(fName));
			while( (line=r.readLine()) != null ) {
				if( line.length() == 0 || line.startsWith("##") ) continue;
				
				split = line.split("\t");
				
				t = split[2];
				if( t.equalsIgnoreCase( tag ) ) {
					current = new Prediction(split, k);
					pred.add( current );
				} else {
					current.addCDS( line );
				}
				
			}
			r.close();
			//System.out.println( pred.size() + "\t" + k + "\t" +current.index);
		}
		protocol.append( "all: " + pred.size() + "\n" );
		
		//initial filter
		int filtered = 0, clustered=0, transcripts=0;
		HashMap<Integer, int[]> counts = new HashMap<Integer, int[]>();
		for( int i = pred.size()-1; i >= 0; i-- ){
			Prediction p = pred.get(i);
			if( p.getRelScore()>=relScoTh && (!complete || (p.hash.get("start").charAt(0)=='M' && p.hash.get("stop").charAt(0)=='*')) ) {
				filtered++;
			} else {
				pred.remove(i);
			}
			int[] c = counts.get(p.index);
			if( c == null ) {
				c = new int[1];
				counts.put(p.index, c);
			}
			c[0]++;
		}
		protocol.append( "filtered: " + filtered + "\n" );
		
		File out = GeMoMa.createTempFile("GAF-filtered");
		if( filtered > 0 ) {
			Iterator<Entry<Integer,int[]>> it = counts.entrySet().iterator();
			/*while( it.hasNext() ) {
				Entry<Integer,int[]> e = it.next();
				System.out.println(e.getKey() + "\t" + e.getValue()[0]);
			}/**/
			
			//sort
			Collections.sort(pred);
	
			//cluster predictions and filter them
			
			BufferedWriter w = new BufferedWriter( new FileWriter(out) );
			w.append("##gff-version 3");
			w.newLine();
			int i = 0;
			Prediction next = pred.get(i);
			while( i < pred.size() ) {
				current = next;
				int alt=i, best=i, end = current.end;
				i++;	
	
				while( i < pred.size() && (next = pred.get(i)).contigStrand.equals(current.contigStrand)
						//&& (end-Integer.parseInt(next.split[3])+1d)/(end-start+1d) > 0.1
						&& end>Integer.parseInt(next.split[3])
				) {
					end = Math.max(end,Integer.parseInt(next.split[4]));
					
					if( next.score > current.score ) {
						best = i;
						current = next;
					}
					i++;
				}
				//Kalkuel: next=pre.get(i), best = argmax_{j \in [alt,i-1]} pred.get(j).score 
				
				int p=create( pred, alt, best, i, w, noTie, tieTh, cbTh, epTh );
				clustered++;
				transcripts+=p;
			}
			w.close();
		}
				
		protocol.append( "clustered: " + clustered + "\n" );
		protocol.append( "transcripts: " + transcripts + "\n" );
		protocol.append( "genes: " + gene + "\n\n" );
		protocol.append( "no tie: " + GeMoMaAnnotationFilter.noTie + "\n" );
		protocol.append( "tie=1: " + tie1 + "\n" );
		
		
		return new ToolResult("", "", null, new ResultSet(new TextResult("filtered predictions", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "gff", getToolName(), null, true)), parameters, getToolName(), new Date());
	}
	
	static int gene=0, noTie=0, tie1;
	
	static int create( ArrayList<Prediction> list, int start, int best, int end, BufferedWriter w, boolean noTie, double tieTh, double cbTh, double epTh ) throws Exception {
		Prediction current = list.get(best);
		ArrayList<Prediction> used = new ArrayList<Prediction>();
		used.add( current );
		current.discard = true;	

		int count = 0, idx=-1;
		if( end - start > 1 ) {
			for( int i = start; i < end; i++ ) {
				Prediction n = list.get(i);
				if( !n.discard ) {
					if( n.end < current.start || current.end < n.start ) {//no overlap
						if( idx<0 || list.get(idx).score<n.score ) {
							idx=i;
						}
						count++;
					} else { //overlapping
						Prediction x = used.get(0); 
						double cb = x.commonBorders(n);
						double min=2d*Math.min(n.cds.size(),x.cds.size());
						double max=2d*Math.max(n.cds.size(),x.cds.size());
				
						int j = 0; 
						if( cb/max == 1d ) { // identical
							x.addAlternative(n);
						} else if( cb/min>=cbTh ) {//hinreichende Überlappung, aber keine perfekte überlappung
							j++;
							while( j < used.size() ) {
								x = used.get(j); 
								cb = x.commonBorders(n);
								max=2d*Math.max(n.cds.size(),x.cds.size());								
								if( cb/max==1d )  {
									if( x.score >= n.score  ) {
										x.addAlternative(n);
									} else {
										n.copyAlternatives(x);
										used.remove(j);
										used.add(n);
									}
									break;
								}
								j++;
							}
						}
						
						if( j == used.size() ) {
							String s = n.hash.get("tie");
							boolean noInfo = s == null || s.equals("NA");
							double tie = noInfo ? Double.NaN : Double.parseDouble(s);
							if( (noInfo && noTie) || tie >= tieTh ) {
								used.add(n);
							}
						}
						n.discard=true;
					}
				}
			}
			
		}
		
		int pred=0;
		maxEvidence=0;
		maxTie=-1;
		complete=0;
		st=Integer.MAX_VALUE;
		en=Integer.MIN_VALUE;
		Prediction n=null;
		
		//sort
		Collections.sort(used, ScoreComparator.DEF);
		
		for( int i = 0; i < used.size(); i++ ) {
			n = used.get(i);
			int cont = n.write(w, epTh);
			pred += cont;
		}
		if( pred>0 ) {
			//write
			w.append(n.split[0] + "\tGAF\tgene\t" + st + "\t" + en  + "\t.\t" + n.split[6] + "\t.\tID=gene_"+gene+";transcripts=" + pred + ";complete="+complete+";maxEvidence="+maxEvidence+";maxTie=" + (maxTie<0?"?":maxTie) );
			w.newLine();
			gene++;
		}
		
		if( count>0  ) {
			pred += create(list, start, idx, end, w, noTie, tieTh, cbTh, epTh);
		}

		return pred;
	}
	
	static int maxEvidence, st, en, complete;
	static double maxTie;
	
	static class ScoreComparator implements Comparator<Prediction> {
		static ScoreComparator DEF = new ScoreComparator();

		@Override
		public int compare(Prediction o1, Prediction o2) {
			return -Integer.compare(o1.score, o2.score);
		}
	}
	
	static class Prediction implements Comparable<Prediction>{
		boolean discard = false;
		String[] split;
		ArrayList<String> cds;
		HashMap<String,String> hash;
		int length, start, end, score;
		String contigStrand;
		HashSet<String> alternative;
		int index;
		boolean[] evidence;
		
		public Prediction( String[] split, int index ) {
			this.index = index;
			
			this.split = split;
			contigStrand = split[0]+split[6];
			start = Integer.parseInt(split[3]);
			end = Integer.parseInt(split[4]);
			
			String s;
			split = split[8].split(";");
			hash = new HashMap<String, String>();
			for( int i = 0; i < split.length; i++ ) {
				int idx = split[i].indexOf('=');
				s = split[i].substring(idx+1);
				if( s.charAt(0)=='?' ) {
					s = "NA";
				}
				hash.put( split[i].substring(0,idx), s );
			}
			score = Integer.parseInt(hash.get("score"));

			cds = new ArrayList<String>();
			length = 0;
			alternative = new HashSet<String>();
		}
		
		void addCDS( String cds ) {
			this.cds.add( cds );
			String[] split = cds.split("\t");
			length+= (Integer.parseInt(split[4])-Integer.parseInt(split[3])+1);
		}

		@Override
		public int compareTo(Prediction o) {
			int res = contigStrand.compareTo(o.contigStrand);
			if( res == 0 ) {
				res = Integer.compare( Integer.parseInt(split[3]), Integer.parseInt(o.split[3]) );
				if( res == 0 ) {
					res = Integer.compare( Integer.parseInt(split[4]), Integer.parseInt(o.split[4]) );
				}
			}
			return res;
		}
		
		public double getRelScore() {
			return score / (length/3d);
		}
		
		public void copyAlternatives( Prediction n ) {
			alternative.clear();
			alternative.addAll(n.alternative);
			
			if( evidence == null ) {
				evidence = new boolean[MAX];
			}
			if( n.evidence!= null ) {
				System.arraycopy(n.evidence, 0, evidence, 0, evidence.length);
			} else {
				Arrays.fill(evidence, false);
			}
			evidence[index]=true;			
			addAlternative(n);
		}
		
		public void addAlternative( Prediction n ) {
			String rg = n.hash.get("ref-gene");
			if( !rg.equals( hash.get("ref-gene") ) ) {
				alternative.add(rg);
				if( index != n.index ) {
					if( evidence == null ) {
						evidence = new boolean[MAX];
						Arrays.fill(evidence, false);
						evidence[index]=true;
					}
					evidence[n.index]=true;
				}
			}
		}
		
		public int write( BufferedWriter w, double epTh ) throws IOException {
			int count=0;
			if( evidence != null ) {
				for( int i = 0; i < evidence.length; i++ ) {
					if( evidence[i] ) {
						count++;
					}
				}
			} else {
				count=1;
			}
			
			if( count/(double) MAX >= epTh ) {
				String t = hash.get("tie");
				if( t == null || t.equals("NA") ) {
					noTie++;
				} else {
					double tie =Double.parseDouble(t);
					if( tie == 1d ) {
						tie1++;
					}
					maxTie = Math.max(tie, maxTie);
				}
				for( int i = 0; i < split.length; i++ ) {
					w.append( (i==0?"":"\t") + split[i] );
				}
				w.write( ";evidence=" + count + ";Parent=gene_"+gene );
			
				maxEvidence = Math.max(count, maxEvidence);
				st=Math.min(st, start);
				en=Math.max(en, end);
				complete += (hash.get("start").charAt(0)=='M' && hash.get("stop").charAt(0)=='*') ? 1 : 0;
				
				if( alternative.size() > 0 ) {
					String[] it = alternative.toArray(new String[0]);
					Arrays.sort(it);
					w.write( ";alternative=\"" + it[0] );
					for( int i = 1; i < it.length; i++ ) {
						w.write("," + it[i] );
					}
					w.write( "\"" );
				}
				w.newLine();
				for( int i = 0; i < cds.size(); i++ ) {
					w.append( cds.get(i) );
					w.newLine();
				}
				return 1;
			} else {
				return 0;
			}
		}
		
		public String toString() {
			String r = "";
			for( int i = 0; i < split.length; i++ ) {
				r += (i==0?"":"\t") + split[i];
			}
			return r;
		}
		
		HashSet<String> s, e;
		
		public double commonBorders(Prediction o) {
			int anz = 0;
			String c;
			String[] split;
			if( s== null ) {
				s = new HashSet<String>();
				e = new HashSet<String>();
				for( int i = 0; i < cds.size(); i++ ) {
					c = cds.get(i);
					split = c.split("\t");
					s.add(split[3]);
					e.add(split[4]);
				}
			}
			for( int i = 0; i < o.cds.size(); i++ ) {
				c = o.cds.get(i);
				split = c.split("\t");
				
				if( s.contains(split[3]) ) {
					anz++;
				}
				if( e.contains(split[4]) ) {
					anz++;
				}
			}
			return anz;
		}
	}


	@Override
	public ParameterSet getToolParameters() {
		try{
			return
				new SimpleParameterSet(
					new SimpleParameter(DataType.STRING,"tag","the tag used to read the GeMoMa annotations",true,"prediction"),
					new SimpleParameter(DataType.DOUBLE,"relative score filter","the initial filter on the relative score (i.e. score devided by length)", true, 0.75 ),
					new SimpleParameter(DataType.BOOLEAN,"complete","only complete predictions (having start and stop codon) pass the initial filter", true, true ),
					new SimpleParameter(DataType.BOOLEAN,"missing intron evidence filter","the filter for single-exon transcripts or if no RNA-seq data is used, decides for overlapping other transcripts whether they should be used (=true) or discarded (=false)", true, false ),
					new SimpleParameter(DataType.DOUBLE,"intron evidence filter","the filter on the intron evidence given by RNA-seq-data for overlapping transcripts", true, new NumberValidator<Double>(0d, 1d), 1d ),
					new SimpleParameter(DataType.DOUBLE,"common border filter","the filter on the common borders of transcripts, the lower the more transcripts will be checked as alternative splice isoforms", true, new NumberValidator<Double>(0d, 1d), 0.75 ),
					new ParameterSetContainer( new ExpandableParameterSet( new SimpleParameterSet(		
							new FileParameter( "gene annotation files", "GFF files containing the gene annotations (predicted by GeMoMa)", "gff",  true )
						), "gene annotation file", "", 1 ) ),
					new SimpleParameter(DataType.DOUBLE,"evidence percentage filter","Each gene annotation file is handle as separate evidence. A prediction is only returned if it is contained at least in this percentage of evidence files.)", true, new NumberValidator<Double>(0d, 1d), 0.5 )
				);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public String getToolName() {
		return "GeMoMa Annotation Filter";
	}

	@Override
	public String getToolVersion() {
		return "1.3.2";
	}

	@Override
	public String getShortName() {
		return "GAF";
	}

	@Override
	public String getDescription() {
		return "filters most reliable gene predictions";
	}

	@Override
	public String getHelpText() {
		return
				"**What it does**\n\nThis tools filters (and combines) gene predictions obtained from GeMoMa. It allows to combine the predictions from multiple reference organisms, but works also using only one reference organism. This tool does not modify the predictions, but filters redundant and low-quality predictions.\n\n"
				+ "**References**\n\nFor more information please visit http://www.jstacs.de/index.php/GeMoMa or contact jens.keilwagen@julius-kuehn.de.\n";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", "filtered predictions"),
		};
	}
}