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
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import javax.management.RuntimeErrorException;

import projects.gemoma.Synteny2.Prediction;

/**
 * This class extracts the information need to run GeMoMa.
 * 
 * @author Jens Keilwagen
 * 
 * @see GeMoMa
 */
public class GFFFilter {
	
	/**
	 * 
	 * @param args
	 * 0 GFF
	 * 1 tag
	 * 
	 * @throws Exception
	 */
	public static void main ( String[] args ) throws Exception {
		
		String line, t;
		String[] split;
		int idx, h;
				
		ArrayList<Prediction> pred = new ArrayList<Prediction>();
		Prediction current = null;

		//read genes
		BufferedReader r = new BufferedReader( new FileReader(args[0]) );
		while( (line=r.readLine()) != null ) {
			if( line.length() == 0 || line.startsWith("##") ) continue;
			
			split = line.split("\t");
			
			t = split[2];
			if( t.equalsIgnoreCase( args[1] ) ) {
				current = new Prediction(split);
				pred.add( current );
			} else {
				current.addCDS( line );
			}
			
		}
		r.close();
		int all = pred.size();
		
		//initial filter
		double th = 0.75;
		int anz = 0, out = 0;
		for( int i = pred.size()-1; i >= 0; i-- ){
			Prediction p = pred.get(i);
			if( p.getRelScore()>=th && p.hash.get("start").charAt(0)=='M' && p.hash.get("stop").charAt(0)=='*' ) {
				//good
				anz++;
			} else {
				pred.remove(i);
				out++;
			}
		}
		//sort
		Collections.sort(pred);

		//filter overlapping and quality, and write
		BufferedWriter w = new BufferedWriter( new FileWriter(args[2]) );
		int i = 0, anz2=0, anz3=0;
		Prediction next = pred.get(i);
		while( i < pred.size() ) {
			current = next;
			int alt=i, best=i, end = current.end;
			i++;	

			while( i < pred.size() && (next = pred.get(i)).contigStrand.equals(current.contigStrand)
					//&& (end-Integer.parseInt(next.split[3])+1d)/(end-start+1d) > 0.1
					&& end>Integer.parseInt(next.split[3])
			) {
				end = Integer.parseInt(next.split[4]);
				
				if( next.score > current.score ) {
					best = i;
					current = next;
				}
				i++;
			}
			
			int p=create( pred, alt, best, i, w );
			anz2++;
			anz3+=p;
		}
		w.close();
		System.out.println();
		System.out.println(anz2);
		System.out.println( "("+ anz3 + " + " + (anz-anz3) + " + " + out +") / " + all );
	}
	
	private static boolean extract( Prediction p ) {
		String s = p.hash.get("tae");
		double tae = s.equals("NA") ? Double.NaN : Double.parseDouble(s);
		s =p. hash.get("tde");
		double tde = s.equals("NA") ? Double.NaN : Double.parseDouble(s);
		s = p.hash.get("tie");
		double tie = s.equals("NA") ? Double.NaN : Double.parseDouble(s);
		String start = p.hash.get("start");
		String stop = p.hash.get("stop");
		
		//System.out.println(p.hash.get("ID")+"\t" + p.length + "\t" + (p.length%3) + "\t"+p.hash.get("score")+"\t"+tae+"\t"+tde+"\t"+tie+ "\t" + (Integer.parseInt(p.hash.get("score"))/(p.length/3)>1) );		
		
		return p.getRelScore()>0.75 && (start.charAt(0)=='M' && stop.charAt(0)=='*');// || (!Double.isNaN(tae) && tie==1);//TODO 
	}
	
	static int create( ArrayList<Prediction> list, int start, int best, int end, BufferedWriter w ) throws Exception {
		Prediction current = list.get(best), n;
		current.write(w);
		current.discard=true;
		String s;
		//= current.hash.get("tie");
		//boolean hasTie = s != null && (s.equals("NA") || s.equals("0.0"));
		
		int pred=1;
		ArrayList<Prediction> used = new ArrayList<Prediction>();
		used.add( current );
		
		
		//TODO same gene?
		
		if( end - start > 1 ) {
			int count  = 0, idx=-1;
			for( int i = start; i < end; i++ ) {
				n = list.get(i);
				if( !n.discard ) {
					if( n.end < current.start || current.end < n.start ) {
						if( count==0 || list.get(idx).score<n.score ) {
							idx=i;
						}
						count++;
					} else {
						s = n.hash.get("tie");
						if( /*!hasTie ||*/ s != null ) {
							double tie = s == null ? 0 : (s.equals("NA") ? Double.NaN : Double.parseDouble(s));
							if( /*!hasTie ||*/ tie >= 1 ) 
							{
								Prediction x = used.get(0); 
								double cb = x.commonBorders(n);
								double min=2d*Math.min(n.cds.size(),x.cds.size());
								double max=2d*Math.max(n.cds.size(),x.cds.size());
								
								int j = 1; 
								
								if( cb/min>0.75 && cb/max < 1 ) {//hinreichende Überlappung, aber keine perfekte überlappung
									while( j < used.size() ) {
										x = used.get(0); 
										cb = x.commonBorders(n);
										//min=2d*Math.min(n.cds.size(),x.cds.size());
										max=2d*Math.max(n.cds.size(),x.cds.size());
										
										if( cb/max==1 )  {
											break;
										}
										j++;
									}
								}
								if( j == used.size() ) {
									n.write(w);
									used.add(n);
									pred++;
								}
							}
						}
						n.discard=true;
					}
				}
			}
			if( count>0  ) {
				pred+=create(list, start, idx, end, w);
			}
		}
		return pred;
	}
	
	
	static class Prediction implements Comparable<Prediction>{
		boolean discard = false;
		String[] split;
		ArrayList<String> cds;
		HashMap<String,String> hash;
		int length, start, end, score;
		String contigStrand;
		
		public Prediction( String[] split ) {
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
		
		public void write( BufferedWriter w ) throws IOException {
			for( int i = 0; i < split.length; i++ ) {
				w.append( (i==0?"":"\t") + split[i] );
			}
			w.newLine();
			for( int i = 0; i < cds.size(); i++ ) {
				w.append( cds.get(i) );
				w.newLine();
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
					
					//System.out.println(c);
				}
			}/* else {
				for( int i = 0; i < cds.size(); i++ ) {
					c = cds.get(i);
					System.out.println(c);
				}
			}*/
			//System.out.println();
			for( int i = 0; i < o.cds.size(); i++ ) {
				c = o.cds.get(i);
				split = c.split("\t");
				
				if( s.contains(split[3]) ) {
					anz++;
				}
				if( e.contains(split[4]) ) {
					anz++;
				}
				
				//System.out.println(c);
			}
			return anz;
		}
	}
}