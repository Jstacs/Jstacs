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
import java.util.Arrays;

/**
 * A simple tool for combining coverage files that might be created in parallel using a compute cluster.
 * 
 * @author Jens Keilwagen
 */
public class CombineCoverageFiles {

	public static void main(String[] args) throws IOException, CloneNotSupportedException {
		System.out.println("files: " + (args.length-1));
		
		BufferedWriter w = new BufferedWriter( new FileWriter(args[0]) );
		//read files
		BufferedReader[] r = new BufferedReader[args.length-1];
		BedEntry[] bed = new BedEntry[args.length-1];
		String line;
		for( int i = 0; i < r.length; i++ ) {
			System.out.println(i + "\t" + args[i+1]);
			r[i] = new BufferedReader( new FileReader( args[i+1] ) );
			//skip header
			line=r[i].readLine();
			if( i==0 ) {
				w.append(line);
				w.newLine();
			}
			bed[i]=new BedEntry(r[i].readLine()); 
		}
		
		BedEntry current = null, next;
		boolean[] sameStart = new boolean[r.length];
		int same=0;
		do {
			//find next
			next=null;
			same=0;
			Arrays.fill(sameStart, false);
			for( int i = 0; i < bed.length; i++ ) {
				if( bed[i]!=null ) {
					if( next==null ) {
						next=bed[i];
						sameStart[i]=true;
						same++;
					} else {
						int comp=next.compareTo(bed[i]);
						if( comp>=0 ) {
							if( comp>0 ) {
								Arrays.fill(sameStart, false);
								next=bed[i];
								same=0;
							}
							sameStart[i]=true;
							same++;
						}
					}
				}
			}
			
			if( same>0 ) {
				//sum
				BedEntry now = next.clone();
				boolean first=true;
				for( int i = 0; i < bed.length; i++ ) {
					if( sameStart[i] ) {
						if( !first ) {
							if( bed[i].end < now.end ) {
								now.end = bed[i].end;
							}
							now.cov += bed[i].cov;
						} else {
							first=false;
						}
					} else {
						if( bed[i] != null ) {
							if( bed[i].chr.equals(now.chr) && bed[i].start<now.end ) {
								now.end=bed[i].start;
							}
						}
					}
				}
				
				
				//update array
				for( int i = 0; i < sameStart.length; i++ ) {
					if( sameStart[i] ) {
						if( bed[i].end == now.end ) {
							line = r[i].readLine();
							if( line != null ) {
								bed[i] = new BedEntry(line);
							} else {
								bed[i] = null;
							}
						} else {
							//shorten
							bed[i].start=now.end;
						}
					}
				}
				
				//write
				if( current != null ) {
					if( current.chr.equals(now.chr) && current.end == now.start && current.cov==now.cov ) {
						//extend
						current.end = now.end;
					} else {
						//write
						current.write(w);
						current=now;
					}
				} else {
					current=now;
				}
			}
		} while( same>0 );
		if( current !=null ) {
			current.write(w);
		}
		w.close();
	}

	static class BedEntry implements Comparable<BedEntry>, Cloneable {
		String chr;
		int start, end, cov;
		
		BedEntry( String line ) {
			String[] split = line.split("\t");
			chr=split[0];
			start = Integer.parseInt(split[1]);
			end = Integer.parseInt(split[2]);
			cov = Integer.parseInt(split[3]);
		}
		
		public BedEntry clone() throws CloneNotSupportedException {
			return (BedEntry) super.clone();
		}

		@Override
		public int compareTo(BedEntry o) {
			int diff = chr.compareTo(o.chr);
			if( diff == 0 ) {
				diff=Integer.compare(start, o.start);
			}
			return diff;
		}
		
		public void write( BufferedWriter w ) throws IOException {
			w.write(chr + "\t" + start + "\t" + end + "\t" + cov );
			w.newLine();
		}
	}
}