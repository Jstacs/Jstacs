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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * This class allows to evaluate the extended best reciprocal hit approach.
 * It utilizes the output of {@link CompareTranscripts} as input.
 *  
 * @author Jens Keilwagen
 * 
 * @see CompareTranscripts
 */
public class BestReciprocalHit {

	/**
	 * Just the main.
	 * 
	 * @param args
	 * 0 bestHit-A-on-B<br/>
	 * 1 bestHit-B-on-A<br/>
	 * 2 transcriptExonAssignments A<br/>
	 * 3 out<br/>
	 * 
	 * @throws Exception if something went wrong
	 */
	public static void main(String[] args) throws Exception {	
		//read and save prediction
		HashMap<String, Hit> hash1 = read( args[0] );
		HashMap<String, Hit> hash2 = read( args[1] );

		//process prediction
		BufferedReader r = new BufferedReader( new FileReader( args[2] ) );
		Hit current, h;
		HashSet<String> gene = new HashSet<String>();
		ArrayList<Hit> best = new ArrayList<Hit>();
		double bestF1;
		BufferedWriter w = new BufferedWriter( new FileWriter(args[3]) );
		String line;
		ArrayList<String[]> transcripts = new ArrayList<String[]>();
		HashMap<String, String> transcriptMap = new HashMap<String, String>();
		String[] split;
		while( (line=r.readLine()) != null ) {
			if( line.startsWith("#") ) {
				continue;
			}
			split = line.split("\t");
			transcripts.add(split);
			transcriptMap.put(split[1], split[0]);
		}
		
		for( int k = 0; k < transcripts.size(); k++ ) {
			split = transcripts.get(k);
			w.append( split[0] + "\t" + split[1] );
			
			h = hash1.get(split[1]);
			bestF1=-1;
			best.clear();
			boolean pred=false;
			
			if( h == null ) {
				for( int i = 0; i < 10; i++ ) {
					w.append("\tNA");
				}
			} else {
				w.append( h.toString() );
				for( int i = 0; i < h.matches.length; i++ ) {
					current = hash2.get( h.matches[i] );
					if( current != null ) {
						pred=true;
						if( bestF1 <= current.f1 ) {
							if( bestF1 < current.f1 ) {
								best.clear();
							}
							if( best.size() == 0 ) {
								bestF1=current.f1;
							}
							best.add( current );
						}
					}
				}
			}
			w.append( "\t" + (pred?best.size():-1) + "\t" + (pred?bestF1:"NA") );
			
			if( best.size() != 0 && h!= null && h.f1 >= 0 ) {
				//transcript org. 2
				gene.clear();
				for( int i = 0; i < best.size(); i++ ) {
					gene.add(best.get(i).id);
				}
				w.append( "\t" + toString( gene.toArray() ) );
				
				//transcript org. 1
				gene.clear();
				for( int i = 0; i < best.size(); i++ ) {
					current = best.get(i);
					for( int j = 0; j < current.matches.length; j++ ) {
						gene.add(current.matches[j]);
					}
				}
				Object[] o=gene.toArray();
				w.append( "\t" + toString( o ) );
w.flush();
				//gene org. 1
				gene.clear();
				for( int i = 0; i <o.length; i++ ) {
					gene.add( transcriptMap.get( o[i] ) );
				}
				w.append( "\t" + toString( gene.toArray() ) );
				
			} else {
				for( int i = 0; i < 3; i++ ) {
					w.append( "\tNA");
				}
			}
			w.newLine();
		}
		r.close();
		w.close();
	}
	
	private static HashMap<String, Hit> read( String fileName ) throws Exception {
		//read and save prediction
		HashMap<String, Hit> hash = new HashMap<String, Hit>();
		BufferedReader r = new BufferedReader( new FileReader( fileName ) );
		String line;
		while( (line=r.readLine()) != null ) {
			if( !line.startsWith("#") ) {
				Hit h = new Hit(line);
				if( hash.containsKey(h.id) ) {
					throw new IllegalArgumentException(h.id);
				}
				hash.put(h.id, h);
			}
		}
		r.close();
		return hash;
	}
	
	private static String toString(Object[] a) {
    	int iMax = a.length - 1;
        StringBuilder b = new StringBuilder();
        for (int i = 0; i < a.length; i++) {
            b.append(String.valueOf(a[i]));
            if (i < iMax)
                b.append(", ");
        }
        return b.toString();
    }

	static class Hit {
		String geneID, id, line;
		double f1;
		String[] matches;
		
		Hit( String line ) {
			this.line = line;
			String[] split = line.split("\t");
			id = split[1];
			geneID = split[0];
			if( split[split.length-2].equals("NA") ) {
				f1 = Double.NaN;
			} else {
				f1 = Double.parseDouble(split[split.length-2]);
			}
			matches = split[split.length-1].split(";");
			for( int i = 0; i < matches.length; i++ ) {
				matches[i] = matches[i].split(",")[0];
			}		
		}
		
		
		
		public String toString() {
			int idx = line.indexOf('\t');
			idx = line.indexOf('\t',idx+1);
			return line.substring(idx);
		}
	}
}
