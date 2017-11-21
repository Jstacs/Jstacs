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

package projects;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

/**
 * This class splits fastA files in a user-specified number of fastA-files.
 * It helps to run tools in parallel for instance in Galaxy workflows.
 * 
 * @author Jens Keilwagen
 */
public class FastaSplitter {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		int num = Integer.parseInt(args[1]);
		BufferedReader r = new BufferedReader( new FileReader(args[0]) );
		String delim=null, id = null;
		HashMap<String, Integer> hash = null;
		if( args.length>2 ) {
			delim=args[2];
			hash = new HashMap<String, Integer>();
		}
		BufferedWriter[] w = new BufferedWriter[num];
		int[] stats = new int[num];
		for( int i = 0; i < num; i++ ) {
			w[i] = new BufferedWriter(new FileWriter("split-"+i+".fasta"));			
		}
		
		String line;
		int current = num-1;
		while( (line=r.readLine()) != null ) {
			
			if( line.length()>0 && line.charAt(0) == '>' ) {
				if( delim==null ) {
					current = (current+1) % num;
				} else {
					int h = line.lastIndexOf(delim);
					if( h < 0 ) {
						h = line.length();
					}
					id = line.substring( 1, h );
					Integer idx = hash.get(id);
					if( idx == null ) {
						//this id occurs the first time => find a split with minimal number of sequences so far
						int min = 0;
						for( int i = 1; i < num; i++ ) {
							if( stats[i] < stats[min] ) {
								min = i;
							}
						}
						current = min;
						hash.put(id, current);
					} else {
						//this id was assigned to split idx before
						current = idx;
					}
				}				
				stats[current]++;
			}
			w[current].write(line);
			w[current].newLine();
		}
		r.close();
		for( int i = 0; i < num; i++ ) {
			w[i].close();	
			if( stats[i] == 0 ) {
				new File("split-"+i+".fasta").deleteOnExit();
			}
		}
	}

}
