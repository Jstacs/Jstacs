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
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import de.jstacs.tools.Protocol;
import de.jstacs.tools.ui.cli.CLI.SysProtocol;

/**
 * A simple tool for combining intron gff files that might be created in parallel using a compute cluster.
 * 
 * @author Jens Keilwagen
 */
public class CombineIntronFiles {

	public static void main(String[] args) throws IOException {
		String[] in  = new String[args.length-1];
		System.arraycopy(args, 1, in, 0, in.length);
		combine(new SysProtocol(), args[0], in);
	}
	
	public static void combine( Protocol protocol, String out, String... in ) throws IOException {
		BufferedWriter w = new BufferedWriter( new FileWriter(out) );
		
		//read files
		BufferedReader r;
		HashMap<String,HashMap<String,int[]>> combined = new HashMap<String,HashMap<String,int[]>>();
		HashMap<String,int[]> current;
		
		String line;
		for( int i = 0; i < in.length; i++ ) {
			protocol.append(i + "\t" + in[i]+"\n");
			r = new BufferedReader( new FileReader( in[i] ) );
			//skip header
			while( (line=r.readLine()) != null && line.charAt(0)=='#' ) {
				if( i==1 ) {
					w.append(line);
					w.newLine();
				}
			}
			
			//add to hash
			while( (line=r.readLine()) != null ) {
				String[] split = line.split("\t");
				current = combined.get(split[0]);
				if( current == null ) {
					current = new HashMap<String,int[]>();
					combined.put(split[0], current);
				}
				int anz = Integer.parseInt(split[5]);
				String key = split[3] + "\t" + split[4] + ";" + split[6];
				int[] a = current.get(key);
				if( a == null ) {
					a = new int[1];
					current.put(key, a);
				}
				a[0] += anz;
			}
			r.close();
		}

		//write
		HashMap<Integer,int[]> intronL = new HashMap<Integer, int[]>();
		long anz = 0;
		String[] chr = combined.keySet().toArray(new String[0]);
		Arrays.sort(chr);
		for( int i = 0; i < chr.length; i++ ) {
			Iterator<Entry<String,int[]>> it = combined.get(chr[i]).entrySet().iterator();//not sorted!?
			while( it.hasNext() ) {
				Entry<String,int[]> e = it.next();
				w.write( chr[i] + "\tRNAseq\tintron\t" + e.getKey().replaceAll(";", "\t"+e.getValue()[0]+"\t") + "\t.\t.");
				w.newLine();
				
				//for statistic
				String h = e.getKey();
				h=h.substring(0, h.indexOf(';'));
				int idx = h.indexOf('\t');
				int l = Integer.parseInt(h.substring(idx+1))-Integer.parseInt(h.substring(0,idx));
				int[] stat = intronL.get(l);
				if( stat == null ) {
					stat = new int[1];
					intronL.put(l, stat);
				}
				stat[0]++;
				anz++;
			}
		}
		w.close();
		
		//print statistic
		Integer[] il = new Integer[intronL.size()];
		intronL.keySet().toArray(il);
		Arrays.sort(il);
		double all = 0;
		for( int j = 0; j < il.length; j++ ) {
			int[] stat = intronL.get(il[j]);
			all += stat[0];
			protocol.append(il[j] + "\t" + stat[0] + "\t" + (all/anz) +"\n");
		}
		
	}
}
