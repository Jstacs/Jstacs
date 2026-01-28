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
import java.util.Arrays;
import java.util.HashMap;

import de.jstacs.tools.Protocol;
import de.jstacs.tools.ui.cli.CLI.SysProtocol;
import de.jstacs.utils.IntList;

/**
 * A simple tool for combining intron gff files that might be created in parallel using a compute cluster.
 * Assumes sorted GFF.
 * 
 * @author Jens Keilwagen
 */
public class CombineIntronFiles {

	public static void main(String[] args) throws IOException {
		String[] in  = new String[args.length-1];
		System.arraycopy(args, 1, in, 0, in.length);
		combine(new SysProtocol(), new File(args[0]), in);
	}
	
	public static void combine( Protocol protocol, File out, String... in ) throws IOException {
		protocol.append("files: " + in.length+"\n");
		
		BufferedWriter w = new BufferedWriter( new FileWriter(out) );
		
		//read files
		BufferedReader[] r = new BufferedReader[in.length];
		Intron[] intron = new Intron[in.length];
		String line;
		for( int i = 0; i < in.length; i++ ) {
			protocol.append(i + "\t" + in[i]+"\n");
			r[i] = new BufferedReader( new FileReader( in[i] ) );
			//skip header
			while( (line=r[i].readLine()) != null && line.charAt(0)=='#' ) {
				if( i==0 ) {
					w.append(line);
					w.newLine();
				}
			}
			if( line != null ) intron[i] = new Intron(line);
		}

		HashMap<Integer,int[]> intronL = new HashMap<Integer, int[]>();
		long anz = 0;
		IntList same = new IntList();
		do {
			//determine next intron
			same.clear();
			int a = 0;
			for( int i = 0; i < in.length; i++ ) {
				if( intron[i]!=null ) {
					if( same.length()==0 ) {
						same.add(i);
					} else {
						int comp = intron[same.get(0)].compareTo(intron[i]);
						if( comp >= 0 ) {
							if( comp>0 ) {
								same.clear();
							}
							same.add(i);
						}
					}
					a++;
				}
			}
			
			//combine & read next
			if( same.length()>0 ) {
				Intron rep = intron[same.get(0)];
				String[] split = rep.split;
				int l = rep.end - rep.start;
				int[] stat = intronL.get(l);
				if( stat == null ) {
					stat = new int[1];
					intronL.put(l, stat);
				}
				stat[0]++;
				anz++;
				
				int sum=0;
				for( int j = 0; j < same.length(); j++ ) {
					int i = same.get(j);
					sum += Integer.parseInt(intron[i].split[5]);
					line = r[i].readLine();
					if( line != null ) {
						intron[i]=new Intron( line );
					} else {
						intron[i]=null;
					}
				}
				split[5] = ""+sum;
				for( int j = 0; j < split.length; j++ ) {
					w.write( (j==0?"":"\t") + split[j] );
				}
				w.newLine();
			}
		} while( same.length()>0 );
		w.close();
		
		//print statistic
		protocol.append("\n");
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

	static class Intron implements Comparable<Intron>{
		String[] split;
		int start, end;
		
		Intron( String line ) {
			split = line.split("\t");
			start = Integer.parseInt(split[3]);
			end = Integer.parseInt(split[4]);
		}

		@Override
		public int compareTo(Intron o) {
			int diff = split[0].compareTo(o.split[0]);
			if( diff == 0 ) {
				diff = split[6].compareTo(o.split[6]);
				if( diff == 0 ) {
					diff = Integer.compare(start,o.start);
					if( diff == 0 ) {
						diff = Integer.compare(end,o.end);	
					}
				}
			}
			return diff;
		}
	}
}