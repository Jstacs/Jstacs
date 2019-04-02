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

/**
 * A simple tool for combining intron gff files that might be created in parallel using a compute cluster.
 * 
 * @author Jens Keilwagen
 */
public class CombineIntronFiles {

	public static void main(String[] args) throws IOException {
		BufferedWriter w = new BufferedWriter( new FileWriter(args[0]) );
		
		//read files
		BufferedReader r;
		HashMap<String,HashMap<String,int[]>> combined = new HashMap<String,HashMap<String,int[]>>();
		HashMap<String,int[]> current;
		
		String line;
		for( int i = 1; i < args.length; i++ ) {
			System.out.println((i-1) + "\t" + args[i]);
			r = new BufferedReader( new FileReader( args[i] ) );
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
			System.out.println(il[j] + "\t" + stat[0] + "\t" + (all/anz) );
		}
		
	}
}
