package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

/**
 * A simple tool for combining coverage files that might be created in parallel using a compute cluster.
 * 
 * @author Jens Keilwagen
 */
public class CombineCoverageFiles {

	public static void main(String[] args) throws IOException {
		System.out.println("files: " + (args.length-1));
		
		BufferedWriter w = new BufferedWriter( new FileWriter(args[0]) );
		//read files
		BufferedReader r;
		HashMap<String,HashMap<Integer,int[]>> combined = new HashMap<String,HashMap<Integer,int[]>>();
		HashMap<Integer,int[]> current;
		
		String line;
		for( int i = 1; i < args.length; i++ ) {
			System.out.println((i-1) + "\t" + args[i]);
			r = new BufferedReader( new FileReader( args[i] ) );
			//skip header
			line=r.readLine();
			if( i==1 ) {
				w.append(line);
				w.newLine();
			}
			
			//add to hash
			String old=null;
			current = null;
			while( (line=r.readLine()) != null ) {
				String[] split = line.split("\t");
				if( !split[0].equals(old) ) {
					old=split[0];
					current = combined.get(split[0]);
					if( current == null ) {
						current = new HashMap<Integer,int[]>();
						combined.put(split[0], current);
					}
				}
				int st = Integer.parseInt(split[1]);
				int en = Integer.parseInt(split[2]);
				int cov = Integer.parseInt(split[3]);
				for( int j = st; j < en; j++ ) {
					int[] a = current.get(j);
					if( a==null ) {
						a=new int[1];
						current.put(j, a);
					}
					a[0] += cov;
				}
			}
			r.close();
		}

		//write
		System.out.println();
		System.out.println("write");
		String[] chr = combined.keySet().toArray(new String[0]);
		Arrays.sort(chr);
		for( int i = 0; i < chr.length; i++ ) {
			current = combined.get(chr[i]);
			Integer[] keys = current.keySet().toArray(new Integer[0]);
			Arrays.sort(keys);
			int start=-1, last=-1, cov=-1;
			for( int j = 0; j<keys.length; j++ ) {
				int[] c = current.get(keys[j]);
				//System.out.println(keys[j] + "\t" + c[0] );
				if( last+1==keys[j] && cov == c[0] ) {
					//proceed
					last++;
				} else {
					if( start> 0 ) {
						w.write( chr[i] + "\t" + start + "\t" + (last+1) + "\t" + cov );
						w.newLine();
					}
					start=last=keys[j];
					cov=c[0];
				}
			}
			if( start> 0 ) {
				w.write( chr[i] + "\t" + start + "\t" + (last+1) + "\t" + cov );
				w.newLine();
			}
		}
		w.close();
	}
}
