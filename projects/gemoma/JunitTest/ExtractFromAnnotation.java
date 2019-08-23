package projects.gemoma.JunitTest;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

/**
 * Tries to extract coverage and introns from annotation for simple (unit) tests.
 * 
 * Assumes that the GFF annotation is (somehow) sorted, i.e. the features of a transcript are listed together.
 * 
 * @author Jens Keilwagen
 *
 */
public class ExtractFromAnnotation {

	static class FeatureComparator implements Comparator<String[]> {

		static FeatureComparator def = new FeatureComparator();
				
		@Override
		public int compare(String[] o1, String[] o2) {
			int d = 0/*o1[2].compareTo(o2[2])*/, i = 3;
			while( d==0 && i<=4 ) {
				int i1 = Integer.parseInt(o1[i]);
				int i2 = Integer.parseInt(o2[i]);
				d = Integer.compare(i1, i2);
				i++;
			}
			return o1[6].charAt(0)=='+'?d:-d;
		}
		
	}
	
	public static void extract( ArrayList<String[]> lines, HashMap<String,HashMap<String,int[]>> combinedIntrons, HashMap<String,HashMap<Integer,int[]>> combinedCov ) {
		if( lines.size() == 0 ) {
			return;
		}
		Collections.sort( lines, FeatureComparator.def );
		
		HashMap<String,int[]> currentIntrons = null;
		HashMap<Integer,int[]> currentCov = null;
		
		String lastType=null, old=null;
		int v=-100;
		for( int i = 0; i < lines.size(); i++ ) {
			String[] split = lines.get(i);
//System.out.println(Arrays.toString(split));
			int st = Integer.parseInt(split[3]);
			int en = Integer.parseInt(split[4]);
			
			String k;
			if( stranded ) {
				k=split[0]+split[6];
			} else {
				k=split[0];
			}
			if( !k.equals(old) ) {
				old=k;
				currentCov = combinedCov.get(old);
				if( currentCov == null ) {
					currentCov = new HashMap<Integer,int[]>();
					combinedCov.put(old, currentCov);
				}
				currentIntrons = combinedIntrons.get(split[0]);
				if( currentIntrons == null ) {
					currentIntrons = new HashMap<String,int[]>();
					combinedIntrons.put(split[0], currentIntrons);
				}
			}
			
			for( int j = st; j <= en; j++ ) {
				int[] a = currentCov.get(j);
				if( a==null ) {
					a=new int[1];
					currentCov.put(j, a);
				}
				a[0]++;
			}
			if( split[2].equals(lastType) ) {
				String key;
				if( split[6].charAt(0)=='+') {
					key = v + "\t" + st + ";";
				} else {
					key = (en+1) + "\t" + v + ";";
				}
				if( stranded ) {
					key += split[6];
				} else {
					key += ".";
				}
				int[] a = currentIntrons.get(key);
				if( a == null ) {
					a = new int[1];
					currentIntrons.put(key, a);
				}
//System.out.println("add intron " + in++ + "\t" + key );
				a[0]++;
				//System.out.println("HIER");
			}
			if( split[6].charAt(0)=='+') {
				v=en+1;
			} else {
				v=st;
			}
			lastType=split[2];
//System.out.println(lastType+"\t" +lastParent);
		}
		
		lines.clear();
//System.out.println();
	}
	
static int in = 0;
	
	static boolean stranded=true;//TODO
	
	/**
	 * @param args
	 * 0 .. annotation
	 * 1 .. type
	 * 2 .. introns out
	 * 3 .. coverage out
	 * 
	 * @throws IOException if we have IO problems
	 */
	public static void main(String[] args) throws IOException {
		HashMap<String,HashMap<String,int[]>> combinedIntrons = new HashMap<String,HashMap<String,int[]>>();
		HashMap<String,HashMap<Integer,int[]>> combinedCov = new HashMap<String,HashMap<Integer,int[]>>();
		
		//read annotation
		BufferedReader r = new BufferedReader( new FileReader( args[0] ));
		String line, type=args[1];
		
		ArrayList<String[]> sameParent = new ArrayList<String[]>();
		String lastParent = null;
		while( (line=r.readLine()) != null ) {
			if( line.equalsIgnoreCase("##FASTA") ) {//http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
				break;  
			}
			if( line.length() == 0 || line.startsWith("#") ) continue;
			String[] split = line.split("\t");
			if( !split[2].endsWith(type) ) {
				continue;
			}
			
			int idx = split[8].indexOf("Parent=");
			if( idx < 0 ) {
				throw new RuntimeException("unknown parent");
			} else {
				if( lastParent == null ) {
					lastParent=split[8].substring(idx, split[8].indexOf(';',idx)+1);
				}
				if( !split[8].contains(lastParent) ) {
System.out.println(lastParent);
					extract( sameParent, combinedIntrons, combinedCov );
					lastParent=split[8].substring(idx, split[8].indexOf(';',idx)+1);
				}
				sameParent.add(split);
			}
		}
		extract( sameParent, combinedIntrons, combinedCov );
		r.close();
		
		BufferedWriter w;
		//write introns
		w = new BufferedWriter( new FileWriter(args[2]) );
		HashMap<Integer,int[]> intronL = new HashMap<Integer, int[]>();
		String[] chr = combinedIntrons.keySet().toArray(new String[0]);
		Arrays.sort(chr);
		for( int i = 0; i < chr.length; i++ ) {
			Iterator<Entry<String,int[]>> it = combinedIntrons.get(chr[i]).entrySet().iterator();//not sorted!?
			while( it.hasNext() ) {
				Entry<String,int[]> e = it.next();
				w.write( chr[i] + "\tAnnnotation\tintron\t" + e.getKey().replaceAll(";", "\t"+e.getValue()[0]+"\t") + "\t.\t.");
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
			}
		}
		w.close();
		
		//write cov
		System.out.println();
		System.out.println(stranded);
		System.out.println(Arrays.toString( combinedCov.keySet().toArray(new String[0])));
		if( stranded ) {
			write( args[3]+"-fwd", combinedCov, '+' );
			write( args[3]+"-rev", combinedCov, '-' );
		} else {
			write( args[3], combinedCov );
		}
		
	}
	
	static void write( String fName, HashMap<String,HashMap<Integer,int[]>> combinedCov, char... f ) throws IOException {
		BufferedWriter w = new BufferedWriter( new FileWriter(fName) );
		String[] chr = combinedCov.keySet().toArray(new String[0]);
		if( f.length > 0 ) {
			ArrayList<String> filtered = new ArrayList<String>();
			for( int i = 0; i < chr.length; i++ ) {
				if( chr[i].charAt(chr[i].length()-1)==f[0] ) {
					filtered.add(chr[i]);
				}
			}
			chr = filtered.toArray(new String[0]);
		}
		Arrays.sort(chr);
		int anz = 0;
		for( int i = 0; i < chr.length; i++ ) {
			String chrom = chr[i];
			if( f.length>0) {
				chrom= chrom.substring(0, chrom.length()-1);
			}
			HashMap<Integer,int[]> currentCov = combinedCov.get(chr[i]);
			Integer[] keys = currentCov.keySet().toArray(new Integer[0]);
			Arrays.sort(keys);
			int start=-1, last=-1, cov=-1;
			for( int j = 0; j<keys.length; j++ ) {
				int[] c = currentCov.get(keys[j]);
				//System.out.println(keys[j] + "\t" + c[0] );
				if( last+1==keys[j] && cov == c[0] ) {
					//proceed
					last++;
				} else {
					if( start> 0 ) {
						w.write( chrom + "\t" + start + "\t" + (last+1) + "\t" + cov );
						w.newLine();
						anz++;
					}
					start=last=keys[j];
					cov=c[0];
				}
			}
			if( start> 0 ) {
				w.write( chrom + "\t" + start + "\t" + (last+1) + "\t" + cov );
				w.newLine();
				anz++;
			}
		}
		w.close();
		System.out.println( "wrote " + Arrays.toString(f) + ":\t " + anz  );
	}
}