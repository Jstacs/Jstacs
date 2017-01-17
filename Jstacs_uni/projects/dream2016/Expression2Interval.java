package projects.dream2016;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import de.jstacs.utils.Time;
import de.jstacs.utils.ToolBox;


/**
 * Compute the 
 * mean expressions,
 * coefficient of variation, and
 * relative difference of the replicates
 * 
 * @author Jens Keilwagen
 */
public class Expression2Interval {
	
	static class Transcript implements Comparable<Transcript> {

		String id;
		int start;
		double[] tpm;
		double mean, relDiff, cv;
		
		public Transcript(String id, int start) {
			this.id = id;
			this.start = start;
			tpm = new double[2];
		}
		
		void set( int idx, double value ) {
			tpm[idx]=value;
		}
		
		void compute() {
			mean = ToolBox.mean(0, tpm.length, tpm);
			relDiff = mean == 0 ? 0 : Math.abs(tpm[0]-tpm[1])/mean;
			cv = mean == 0 ? 0 : ToolBox.sd(0, tpm.length, tpm) / mean; //coefficient of variation
		}
		
		@Override
		public int compareTo(Transcript o) {
			return Integer.compare(start, o.start);
		}		
	}

	static void add( String fName, int idx, HashMap<String, Transcript> trans ) throws IOException {//TODO
		BufferedReader r = new BufferedReader( new FileReader( fName ) );
		String line;
		r.readLine();
		while( (line=r.readLine()) != null ) {
			String[] split = line.split("\t");
			String[] id = split[1].split(",");
			double v = Double.parseDouble(split[5]);
			for( String name : id ) {
				Transcript t = trans.get(name);
				t.set(idx, v);
			}
		}
		r.close();
	}
	
	/**
	 * Creates genomewide expression interval files.
	 * 
	 * @param args
	 * 0 .. genome annotation
	 * 1 .. expression file 1
	 * 2 .. expression file 2
	 * 3 .. distance
	 * 4 .. reference file
	 * 
	 * @throws IOException
	 */
	public static void main( String[] args ) throws IOException {
		//read anno
		HashMap<String, ArrayList<Transcript>> anno = new HashMap<String, ArrayList<Transcript>>();
		HashMap<String, Transcript> trans = new HashMap<String, Transcript>();
		GZIPInputStream stream = new GZIPInputStream(new FileInputStream(args[0]));
		BufferedReader read = new BufferedReader(new InputStreamReader(stream));
		String line;
		while( (line=read.readLine()) != null ) {
			if( line.charAt(0)!='#'){
				String[] split = line.split("\t");
				if( split[2].equals("transcript") ) {
					int start = split[8].indexOf("ID=")+3;
					int end = split[8].indexOf(';',start);
					String id = split[8].substring(start, end);
					Transcript t = new Transcript(id, Integer.parseInt(split[(split[6].charAt(0)=='+'?3:4)]));
					
					trans.put(id, t);
					ArrayList<Transcript> list = anno.get(split[0]);
					if( list == null ) {
						list = new ArrayList<Transcript>();
						anno.put(split[0], list);
					}
					list.add(t);
				}
			}
		}
		read.close();
		stream.close();

		//read expression
		add(args[1], 0, trans);
		add(args[2], 1, trans);
		
		//sort
		Iterator<String> it = anno.keySet().iterator();
		while( it.hasNext() ) {
			String chr = it.next();
			ArrayList<Transcript> list = anno.get(chr);
			for( int i = 0; i < list.size(); i++ ) {
				list.get(i).compute();
			}
			Collections.sort(list);
		}
		
		String chr=null;
		int s = -100, idx=-100, dist = Integer.parseInt(args[3]);
		stream = new GZIPInputStream(new FileInputStream(args[4]));//ref
		read = new BufferedReader(new InputStreamReader(stream));
		int last = args[1].substring(0, args[1].length()-4).lastIndexOf('.');
		String fName = args[1].substring(0, last);
		System.out.println(fName);
		//System.exit(1);
		FileOutputStream output = new FileOutputStream(fName+"-expression-interval.txt.gz");
		try {
			BufferedWriter w = new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8") );
			Time t = Time.getTimeInstance(null);
			try {
				line = read.readLine();
				do {
					System.out.println(line);
					chr = line.substring(1, line.length()-1);
					ArrayList<Transcript> transChr = anno.get(chr);
					idx=0;
					s=0;
					w.append(line);
					w.newLine();
					
					while( (line=read.readLine()) != null && line.charAt(0)!='[' ) {
						double maxMean = 0, x=0, cv=0;
						if( transChr != null ) {
							int i = idx;
							Transcript tr;
							while( i < transChr.size() && (tr=transChr.get(idx)).start-dist <= s ) {
								if(  tr.start+dist < s ) {
									idx++;
								} else {
									if( tr.start+dist >= s && maxMean < tr.mean ) {
										maxMean = tr.mean;
										x = tr.relDiff;
										cv = tr.cv;
									}
								}
								i++;
							}			
						}
						w.append(maxMean+"\t"+ cv + "\t" + x); //TODO konservierung: min and coefficient of variation
						w.newLine();
						s+=50;
					}
					System.out.println(t.getElapsedTime());
				} while( line != null );
			} finally {
				w.close();
			}
		} finally {
			output.close();
		}
		read.close();
		stream.close();
	}	
}


