package projects.dream2016;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import de.jstacs.utils.IntList;

/**
 * Creates a DNase peak feature file from two peak files.
 *  
 * @author Jens Keilwagen
 */
public class Peak2Interval {

	/**
	 * 
	 * @param args
	 * 0 ... prefix of feature files
	 * 1 ... reference file
	 * 
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {
		HashMap<String,IntList> conservative = getPeaks( args[0] + ".conservative.narrowPeak.gz");
		HashMap<String,IntList> relaxed = getPeaks( args[0] + ".relaxed.narrowPeak.gz");
		
		GZIPInputStream stream = new GZIPInputStream(new FileInputStream(args[1]));
		BufferedReader read = new BufferedReader(new InputStreamReader(stream));
		
		String line, chr=null;
		FileOutputStream output = new FileOutputStream(args[0]+"-DNase-peak2interval.txt.gz");
		try {
			BufferedWriter w = new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8") );
			IntList co, re;
			co=re=null;
			int b, c, r, d1, d2;
			b=c=r=-999999999;
			try {
				while( (line=read.readLine()) != null ) {
					if(line.startsWith("[")){
						System.out.println(line);
						w.append(line);
						w.newLine();
						chr = line.substring(1, line.length()-1);
						b=c=r=0;
						co = conservative.get(chr);
						re = relaxed.get(chr);
					}else{
						//conservative distance (in bp)
						while( co!=null && c < co.length() && b > co.get(c) )  {
							c++;
						}
						d1 = d2 = Integer.MAX_VALUE; 
						if( c > 0 ) {
							d1=b-co.get(c-1);
						}
						if( co != null && c < co.length() ) {
							d2=Math.max(0, co.get(c)-(b+50));
						}
						int D = Math.min(d1, d2);
						
						//relaxed distance (in bp)
						while( re!=null && r < re.length() && b > re.get(r) )  {
							r++;
						}
						d1 = d2 = Integer.MAX_VALUE; 
						if( r > 0 ) {
							d1=b-re.get(r-1);
						}
						if( re != null && r < re.length() ) {
							d2=Math.max(0,re.get(r)-(b+50));
						}
						
						w.append( D + "\t" + Math.min(d1, d2) );
						w.newLine();
						b+=50;
					}
				}
			} finally {
				w.close();
			}
		} finally {
			output.close();
		}
		read.close();
		stream.close();
	}
	
	// read peak file
	public static HashMap<String,IntList> getPeaks( String fileName ) throws FileNotFoundException, IOException {
		HashMap<String,IntList> res = new HashMap<String,IntList>();
		IntList current;
		GZIPInputStream stream = new GZIPInputStream(new FileInputStream(fileName));
		BufferedReader r = new BufferedReader(new InputStreamReader(stream));
		String line;
		while( (line = r.readLine()) != null ) {
			String[] split = line.split("\t");
			current = res.get(split[0]);
			if( current == null ) {
				current = new IntList();
				res.put(split[0], current);
			}
			current.add( Integer.parseInt(split[1]) + Integer.parseInt(split[9]));
		}
		r.close();
		IntList il = new IntList();
		il.sort();
		int i = 0, a = 0;
		Iterator<String> it = res.keySet().iterator();
		while( it.hasNext() ) {
			current = res.get(it.next());
			current.sort();
			i++;
			a+=current.length();
		}
		System.out.println(fileName+ "\t" + i + "\t" + a);
		return res;
	}
}
