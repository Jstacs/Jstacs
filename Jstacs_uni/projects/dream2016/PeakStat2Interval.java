package projects.dream2016;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.IntList;

/**
 * Creates a DNase peak feature file from two peak files.
 *  
 * @author Jens Keilwagen
 */
public class PeakStat2Interval {

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
		HashMap<String,ArrayList<ComparableElement<double[],Integer>>> conservative = getPeaks( args[0] + ".conservative.narrowPeak.gz");
		HashMap<String,ArrayList<ComparableElement<double[],Integer>>> relaxed = getPeaks( args[0] + ".relaxed.narrowPeak.gz");
		
		GZIPInputStream stream = new GZIPInputStream(new FileInputStream(args[1]));
		BufferedReader read = new BufferedReader(new InputStreamReader(stream));
		
		String line, chr=null;
		FileOutputStream output = new FileOutputStream(args[0]+"-DNase-peakStat2interval.txt.gz");
		try {
			BufferedWriter w = new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8") );
			ArrayList<ComparableElement<double[],Integer>> co, re;
			co=re=null;
			ArrayList<ComparableElement<double[],Integer>> currCo = new ArrayList<>(), currRe = new ArrayList<>();
			int b=0, c=0, r=0;
			try {
				while( (line=read.readLine()) != null ) {
					if(line.startsWith("[")){
						System.out.println(line);
						w.append(line);
						w.newLine();
						chr = line.substring(1, line.length()-1);
						co = conservative.get(chr);
						re = relaxed.get(chr);
						currCo.clear();
						currRe.clear();
						b = c = r = 0;
					}else{
						
						double coStat1 = 0.0, coStat2 = 0.0;
						if(co != null){
							while(c<co.size() && co.get(c).getWeight() <= b+50){
								currCo.add(co.get(c));
								c++;
							}

							for(int i=currCo.size()-1;i>=0;i--){
								if(currCo.get(i).getElement()[1] < b){
									currCo.remove(i);
								}
							}

							if(currCo.size() > 0){
								for(int i=0;i<currCo.size();i++){
									double[] curr = currCo.get(i).getElement();
									int start = Math.max(b, (int)curr[0]);
									int end = Math.min(b+50, (int)curr[1]);
									coStat1 += curr[2]*(end-start);
									coStat2 += curr[3]*(end-start);
								}
								coStat1 /= currCo.size()*50;
								coStat2 /= currCo.size()*50;
							}
						}
						
						
						double reStat1 = 0.0, reStat2 = 0.0;
						if(re != null){
							while(r<re.size() && re.get(r).getWeight() <= b+50){
								currRe.add(re.get(r));
								r++;
							}

							for(int i=currRe.size()-1;i>=0;i--){
								if(currRe.get(i).getElement()[1] < b){
									currRe.remove(i);
								}
							}

							if(currRe.size() > 0){
								for(int i=0;i<currRe.size();i++){
									double[] curr = currRe.get(i).getElement();
									int start = Math.max(b, (int)curr[0]);
									int end = Math.min(b+50, (int)curr[1]);
									reStat1 += curr[2]*(end-start);
									reStat2 += curr[3]*(end-start);
									//System.out.println(start+" "+end+" "+Arrays.toString(curr));
								}
								reStat1 /= currRe.size()*50;
								reStat2 /= currRe.size()*50;
							}
						}
						w.append( coStat1+"\t"+reStat1+"\t"+coStat2+"\t"+reStat2 );
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
	public static HashMap<String,ArrayList<ComparableElement<double[],Integer>>> getPeaks( String fileName ) throws FileNotFoundException, IOException {
		HashMap<String,ArrayList<ComparableElement<double[],Integer>>> res = new HashMap<>();
		ArrayList<ComparableElement<double[],Integer>> current;
		GZIPInputStream stream = new GZIPInputStream(new FileInputStream(fileName));
		BufferedReader r = new BufferedReader(new InputStreamReader(stream));
		String line;
		while( (line = r.readLine()) != null ) {
			String[] split = line.split("\t");
			current = res.get(split[0]);
			if( current == null ) {
				current = new ArrayList<>();
				res.put(split[0], current);
			}
			
			int start = Integer.parseInt(split[1]);
			int end = Integer.parseInt(split[2]);
			double stat = Double.parseDouble(split[7]);
			double stat2 = Double.parseDouble(split[6]);
			
			ComparableElement<double[], Integer> el = new ComparableElement<double[], Integer>(new double[]{start,end,stat2,stat}, start);
			
			current.add( el );
		}
		r.close();
		IntList il = new IntList();
		il.sort();
		int i = 0, a = 0;
		Iterator<String> it = res.keySet().iterator();
		while( it.hasNext() ) {
			current = res.get(it.next());
			Collections.sort(current);
			i++;
			a+=current.size();
		}
		System.out.println(fileName+ "\t" + i + "\t" + a);
		return res;
	}
}
