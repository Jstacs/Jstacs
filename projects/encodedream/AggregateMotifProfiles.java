package projects.encodedream;
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
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.HashMap;

import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;

public class AggregateMotifProfiles {

	public static void run(BufferedReader read, PrintStream out, int bin, String faiFile) throws IOException {
		
		//BufferedReader read = new BufferedReader(new InputStreamReader(System.in));
		
		BufferedReader faidx = new BufferedReader(new InputStreamReader(new FileInputStream(faiFile)));
		String str = null;
		
		HashMap<String,Integer> sizes = new HashMap<>();
		
		while( (str = faidx.readLine()) != null ){
			String[] parts = str.split("\t");
			String chr = parts[0];
			int len = Integer.parseInt(parts[1]);
			sizes.put(chr,len);
		}
		
		String lastChr = "";
		str = null;
		DoubleList scores = new DoubleList();
		int start = 0;
		int last = 0;
		while( (str = read.readLine()) != null ){
			String[] parts = str.split("\t");
			//System.out.println(parts.length);
			if(parts.length==4){
				String chr = parts[0];
				
				if(!chr.equals(lastChr)){
					last = start - start%bin;
					if(lastChr.length()>0 && last+bin<=sizes.get(lastChr)){
						print(lastChr,last,scores, out);
					}
					if(lastChr.length()>0 && last+bin<=sizes.get(lastChr)){
						print(lastChr,last+bin,sizes.get(lastChr)-bin,bin,out);
					}
					scores.clear();
				}
				start = Integer.parseInt(parts[1]);
				double fwd = Double.parseDouble(parts[2]);
				double rev = Double.parseDouble(parts[3]);
				
				if(scores.length()==0){
					if(start>bin){
						print(chr,0,start-bin,bin,out);
						last = start;
					}
				}
				
				if(start%bin==0){
					if(last+bin != start-bin){
						print(chr,last+bin,start-2*bin,bin,out);
					}
					last = start -  bin;
					print(chr,last,scores, out);
					scores.clear();
				}
				scores.add(fwd);
				scores.add(rev);
				lastChr = chr;
			}
		}
		
		last = start - (start%bin==0 ? bin : start%bin);
		if(last+bin<=sizes.get(lastChr)){
			print(lastChr,last,scores, out);
		}
		if(last+bin<=sizes.get(lastChr)){
			print(lastChr,last+bin,sizes.get(lastChr)-bin,bin,out);
		}
		//read.close();
	}

	private static void print(String chr, int start, int end, int bin, PrintStream out){
		for(int i=start;i<=end;i+=bin){
			out.println(chr+"\t"+i+"\tNA\tNA");
		}
	}
	
	private static void print(String chr, int off, DoubleList scores, PrintStream out) {
		if(scores.length()==0){
			return;
		}
		double[] temp = scores.toArray();
		double max = ToolBox.max(temp);
		double ls = Normalisation.getLogSum(temp);
		
		out.println(chr+"\t"+off+"\t"+max+"\t"+ls);
		
	}

}
