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
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class DeriveLabels {

	private static class Label{
		private int start;
		private int end;
		private char label;
		public Label(int start, int end){
			this.start = start;
			this.end = end;
		}
	}
	
	private static class Peak{
		private int start;
		private int end;
		private int summit;
		
		public Peak(int start, int end){
			this(start,end,(end-start)/2);
		}
		
		public Peak(int start, int end, int summit){
			this.start = start;
			this.end = end;
			this.summit = summit;
		}
	}
	
	
	private static HashMap<String, ArrayList<Peak>> readPeaks(String file) throws IOException{
		BufferedReader peaks = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
		String str = null;
		HashMap<String, ArrayList<Peak>> list = new HashMap<>();
		while( (str = peaks.readLine()) != null ){
			String[] parts = str.split("\t");
			Peak p = null;
			if(parts.length < 10){
				p = new Peak(Integer.parseInt(parts[1]), Integer.parseInt(parts[2]));
			}else{
				p = new Peak(Integer.parseInt(parts[1]), Integer.parseInt(parts[2]), Integer.parseInt(parts[9]));
			}
			list.putIfAbsent(parts[0], new ArrayList<>());
			list.get(parts[0]).add(p);
		}
		peaks.close();
		for(String key : list.keySet()){
			list.get(key).sort((a,b)->{return Integer.compare(a.start,b.start);});
		}
		return list;
	}
	
	public static void run(String consFile, String relFile, String faiFile, PrintWriter wr, int bin, int region, boolean bySummit) throws FileNotFoundException, IOException {
			
		HashMap<String, ArrayList<Peak>> relaxed = readPeaks(relFile);
		HashMap<String, ArrayList<Peak>> conservative = readPeaks(consFile);
		
		BufferedReader faidx = new BufferedReader(new InputStreamReader(new FileInputStream(faiFile)));
		
		String str = null;
		
		while( (str = faidx.readLine()) != null ){
			String[] parts = str.split("\t");
			String chr = parts[0];
			int len = Integer.parseInt(parts[1]);
			
			ArrayList<Peak> li_relaxed = relaxed.get(chr);
			ArrayList<Peak> li_conservative = conservative.get(chr);
			for(int i=0;i+bin<len;i+=bin){
				int start = i;
				int end = i+bin;
				int endRegion = i+region;
				char label = getLabel(start,end, endRegion,li_relaxed,li_conservative, bySummit);
			
				wr.println(chr+"\t"+start+"\t"+label);
			}
		}
		
		faidx.close();
	}

	private static char getLabel(int start, int end, int endRegion, ArrayList<Peak> li_relaxed, ArrayList<Peak> li_conservative, boolean bySummit) {
		
		char label = 'U';
		if(overlaps(start,endRegion,li_relaxed,1,false)){
			label = 'A';
		}
		if(overlaps(start,endRegion,li_conservative,(end-start)/2,bySummit)){
			label = 'B';
		}
		if(overlaps(start,end,li_conservative,(end-start)/2,true)){
			label = 'S';
		}
		return label;
		
	}

	private static boolean overlaps(int start, int end, ArrayList<Peak> li, int minOverlap, boolean bySummit) {
		if(li == null){
			return false;
		}
		for(int i=0;i<li.size();i++){
			Peak p = li.get(i);
			if(bySummit){
				if(start<=p.start+p.summit && p.start+p.summit <= end){
					return true;
				}
			}else{
				int locStart = Math.max(start, p.start);
				int locEnd = Math.min(end, p.end);
				if(locEnd-locStart>=minOverlap){
					return true;
				}
			}
			if(p.start>end){
				break;
			}
		}
		return false;
	}

}
