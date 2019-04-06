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
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Pileup {

	
	public static byte[] map;
	public static char[] invMap = new char[]{'A','C','G','T','+','-'};
	
	static{
		map = new byte[256];
		Arrays.fill(map, (byte)-1);
		map['A'] = 0;
		map['C'] = 1;
		map['G'] = 2;
		map['T'] = 3;
	}
	
	
	public static class CovPile{
		String chr;
		int pos;
		int five_prime;
		
		public CovPile(String chr, int pos, int five_prime){
			this.chr = chr;
			this.pos = pos;
			this.five_prime = five_prime;
		}
		
		public String toString(){
			return chr+"\t"+pos+"\t"+five_prime;
		}

		public String getChr() {
			return chr;
		}
		
		public int getFivePrime(){
			return five_prime;
		}
		
		public int getPos(){
			return pos;
		}
	}
	
	public static class Pile extends CovPile{
		
		
		int a,c,g,t,n;
		int in,del;
		
		HashMap<String,Integer> inserts;
		
		public Pile(String chr, int pos, int a, int c, int g, int t, int n, int in, int del, int five_prime, HashMap<String,Integer> inserts) {
			super(chr,pos,five_prime);

			this.a = a;
			this.c = c;
			this.g = g;
			this.t = t;
			this.n = n;
			this.in = in;
			this.del = del;
			this.inserts = inserts;
		}
		
		public String toString(){
			return super.toString()+"\t"+a+"\t"+c+"\t"+g+"\t"+t+"\t"+in+"\t"+del+"\t"+(inserts==null ? "" : inserts.toString());
		}

		public int getCount(int cas) {
			switch(cas){
				case 0: return a;
				case 1: return c;
				case 2: return g;
				case 3: return t;
				case 4: return in;
				case 5: return del;
				case 6: return n;
				default: return -1;
			}
		}

		public int[] getVariants(int ref) {
			double[] temp = new double[]{a,c,g,t,-1,del};
			temp[ref] = -1;
			int[] o = ToolBox.order(temp, true);
			int[] o2 = new int[]{o[0],o[1]};
			return o2;
			/*int max = temp[0];
			int maxIdx = 0;
			for(int i=1;i<temp.length;i++){
				if(temp[i]>max){
					max = temp[i];
					maxIdx = i;
				}
			}
			return maxIdx;*/
		}

		public int getSum() {
			return a+c+g+t+del+n;
		}
		
		public Pair<String,int[]> getMostFrequentInsert(){
			int nIns = 0;
			int n = 0;
			String insStr = null;
			if(inserts != null){

				Iterator<String> it = inserts.keySet().iterator();
				while(it.hasNext()){
					String key = it.next();
					int val = inserts.get(key);
					if(val > nIns){
						nIns = val;
						insStr = key;
					}
					n += val;
				}
			}
			return new Pair<String,int[]>(insStr, new int[]{nIns,n-nIns});
		}
		
		/*public StatEl getStats(int ref){
			int r = getCount(ref);
			int[] o = getVariants(ref);
			int v = getCount(o[0]);
			int n1 = getCount(o[1]);
			int n2 = getSum()-r-v-n1;
			int[] nInsOth = getMostFrequentInsert().getSecondElement();		
			if(r+v+n1+n2<nInsOth[0]){
				System.out.println(r+" "+v+" "+n1+" "+n2+" "+" "+chr+" "+pos+" "+Arrays.toString(nInsOth));
				System.exit(1);
			}
			return new StatEl(r, v, n1, n2,nInsOth[0],nInsOth[1]);
		}*/
		
	}
	
	public static void main(String[] args) throws IOException{
		
		ObjectStream<Pile> ps = new ObjectStream<>(10000);
		
		new Thread( ()->{
			try {
				pileup(args[0], ps, true, false);
				ps.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}).start();
		
		ps.print(System.out);
	}
	
	
	public static void pileup(String bam, ObjectStream<? extends CovPile> piles, boolean variants, boolean coverage) throws IOException {
		
		SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency( ValidationStringency.SILENT );
		
		SamReader sr = srf.open(new File(bam));
		
		SAMRecordIterator samIt = sr.iterator();
		samIt = samIt.assertSorted(SortOrder.coordinate);
		
		byte[] map = new byte[256];
		Arrays.fill(map, (byte)-1);
		map['A'] = 0;
		map['C'] = 1;
		map['G'] = 2;
		map['T'] = 3;
		
		int[][] counts = null;
		HashMap<Integer, LinkedList<String>> insertions = null;
		if(variants){
			counts = new int[10000][4+1+1+1];
			insertions = new HashMap<Integer, LinkedList<String>>();
		}
		int[] totalCount = null;
		if(coverage){
			totalCount = new int[10000];
		}
		int refOff = 0;
		
		int maxEnd = 0;
		
		String oldChr = null;
		
		
		while(samIt.hasNext()){
			SAMRecord rec = samIt.next();
			if(!rec.getReadUnmappedFlag()&&rec.getMappingQuality()>0){
				int refStart = rec.getUnclippedStart();
				String chr = rec.getReferenceName();
				if( (oldChr != null && !chr.equals(oldChr))){
					collect((ObjectStream<CovPile>)piles,counts,totalCount,oldChr,refOff,maxEnd,maxEnd,insertions);
					refOff = 0;
					maxEnd = 0;
					if(variants){
						counts = new int[10000][4+1+1+1];
						insertions.clear();
					}
					if(coverage){
						totalCount = new int[10000];
					}
					
				}
				oldChr = chr;
				
				collect((ObjectStream<CovPile>)piles,counts,totalCount,chr,refOff,refStart,maxEnd,insertions);
				refOff = refStart;
				
				int refEnd = rec.getAlignmentEnd()+1;
				int len = rec.getReadLength();
				
				byte[] bases = rec.getReadBases();
				//System.out.println(rec.getReadString());
				//System.out.println(Arrays.toString(bases));
				int lastRefPos = -1;
				
				if(variants){
					for(int i=1;i<=len;i++){
						int pos = rec.getReferencePositionAtReadPosition(i);
						int insStart = i;
						while(pos==0&&i<len){//insert to reference
							i++;
							pos = rec.getReferencePositionAtReadPosition(i);
						}
						if(i>insStart){
							int idx = (lastRefPos-refOff);
							counts[idx][4]++;

							insertions.putIfAbsent(lastRefPos, new LinkedList<String>());
							insertions.get(lastRefPos).add(rec.getReadString().substring(insStart, i));
						}

						if(lastRefPos >-1 && pos-1 != lastRefPos){//deletion from reference
							for(int j=lastRefPos+1;j<pos;j++){
								int idx = (j-refOff);

								counts[idx][5]++;
							}
						}
						//aligned
						int idx = (pos-refOff);

						if(map[bases[i-1]]>=0){
							counts[idx][map[bases[i-1]]]++;
						}else{
							counts[idx][6]++;
						}

						if(pos > 0){
							lastRefPos = pos;
						}
					}
				}
				
				if(coverage && !rec.getReadNegativeStrandFlag()){//count only 5' end
					int idx = rec.getUnclippedStart()-refOff;
					totalCount[idx]++;
				}else if(coverage && rec.getReadNegativeStrandFlag()){
					int idx = rec.getUnclippedEnd()-refOff;
					totalCount[idx]++;
				}
				
				
				if(refEnd > maxEnd){
					maxEnd = refEnd;
				}
			}
		}
		collect((ObjectStream<CovPile>)piles,counts,totalCount,oldChr,refOff,maxEnd,maxEnd,insertions);
		sr.close();

	}

	private static void collect(ObjectStream<CovPile> out, int[][] counts, int[] totalCounts, String chr, int refOff, int refStart, int maxEnd,HashMap<Integer, LinkedList<String>> inserts) {
		for(int i=refOff;i<Math.min(refStart, maxEnd);i++){
			int locIdx = i-refOff;
			
			HashMap<String,Integer> insMap = null;
			
			if(inserts != null && inserts.containsKey(i)){
				insMap = new HashMap<>();
				LinkedList<String> myIns = inserts.remove(i);
				Iterator<String> insIt = myIns.iterator();
				while(insIt.hasNext()){
					String key = insIt.next();
					if(insMap.containsKey(key)){
						insMap.put(key, insMap.get(key)+1);
					}else{
						insMap.put(key, 1);
					}
				}
				
			}
			if(counts != null){
				Pile pile = new Pile(chr, i, counts[locIdx][0], counts[locIdx][1], counts[locIdx][2], counts[locIdx][3], counts[locIdx][6], counts[locIdx][4], counts[locIdx][5], totalCounts != null ? totalCounts[locIdx] : -1, insMap);
				out.add(pile);
			}else if(totalCounts != null){
				CovPile pile = new CovPile(chr, i, totalCounts[locIdx]);
				out.add(pile);
			}
		}
		
		if(counts != null){
			if(refStart-refOff>0 && refStart-refOff < counts.length){
				System.arraycopy(counts, refStart-refOff, counts, 0, counts.length-(refStart-refOff));
			}
			for(int i=Math.max(0, counts.length-(refStart-refOff));i<counts.length;i++){
				counts[i] = new int[counts[i].length];
			}
		}
		if(totalCounts != null){
			if(refStart-refOff>0 && refStart-refOff < totalCounts.length){
				System.arraycopy(totalCounts, refStart-refOff, totalCounts, 0, totalCounts.length-(refStart-refOff));
			}
			for(int i=Math.max(0, totalCounts.length-(refStart-refOff));i<totalCounts.length;i++){
				totalCounts[i] = 0;
			}
		}


	}

}
