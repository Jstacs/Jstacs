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
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

import org.broad.igv.bbfile.BBFileHeader;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import de.jstacs.utils.DoubleList;
import de.jstacs.utils.ToolBox;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import projects.encodedream.Pileup.CovPile;

public class Coverage {
	
	
	public static void main(String[] args) {
		
		
		ObjectStream<CovPile> ps2 = new ObjectStream<>(10000);
		
		new Thread( ()->{
			try {
				Pileup.pileup(args[0], ps2,false,true);
				ps2.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}).start();

		double lambdaBG = estimateLambdaBG(ps2,args[0]);
		//System.out.println(lambdaBG);
		
		ObjectStream<CovPile> ps = new ObjectStream<CovPile>(10000);
		new Thread( ()->{
			try {
				Pileup.pileup(args[0], ps, false, true);
				ps.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}).start();
		
		coverage(ps, lambdaBG, args[0], System.out, 50);
		
	}
	
	public static double estimateLambdaBG(ObjectStream<CovPile> ps, String bam){
		
		SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency( ValidationStringency.SILENT );
		SamReader sr = srf.open(new File(bam));
		SAMSequenceDictionary dict = sr.getFileHeader().getSequenceDictionary();
		
		
		double sum = 0;
		double n = 0;
		
		String lastChr = "";
		while(ps.hasNext()){
			CovPile pile = ps.next();
			String chr = pile.getChr();
			double count = pile.getFivePrime();
			
			if(!chr.equals(lastChr)){
				n+= dict.getSequence(chr).getSequenceLength();
			}
			
			sum += count;
			
			lastChr = chr;
		}
		
		return sum/n;
	}
	
	public static void coverage(String faiFile, String bigwig, PrintStream out, int bin) throws NumberFormatException, IOException{
		BigWigAccessor bwa = new BigWigAccessor(bigwig);
		
		BufferedReader read = new BufferedReader(new FileReader(faiFile));
		String str = null;
		while( (str = read.readLine()) != null ){
			String[] parts = str.split("\t"); 
			String chr = parts[0];
			int len = Integer.parseInt(parts[1]);
			double[] orig = bwa.getProfileInRegion(chr, 0, len);
			
			double[][] temp = aggregate(orig, bin);
			print(chr,temp, out, bin);
			
		}
		
	}
	
	public static void coverageMeth(String faiFile, String bedFileGz, PrintStream out, int bin) throws NumberFormatException, FileNotFoundException, IOException {
        SortedBedAccessor sba = new SortedBedAccessor(bedFileGz);
        BufferedReader read = new BufferedReader(new FileReader(faiFile));
        String str = null;
        while ((str = read.readLine()) != null) {
            String[] parts = str.split("\t");
            String chr = parts[0];
            int len = Integer.parseInt(parts[1]);
            double[] orig = sba.getProfileForChromosome(chr, len);
            double[][] temp = aggregateMeth(orig, bin);
            print(chr, temp, out, bin);
        }
        read.close();
    }
	
	public static void coverage(ObjectStream<CovPile> ps, double lambdaBG, String bam, PrintStream out, int bin){
		SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency( ValidationStringency.SILENT );
		SamReader sr = srf.open(new File(bam));
		SAMSequenceDictionary dict = sr.getFileHeader().getSequenceDictionary();
		
		//List<SAMSequenceRecord> sequences = dict.getSequences();
		
		double[] orig = null;
		
		String lastChr = "";
		while(ps.hasNext()){
			CovPile pile = ps.next();
			String chr = pile.getChr();
			int pos = pile.getPos();
			double count = pile.getFivePrime();
			
			if(!chr.equals(lastChr)){
				
				if(orig != null){
					double[][] temp = process(orig,lambdaBG, bin);
					print(lastChr,temp, out, bin);
				}
				
				if(lastChr.length() > 0 && dict.getSequenceIndex(chr) != dict.getSequenceIndex(lastChr)+1){
					for(int k=dict.getSequenceIndex(lastChr)+1;k<dict.getSequenceIndex(chr);k++){
						orig = new double[dict.getSequence(k).getSequenceLength()];
						double[][] temp = process(orig,lambdaBG,bin);
						print(dict.getSequence(k).getSequenceName(),temp,out,bin);
					}
				}
				
				orig = new double[dict.getSequence(chr).getSequenceLength()];
				//System.out.println(chr);
			}
			orig[pos-1] = count;
			
			lastChr = chr;
		}
		double[][] temp = process(orig,lambdaBG, bin);
		print(lastChr, temp, out, bin);
		for(int i = dict.getSequenceIndex(lastChr) + 1; i < dict.size(); ++i) {
            orig = new double[dict.getSequence(i).getSequenceLength()];
            temp = process(orig, lambdaBG, bin);
            print(dict.getSequence(i).getSequenceName(), temp, out, bin);
        }
	}

	private static void print(String chr, double[][] temp, PrintStream out, int bin) {
		for(int i=0;i<temp.length;i++){
			out.print(chr);
			out.print("\t");
			out.print(i*bin);
			for(int j=0;j<temp[i].length;j++){
				out.print("\t");
				out.print(temp[i][j]);
			}
			out.println();
		}
	}

	private static double[][] process(double[] orig, double lambdaBG, int bin) {
		//"poisson" normalization cf. MACS
		double[] processed = new double[orig.length];
		{
			double[] sum = new double[]{0.0,0.0};
			double[] n = new double[]{0.0,0.0};
			int[] l = new int[]{Math.min(5000,orig.length),Math.min(10000,orig.length)};

			for(int j=0;j<sum.length;j++){
				for(int i=0;i<l[j];i++){
					sum[j] += orig[i];
					n[j]++;
				}
			}

			for(int i=0;i<orig.length;i++){
				for(int j=0;j<l.length;j++){
					if(i+l[j]<orig.length){
						sum[j] += orig[i+l[j]];
					}else{
						n[j]--;
					}
					if(i-l[j]-1>=0){
						sum[j] -= orig[i-l[j]-1];
					}else{
						n[j]++;
					}
				}
				double lambda = lambdaBG;
				for(int j=0;j<sum.length;j++){
					double temp = sum[j]/n[j];
					if(temp>lambda){
						lambda = temp;
					}
				}
				processed[i] = orig[i]/lambda;
			}
		}
		
		//smooth
		double sum = 0;
		double n = 0;
		int l = 75;
		for(int i=0;i<l;i++){
			sum += processed[i];
			n++;
		}
		for(int i=0;i<processed.length;i++){
			if(i+l<processed.length){
				sum += processed[i+l];
			}else{
				n--;
			}
			if(i-l-1>=0){
				sum -= processed[i-l-1];
			}else{
				n++;
			}
			orig[i] = sum/n;
		}
		
		return aggregate(orig, bin);
	}
	
	private static double[][] aggregateMeth(double[] orig, int bin) {
		
        int broad = 1000;
        
        DoubleList temp = new DoubleList(bin + 1);
        double[][] res = new double[orig.length / bin][];
        for(int i = 0; i + bin <= orig.length; i += bin) {
            temp.clear();
            for (int j = i; j < i + bin; ++j) {
                if (!Double.isNaN(orig[j])) {
                    temp.add(orig[j]);
                }
            }
            double min = 0.0;
            double average = 0.0;
            double max = 0.0;
            if (temp.length() > 0) {
                min = temp.min(0, temp.length());
                average = temp.mean(0, temp.length());
                max = temp.max(0, temp.length());
            }
            (res[i / bin] = new double[8])[0] = min;
            res[i / bin][1] = average;
            res[i / bin][2] = max;
        }
        for(int i = 0; i < orig.length; ++i) {
            if (Double.isNaN(orig[i])) {
                orig[i] = 0.0;
            }
        }
        for(int i = 0; i + bin <= orig.length; i += bin) {
            double mean = ToolBox.mean(i, i + bin, orig);
            int broadStart = Math.max(0, i - broad);
            int broadEnd = Math.min(orig.length, i + broad);
            double bMeanAfter = ToolBox.mean(i, broadEnd, orig);
            double bMeanBefore = ToolBox.mean(broadStart, i + 1, orig);
            double bMaxAfter = ToolBox.max(i, broadEnd, orig);
            double bMaxBefore = ToolBox.max(broadStart, i + 1, orig);
            res[i / bin][3] = mean;
            res[i / bin][4] = bMeanAfter;
            res[i / bin][5] = bMeanBefore;
            res[i / bin][6] = bMaxAfter;
            res[i / bin][7] = bMaxBefore;
        }
        return res;
    }
	
	
	
	private static double[][] aggregate(double[] orig, int bin) {
		//aggregate

		int broad = 1000;

		double[][] res = new double[orig.length/bin][];
		for(int i=0;i+bin<=orig.length;i+=bin){
			double min = ToolBox.min(i, i+bin, orig);
			double median = ToolBox.median(i,i+bin,orig);

			int broadStart = Math.max(0, i-broad);
			int broadEnd = Math.min(orig.length, i+broad);
			double bMinAfter = ToolBox.min(i,broadEnd,orig);
			double bMinBefore = ToolBox.min(broadStart,i+1,orig);
			double bMaxAfter = ToolBox.max(i,broadEnd,orig);
			double bMaxBefore = ToolBox.max(broadStart,i+1,orig);

			double orange = orange(i,i+bin,orig);
			double stepsUp = mostMonotonSteps(i, i+bin, orig, 1.0);
			double stepsDown = mostMonotonSteps(i, i+bin, orig, -1.0);
			res[i/bin] = new double[]{min,median,bMinAfter,bMinBefore,bMaxAfter,bMaxBefore,orange,stepsUp,stepsDown};
		}
		return res;
	}

	private static int orange( int start, int end, double[] profile ) {
		int different = 0;
		start = Math.max(0, start);
		end = Math.min(profile.length, end);
		for( int s = start+1; s < end; s++ ) {
			if( profile[s-1] != profile[s] ) {
				different++;
			}
		}
		return different;
	}
	
	private static int mostMonotonSteps( int start, int end, double[] profile, double vz ) {
		int num=0, max=0;
		start = Math.max(0, start);
		end = Math.min(profile.length, end);
		for( int s = start+1; s < end; s++ ) {
			if( vz*profile[s-1] < vz*profile[s] ) {
				num++;
			} else if( vz*profile[s-1] > vz*profile[s] ) {
				if( num > max ) {
					max=num;
				}
				num=0;
			}
		}
		return max;
	}	

	
	private static class BigWigAccessor {

		private BBFileReader reader;
		
		public BigWigAccessor(String bigWigFile) throws IOException {
			reader = new BBFileReader(bigWigFile);
			BBFileHeader header = reader.getBBFileHeader();
			
			if(!header.isHeaderOK()){
				throw new RuntimeException("Header not OK");
			}
			
			if(!header.isBigWig()){
				throw new RuntimeException("No Bigwig");
			}
		}
		
		public double[] getProfileInRegion(String chr, int start, int end){
			double[] res = new double[end-start];
			fillProfileInRegion( chr, start, end, res );
			return res;
		}
		
		public void fillProfileInRegion(String chr, int start, int end, double[] res) {
			BigWigIterator it = reader.getBigWigIterator(chr,start,chr,end,false);
			
			while(it.hasNext()){
				WigItem item = it.next();
				Arrays.fill( res, Math.max(0, item.getStartBase()-start), Math.min(item.getEndBase(),end)-start, item.getWigValue() );
				/*
				int s = item.getStartBase();
				int e = item.getEndBase();
				double v = item.getWigValue();
				for(int i=Math.max(start, s);i<Math.min(e, end);i++){
					//System.out.println(i+"\t"+(i-start)+"\t"+v);
					res[i-start] = v;
				}*/
			}		
		}	
	}
	
	
	private static class SortedBedAccessor{
	    private String bedFileGz;
	    
	    public SortedBedAccessor(final String bedFileGz) {
	        this.bedFileGz = bedFileGz;
	    }
	    
	    public double[] getProfileForChromosome(String chr, int len) throws FileNotFoundException, IOException {
	        final double[] res = new double[len];
	        Arrays.fill(res, Double.NaN);
	        BufferedReader read = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(this.bedFileGz))));
	        String chr2 = chr + "\t";
	        String str = null;
	        while ((str = read.readLine()) != null) {
	            if (str.startsWith(chr2)) {
	                String[] parts = str.split("\t");
	                if (!chr.equals(parts[0])) {
	                	read.close();
	                    throw new RuntimeException();
	                }
	                int pos = Integer.parseInt(parts[1]);
	                double perc = Double.parseDouble(parts[10]) / 100.0;
	                double n = Double.parseDouble(parts[9]);
	                double methN = n * perc;
	                res[pos] = methN / (n + 1.0) * 100.0;
	            }
	        }
	        read.close();
	        return res;
	    }
	}

	
}
