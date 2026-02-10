package projects.sigma;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.utils.DoubleList;
import de.jstacs.utils.ToolBox;

public class NonMaximumSuppression {

	public static void main(String[] args) throws IOException {
		
		BufferedReader read = new BufferedReader(new FileReader(args[0]));
		
		int width = Integer.parseInt(args[1]);
		
		
		PrintWriter wr = new PrintWriter(args[0]+"_w"+width+".tsv");
		
		String lastChr = "";
		String lastStrand = "";
		
		double[] region = new double[width];
		
		int off = width/2;
		
		int idx = 0;
		
		String str = null;
		while( (str = read.readLine()) != null ){
			String[] parts = str.split("\t");
			
			String chr = parts[0];
			String strand = parts[2];
			
			int pos = Integer.parseInt(parts[1]);
			double score = Double.parseDouble(parts[3]);
			
			if(!chr.equals(lastChr) || !strand.equals(lastStrand)){
				Arrays.fill(region, Double.NEGATIVE_INFINITY);
			}
			
			region[idx] = score;
			
			int pos2 = pos-off;
			if(strand.equals("-")){
				pos2 = pos+off;
			}

			double max = ToolBox.max(region);
			int idx2 = idx-off;
			if(idx2<0){
				idx2 = region.length+idx2;
			}
			if(region[idx2] != max){
				wr.println(chr+"\t"+(pos2)+"\t"+strand+"\tNA");
			}else{
				wr.println(chr+"\t"+(pos2)+"\t"+strand+"\t"+max);
			}
			
			
			idx++;
			if(idx >= region.length){
				idx = 0;
			}
			
			lastChr = chr;
			lastStrand = strand;
		}
		
		read.close();
		wr.close();
	}

}
