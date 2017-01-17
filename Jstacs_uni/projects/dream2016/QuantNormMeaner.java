package projects.dream2016;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class QuantNormMeaner {

	public static void main(String[] args) throws Exception {
		
		BufferedReader orig = new BufferedReader( new InputStreamReader(new GZIPInputStream(new FileInputStream(args[0]))) );
		
		BufferedReader normed = new BufferedReader( new InputStreamReader(new GZIPInputStream(new FileInputStream(args[1]))) );
		BufferedReader normed2 = new BufferedReader( new InputStreamReader(new GZIPInputStream(new FileInputStream(args[1]))) );
		
		BufferedWriter wr = new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args[1]+"-mean.gz"))) );
		
		
		
		int col = Integer.parseInt(args[2]);
		
		String normedStr = null;
		int i=0;
		
		String origStr = orig.readLine();
		String[] split = origStr.split("\t");
		double origVal = Double.parseDouble(split[col]);
		
		while( true ){
			
			double nextOrigVal = origVal;
			
			double normedMean = 0;
			int j=i;

			while( origVal==nextOrigVal ){
				
				origStr = orig.readLine();
				if(origStr != null){
					split = origStr.split("\t");
					nextOrigVal = Double.parseDouble(split[col]);
				}
				normedStr = normed.readLine();
				String[] normedSplit = normedStr.split("\t");
				double val = Double.parseDouble(normedSplit[col]);
				normedMean += val;
				j++;
				
				if(origStr == null){
					break;
				}
			}
			normedMean /= (j-i);

			
			origVal = nextOrigVal;
			
			
			for(;i<j;i++){
				String curr = normed2.readLine();
				String[] currSplit = curr.split("\t");
				currSplit[col] = ""+normedMean;

				wr.append(String.join("\t", currSplit) );
				wr.newLine();
			}
			
			if(origStr == null){
				break;
			}
			
		}
		orig.close();
		normed.close();
		normed2.close();
		wr.close();

	}

}
