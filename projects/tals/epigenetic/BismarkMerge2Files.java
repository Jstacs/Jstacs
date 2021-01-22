package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class BismarkMerge2Files {
	
	public static void main(String[] args) throws IOException {
		String file1=args[0];
		String file2=args[1];
		String outGZFile=args[2];
		
		HashMap<String, HashMap<Integer,String>> tempBismark=new HashMap<>();
		
		BufferedReader BR=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(file1)))));
		String line="";
		String[] splitLine;
		while ((line = BR.readLine()) != null){
			//<chromosome>	<start position>	<end position>	<methylation percentage>	<count methylated>	<count unmethylated>
			//chromosome02    359     359     83.3333333333333        5       1
			//CM/(CM+CU)
				splitLine=line.split("\t");
				
				HashMap <Integer,String> temp=null;
				if(tempBismark.containsKey(splitLine[0])){
					temp=tempBismark.get(splitLine[0]);
				}else {
					temp=new HashMap<>();
				}
				temp.put(Integer.parseInt(splitLine[1]), line);
				
				tempBismark.put(splitLine[0], temp);

		}
		BR.close();

		BufferedReader BR2=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(file2)))));
		line="";
		String[] splitLineF1;
		while ((line = BR2.readLine()) != null){
			//<chromosome>	<start position>	<end position>	<methylation percentage>	<count methylated>	<count unmethylated>
			//chromosome02    359     359     83.3333333333333        5       1
			//CM/(CM+CU)
				splitLine=line.split("\t");
				HashMap <Integer,String> temp=tempBismark.get(splitLine[0]);
				if(!temp.containsKey(Integer.parseInt(splitLine[1]))){
					temp.put(Integer.parseInt(splitLine[1]), line);
				}else{
					splitLineF1=temp.get(Integer.parseInt(splitLine[1])).split("\t");
					double count_methyl=Integer.parseInt(splitLineF1[4])+Integer.parseInt(splitLine[4]);
					double count_unmethyl=Integer.parseInt(splitLineF1[5])+Integer.parseInt(splitLine[5]);
					double methylationLevel=0.0;
					if(count_methyl>0.0){
						methylationLevel=count_methyl/(count_methyl+count_unmethyl)*100;
					}
					
					String newLine=splitLineF1[0]+"\t"+splitLineF1[1]+"\t"+splitLineF1[2]+"\t"+methylationLevel+"\t"+((int)count_methyl)+"\t"+((int)count_unmethyl);
					temp.put(Integer.parseInt(splitLine[1]), newLine);
				}
				
				tempBismark.put(splitLine[0], temp);
		}
		BR2.close();
		
		BufferedWriter BW=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(outGZFile)))));
		for(String chrom : tempBismark.keySet()){
			HashMap <Integer,String> temp=tempBismark.get(chrom);
			List<Integer> sortedByKey = new ArrayList<>(temp.keySet());
			Collections.sort(sortedByKey);
			for (Integer i : sortedByKey) {
				
				BW.write(temp.get(i)+"\n");
			}			
		}
		BW.close();
	}
	
}
