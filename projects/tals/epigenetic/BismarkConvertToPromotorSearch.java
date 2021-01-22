package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;


public class BismarkConvertToPromotorSearch {
	public static void main(String[] args) throws Exception {
		String bismarkFile=args[0];
		String promotorFasta=args[1];
		String outGZFile=args[2];
		
		HashMap<String, HashMap<Integer,String>> tempBismark=new HashMap<>();
		
		BufferedReader BR=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(bismarkFile)))));
		String line="";
		String[] splitLine;
		
		while ((line = BR.readLine()) != null){
			//<chromosome>	<start position>	<end position>	<methylation percentage>	<count methylated>	<count unmethylated>
			//chromosome02    359     359     83.3333333333333        5       1
			//CM/(CM+CU)
			if(line.startsWith("Chr")){
				splitLine=line.split("\t");
				
				HashMap <Integer,String> temp=null;
				if(tempBismark.containsKey(splitLine[0])){
					temp=tempBismark.get(splitLine[0]);
				}else {
					//System.out.println(splitLine[0]);
					
					temp=new HashMap<>();
				}
				temp.put(Integer.parseInt(splitLine[1]), line);
				
				tempBismark.put(splitLine[0], temp);
			}
		}
		BR.close();
		
		BufferedReader FA=new BufferedReader(new FileReader(promotorFasta));
		BufferedWriter BW=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(outGZFile)))));
		line="";
		
		String gene="";
		String chrom="";
		int startPos=-1;
		int endPos=-1;
		String[] splitHeader;
		String[] splitArea;
		String[] splitPos;
		String[] splitBismark;

		int startBismark=-1;
		boolean strand;
		while ((line = FA.readLine()) != null){
			if(line.startsWith(">")){
				splitHeader=line.split(" ");
				splitArea=splitHeader[2].split(":");

				if(splitArea[2].equals("+")){
					strand=true;
				}else{
					strand=false;
				}

				splitPos=splitArea[1].split("-");
				gene=splitHeader[1];
				
				chrom=splitArea[0];

				startPos=Integer.parseInt(splitPos[0]);
				
				endPos=Integer.parseInt(splitPos[1]);
				
				int k=0;
				int pos=-1;
				for(int i=startPos;i<endPos;i++){
					if(strand==true){
						if(tempBismark.containsKey(chrom)){
							if(tempBismark.get(chrom).containsKey(i+1)){
								splitBismark=tempBismark.get(chrom).get(i+1).split("\t");
								startBismark=Integer.parseInt(splitBismark[1])-startPos;
								BW.write(gene+"\t"+startBismark+"\t"+startBismark+"\t"+splitBismark[3]+"\t"+splitBismark[4]+"\t"+splitBismark[5]+"\n");
							}
							
						}
					}
					if(strand==false){
						pos=endPos-k;
						
						if(tempBismark.containsKey(chrom)){
							if(tempBismark.get(chrom).containsKey(pos)){
								splitBismark=tempBismark.get(chrom).get(pos).split("\t");
								startBismark=k+1;
								BW.write(gene+"\t"+startBismark+"\t"+startBismark+"\t"+splitBismark[3]+"\t"+splitBismark[4]+"\t"+splitBismark[5]+"\n");
							}
						}
					}
					k++;
				}
			}
		}
		FA.close();
		BW.close();
	}
}
