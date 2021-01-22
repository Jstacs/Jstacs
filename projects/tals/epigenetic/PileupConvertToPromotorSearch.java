package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;


public class PileupConvertToPromotorSearch {
	public static void main(String[] args) throws Exception {
		//input is pileup-File
		String pileupFile=args[0];
		String promotorFasta=args[1];
		String outFile=args[2];
		
		HashMap<String, HashMap<Integer,String>> tempPileup=new HashMap<>();
		
		BufferedReader BR=new BufferedReader(new FileReader( pileupFile));
		BufferedWriter BW=new BufferedWriter(new FileWriter(new File(outFile)));
		String line="";
		String[] splitLine;
		while ((line = BR.readLine()) != null){
			if(line.matches(".*\t(\\d*)\t.*")){
//			Chr1    1015    1
//			Chr1    1016    3
				splitLine=line.split("\t");
				if(Double.parseDouble(splitLine[2])>0.0){
					HashMap <Integer,String> temp=null;
					String chrom=splitLine[0];
					int pos=Integer.parseInt(splitLine[1]);
					if(tempPileup.containsKey(chrom)){
						temp=tempPileup.get(chrom);
					}else {
						
						temp=new HashMap<>();
					}
					temp.put(pos, line);
					
					tempPileup.put(chrom, temp);
				}
				
			}		
		}
		BR.close();
		
		BufferedReader FA=new BufferedReader(new FileReader(promotorFasta));
		
		line="";
		
		String gene="";
		String chrom="";
		int startPos=-1;
		int endPos=-1;
		String[] splitHeader;
		String[] splitArea;
		String[] splitPos;
		String[] splitBAM;

		int startBAM=-1;
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
				
				for(int i=startPos;i<endPos;i++){
						if(tempPileup.containsKey(chrom)){
							if(tempPileup.get(chrom).containsKey(i+1)){
								String lineBAM =tempPileup.get(chrom).get(i+1);
								splitBAM=lineBAM.split("\t");
								startBAM=Integer.parseInt(splitBAM[1])-startPos;
								//SRR2981221.167744       83      Chr1    11115   44      38M     =       10942   -211    CGTTTATGTGGCATTAGAATTAAAAATATATGTGGAGC  IIIIGIIIIIGIIIIIIIIGIIIIIIIIIIIIIIIIII  AS:i:76 XS:i:64 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:38 YS:i:76 YT:Z:CP
								String out=gene+"\t"+startBAM+"\t"+splitBAM[2];
								BW.write(out+"\n");
								
							}
							
						}
				}
			}
		}
		FA.close();
		BW.close();	
	}
}
