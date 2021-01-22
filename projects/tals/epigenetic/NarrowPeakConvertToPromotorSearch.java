package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;

public class NarrowPeakConvertToPromotorSearch {
	public static void main(String[] args) throws Exception {
		String narrowPeakFile=args[0];
		String promotorFasta=args[1];
		String outFile=args[2];
		
		HashMap<String, HashMap<Integer,NarrowPeak>> NarrowPeakHash=new HashMap<>();
		
		BufferedReader BR=new BufferedReader(new FileReader( narrowPeakFile));
		BufferedWriter BW=new BufferedWriter(new FileWriter(new File(outFile)));
		String line="";
		String[] splitLine;
		HashMap <Integer,NarrowPeak> temp=null;
		while ((line = BR.readLine()) != null){
			//Chr3	31811403	31817656	Chr3.31811403	10000	.	542552.780333149	302.47173548983	299.811916769022	2401
			//Chr9	14505446	14507122	Chr9.14505446	7984.24190504116	.	433187.264443252	302.47173548983	299.811916769022	248

			splitLine=line.split("\t");
			
			String chrom=splitLine[0];
			Integer startPos=Integer.parseInt(splitLine[1]);
			if(NarrowPeakHash.containsKey(chrom)){
				temp=NarrowPeakHash.get(chrom);
			}else {
				temp=new HashMap<Integer,NarrowPeak>();
			}
			temp.put(startPos, new NarrowPeak(splitLine[0], startPos, Integer.parseInt(splitLine[2]), Float.parseFloat(splitLine[4]), Float.parseFloat(splitLine[6])));
			
			NarrowPeakHash.put(chrom, temp);	
					
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

		NarrowPeak aktNarrowPeak=null;
		String out="";
		while ((line = FA.readLine()) != null){
			if(line.startsWith(">")){
				splitHeader=line.split(" ");
				splitArea=splitHeader[2].split(":");
				
				splitPos=splitArea[1].split("-");
				gene=splitHeader[1];
				
				chrom=splitArea[0];
				
				startPos=Integer.parseInt(splitPos[0]);
				
				endPos=Integer.parseInt(splitPos[1]);

				if(NarrowPeakHash.containsKey(chrom)){
					temp=NarrowPeakHash.get(chrom);
					for(int startPosPeak : temp.keySet()){
						aktNarrowPeak =temp.get(startPosPeak);

							int endPosPeak=aktNarrowPeak.getPeakEndPos();
							
							if((startPosPeak<=endPos)&(endPosPeak>=startPos)){
								if((startPosPeak>=startPos)&(endPosPeak<=endPos)){
									out=gene+"\t"+(startPosPeak-startPos)+"\t"+(endPosPeak-startPos)+"\t"+"."+"\t"+aktNarrowPeak.getPeakScore()+"\t"+"."+"\t"+aktNarrowPeak.getPeakValue();
									BW.write(out+"\n");
								}else if((startPosPeak<=startPos)&(endPosPeak<=endPos)){
									out=gene+"\t"+0+"\t"+(endPosPeak-startPos)+"\t"+"."+"\t"+aktNarrowPeak.getPeakScore()+"\t"+"."+"\t"+aktNarrowPeak.getPeakValue();
									BW.write(out+"\n");
								}else if((startPosPeak>=startPos)&(endPosPeak>=endPos)){
									out=gene+"\t"+(startPosPeak-startPos)+"\t"+(endPos-startPos)+"\t"+"."+"\t"+aktNarrowPeak.getPeakScore()+"\t"+"."+"\t"+aktNarrowPeak.getPeakValue();
									BW.write(out+"\n");
								}else if((startPosPeak<=startPos)&(endPosPeak>=endPos)){
									out=gene+"\t"+0+"\t"+(endPos-startPos)+"\t"+"."+"\t"+aktNarrowPeak.getPeakScore()+"\t"+"."+"\t"+aktNarrowPeak.getPeakValue();
									BW.write(out+"\n");
								}
							}
					}
					
				}
			
			}
		}
		FA.close();
		BW.close();
	}
}
