package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;

public class PileupCoverageprofilHashMap {

	private HashMap<String, PileupCoverageprofil> pileupCoverageProfiles=new HashMap<>();
	
	public PileupCoverageprofilHashMap(String faiPath,String PathtoFile) throws Exception {
		this(faiPath,PathtoFile,300,50,false);
	}
	
	public PileupCoverageprofilHashMap(String faiPath,String PathtoFile, int before,int after,boolean calculateAlwaysOnCompleteSeq) throws Exception {
		HashMap<String, Integer> chromLengthHash=new HashMap<>();
		
		BufferedReader readFai = new BufferedReader(new FileReader(faiPath));
		String line="";
		String[] splitLine;
		int chromLength=-1;
		while ((line = readFai.readLine()) != null){

			splitLine=line.split("\t");
			chromLength=Integer.parseInt(splitLine[1]);
			chromLengthHash.put(splitLine[0],chromLength);
			
		}
		readFai.close();
		
		BufferedReader BR=new BufferedReader(new InputStreamReader(new FileInputStream(new File(PathtoFile))));
		line="";
		HashMap <String,float[]>covValueHash=new HashMap<>();
		float[] covValue;
		HashMap <String,int[][]>numberOfCoveragePositionsHash=new HashMap<>();
		int[][] numberOfCoveragePositions;
		for(String chrom :chromLengthHash.keySet() ){
			numberOfCoveragePositions=new int[2][chromLengthHash.get(chrom)];
			covValue=new float[chromLengthHash.get(chrom)];
			Arrays.fill(numberOfCoveragePositions[0],0);
			Arrays.fill(numberOfCoveragePositions[1],0);
			Arrays.fill(covValue,0.0f);
			numberOfCoveragePositionsHash.put(chrom, numberOfCoveragePositions.clone());
			covValueHash.put(chrom, covValue.clone());
		}
		
		while ((line = BR.readLine()) != null){
			String[] splitLine2=line.split("\t");
			String chrom=splitLine2[0];
			covValue=covValueHash.get(chrom);
			covValue[Integer.parseInt(splitLine2[1])-1]=Float.parseFloat(splitLine2[2]);
			covValueHash.put(chrom,covValue);
		}
		BR.close();
		for(String chrom :chromLengthHash.keySet() ){
			numberOfCoveragePositions=numberOfCoveragePositionsHash.get(chrom);
			covValue=covValueHash.get(chrom);
			chromLength=chromLengthHash.get(chrom);
			//calculate number of coverage positions on complete promotor
			if(calculateAlwaysOnCompleteSeq){
				int tempNumberOfCoveragePositions=0;
				for(int i=0;i<chromLength;i++){
					if(covValue[i]>0.0f){
						tempNumberOfCoveragePositions++;
					}
				}
				Arrays.fill(numberOfCoveragePositions[0], tempNumberOfCoveragePositions);
				Arrays.fill(numberOfCoveragePositions[1], tempNumberOfCoveragePositions);
				
			}else{ //calculate number of coverage positions in defined area surround binding side
				int lastPos=chromLengthHash.get(chrom)-1;
				boolean isFirst=true;
				for(int i=0;i<chromLength;i++){
					int windowStart=((i-before<0)?0:(i-before));
					int windowend=((i+after>lastPos)?lastPos:(i+after));
					
					if(isFirst){
						for(int s=windowStart;s<=windowend;s++){
							if(covValue[s]>0.0f){
								numberOfCoveragePositions[0][i]++;
							}
						}
					}else{
						if(covValue[windowStart]>0.0f){
							numberOfCoveragePositions[0][i]--;
						}
						if(covValue[windowend-1]>0.0f){
							numberOfCoveragePositions[0][i]++;
						}
					}	
				}
				
				isFirst=true;
				for(int i=0;i<chromLength;i++){
					int windowStart=((i-after<0)?0:(i-after));
					int windowend=((i+before>lastPos)?lastPos:(i+before));
					
					if(isFirst){
						for(int s=windowStart;s<=windowend;s++){
							if(covValue[s]>0.0f){
								numberOfCoveragePositions[1][i]++;
							}
						}
					}else{
						if(covValue[windowStart]>0.0f){
							numberOfCoveragePositions[1][i]--;
						}
						if(covValue[windowend-1]>0.0f){
							numberOfCoveragePositions[1][i]++;
						}
					}	
				}
			}
			
			numberOfCoveragePositionsHash.put(chrom, numberOfCoveragePositions);
			pileupCoverageProfiles.put(chrom, new PileupCoverageprofil(chrom, chromLengthHash.get(chrom), numberOfCoveragePositionsHash.get(chrom), before,after,calculateAlwaysOnCompleteSeq));
		}
		
	}
	
	public PileupCoverageprofil getPileupCoverageprofil(String chrom){
		return this.pileupCoverageProfiles.get(chrom);
	}
	
}
