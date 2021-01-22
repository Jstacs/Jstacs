package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

public class MethylationprofilHashMap {

	private HashMap<String, Methylationprofil> methylationProfiles=new HashMap<>();
	
	public MethylationprofilHashMap(String faiPath,String pathToBismarkFile) throws Exception {
		this(faiPath,pathToBismarkFile,0f,0f);
	}
	
	public MethylationprofilHashMap(String faiPath,String pathToBismarkFile, float PseudoCounts, float probCMethylated) throws Exception {
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
	//	methylationProfiles =new
		
		BufferedReader BR=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(pathToBismarkFile)))));
		line="";
		HashMap <String,float[]>MethylationProbHash=new HashMap<>();
		float[] MethylationProb;
		for(String chrom :chromLengthHash.keySet() ){
			MethylationProb=new float[chromLengthHash.get(chrom)];
			Arrays.fill(MethylationProb,PseudoCounts*probCMethylated);
			MethylationProbHash.put(chrom, MethylationProb.clone());
		}
		while ((line = BR.readLine()) != null){
			
			//<chromosome>	<start position>	<end position>	<methylation percentage>	<count methylated>	<count unmethylated>
			//chromosome02    359     359     83.3333333333333        5       1
			//CM/(CM+CU)
			
			float count_methylated=0f;
			float count_unmethylated=0f;
				String[] splitLine2=line.split("\t");
				String chrom=splitLine2[0];
				
				if(!splitLine2[1].equals(splitLine2[2])){
					BR.close();
					throw new Exception("Different <start position> and <end position> in "+pathToBismarkFile);
				}
				
				count_methylated=Float.parseFloat(splitLine2[4])+PseudoCounts*probCMethylated;
				count_unmethylated=Float.parseFloat(splitLine2[5])+PseudoCounts*(1-probCMethylated);
				MethylationProb=MethylationProbHash.get(chrom);
				MethylationProb[Integer.parseInt(splitLine2[1])-1]=count_methylated/(count_methylated+count_unmethylated);
				MethylationProbHash.put(chrom, MethylationProb);	
				
		}
		BR.close();
		for(String chrom : MethylationProbHash.keySet() ){
			methylationProfiles.put(chrom, new Methylationprofil(0, chrom, chromLengthHash.get(chrom), MethylationProbHash.get(chrom), true, PseudoCounts,probCMethylated));
		}
	}
	
	public Methylationprofil getMethylationprofil(String chrom){
		return this.methylationProfiles.get(chrom);
	}
	
	
}
