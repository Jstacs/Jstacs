package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;

public class NormalizePielupOutput {
	public static void main(String[] args) throws Exception {
		String pileupFile=args[0];
		String outFile=args[1];
		
		HashMap<String, HashMap<Integer,String>> tempPileup=new HashMap<>();
		HashMap<String, HashMap<Integer,Integer>> tempPileupCov=new HashMap<>();
		
		BufferedReader BR=new BufferedReader(new FileReader( pileupFile));
		BufferedWriter BW=new BufferedWriter(new FileWriter(new File(outFile)));
		String line="";
		String[] splitLine;
		while ((line = BR.readLine()) != null){
			if(line.matches("Chr([\\d]+)\\t([\\d]+)\\t([\\d]+)")){
//				Chr1    1015    1
//				Chr1    1016    3
				splitLine=line.split("\t");
				
				HashMap <Integer,String> temp=null;
				HashMap <Integer,Integer> tempCov=null;
				String chrom=splitLine[0];

				int pos=Integer.parseInt(splitLine[1]);
				if(tempPileup.containsKey(chrom)){
					temp=tempPileup.get(chrom);
					tempCov=tempPileupCov.get(chrom);
				}else {
					temp=new HashMap<>();
					tempCov=new HashMap<>();
				}
				temp.put(pos, line);
				tempCov.put(pos, Integer.parseInt(splitLine[2]));
				
				tempPileup.put(chrom, temp);	
				tempPileupCov.put(chrom, tempCov);	

			}

		}
		BR.close();

		double window=10000.0;
		int half=(int)(window/2);
		
		
		for (String chrom : tempPileupCov.keySet()) {
			
			
			Integer[] tempPileupCovKeySetArray = new Integer[tempPileupCov.get(chrom).keySet().size()];
			tempPileupCov.get(chrom).keySet().toArray(tempPileupCovKeySetArray);
			Arrays.sort(tempPileupCovKeySetArray);
			Integer lastPos=tempPileupCovKeySetArray[tempPileupCovKeySetArray.length-1];
			Double[] originalCov=new Double[lastPos];
			Double[] normalizeCov=new Double[lastPos];
			Arrays.fill(originalCov, 0.0);
			Arrays.fill(normalizeCov, 0.0);
			
			for(int i=0;i<originalCov.length;i++){
				if(tempPileupCov.get(chrom).containsKey(i)){
					originalCov[i]=tempPileupCov.get(chrom).get(i).doubleValue();
				}
			}
			double tempSum=0;
			boolean isFirst=true;
			for(int i=0;i<originalCov.length;i++){
				int windowStart=((i-half<0)?0:(i-half));
				int windowend=((i+half>lastPos)?lastPos:(i+half));
				if(isFirst){
					for(int j=windowStart;j<windowend;j++){
						tempSum+=originalCov[j];
					}
					isFirst=false;
				}else{
					tempSum-=originalCov[windowStart];
					tempSum+=originalCov[windowend-1];
				}

				normalizeCov[i]=originalCov[i]-(tempSum/(windowend-windowStart+1));

				if(tempPileup.get(chrom).containsKey(i)){
					BW.write(chrom+"\t"+i+"\t"+normalizeCov[i]+"\n");
				}
			}
		}
		BW.close();
	}
}
