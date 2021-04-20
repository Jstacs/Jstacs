package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

public class NarrowpeakprofilHashMap {

	private HashMap<String, Narrowpeakprofil> narrowpeakProfiles=new HashMap<>();
	
	public NarrowpeakprofilHashMap(String faiPath,String PathtoFile,String type) throws Exception {
		this(faiPath,PathtoFile,300,50,type);
	}
	
	public NarrowpeakprofilHashMap(HashMap<String, Integer> chromLengthHash,String PathtoFile, int before,int after) throws Exception {
		BufferedReader BR;
		
		if(PathtoFile.endsWith("gz")){
			BR=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(PathtoFile)))));
		}else{
			BR=new BufferedReader(new InputStreamReader(new FileInputStream(new File(PathtoFile))));
		}
		
		String line="";
		HashMap <String,float[]>peakValueHash=new HashMap<>();
		float[] peakValue;
		HashMap <String,boolean[][]>peakBoolHash=new HashMap<>();
		boolean[][] peakBool;
		for(String chrom :chromLengthHash.keySet() ){
			peakBool=new boolean[2][chromLengthHash.get(chrom)];
			peakValue=new float[chromLengthHash.get(chrom)];
			Arrays.fill(peakBool[0],false);
			Arrays.fill(peakBool[1],false);
			Arrays.fill(peakValue,0.0f);
			peakBoolHash.put(chrom, peakBool.clone());
			peakValueHash.put(chrom, peakValue.clone());
		}
		
		while ((line = BR.readLine()) != null){
			String[] splitLine2=line.split("\t");
			String chrom=splitLine2[0];
		
			peakValue=peakValueHash.get(chrom);
			
			int startPeakPos=Integer.parseInt(splitLine2[1]);
			int endPeakPos=Integer.parseInt(splitLine2[2])-1;
			float aktPeak=Float.parseFloat(splitLine2[6]);
			for(int p=startPeakPos;p<=endPeakPos;p++){
				peakValue[p]=aktPeak;
			}	
			
			peakValueHash.put(chrom,peakValue);
		}
		BR.close();
		for(String chrom :chromLengthHash.keySet() ){
			peakBool=peakBoolHash.get(chrom);
			peakValue=peakValueHash.get(chrom);
			int chromLength=chromLengthHash.get(chrom);
			int lastPos=chromLengthHash.get(chrom)-1;
			for(int i=0;i<chromLength;i++){
				int windowStart=((i-before<0)?0:(i-before));
				int windowend=((i+after>lastPos)?lastPos:(i+after));
				for(int s=windowStart;s<=windowend;s++){
					if(peakValue[s]>0.0f){
						peakBool[0][i]=true;
						break;
					}
				}	
			}
			
			for(int i=0;i<chromLength;i++){
				int windowStart=((i-after<0)?0:(i-after));
				int windowend=((i+before>lastPos)?lastPos:(i+before));
				for(int s=windowStart;s<=windowend;s++){
					if(peakValue[s]>0.0f){
						peakBool[1][i]=true;
						break;
					}
				}	
			}
			peakBoolHash.put(chrom, peakBool);
			narrowpeakProfiles.put(chrom, new Narrowpeakprofil(chrom, chromLengthHash.get(chrom), peakBoolHash.get(chrom), before,after));
		}
	}
	
	public NarrowpeakprofilHashMap(String faiPath,String PathtoFile, int before,int after,String type) throws Exception {
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
		HashMap <String,float[]>peakValueHash=new HashMap<>();
		float[] peakValue;
		HashMap <String,boolean[][]>peakBoolHash=new HashMap<>();
		boolean[][] peakBool;
		for(String chrom :chromLengthHash.keySet() ){
			peakBool=new boolean[2][chromLengthHash.get(chrom)];
			peakValue=new float[chromLengthHash.get(chrom)];
			Arrays.fill(peakBool[0],false);
			Arrays.fill(peakBool[1],false);
			Arrays.fill(peakValue,0.0f);
			peakBoolHash.put(chrom, peakBool.clone());
			peakValueHash.put(chrom, peakValue.clone());
		}
		
		boolean typeIsGenome=false;
		if(type.equals("genome")){
			typeIsGenome=true;
		}else if (!type.equals("promotor")) {
			System.err.println("Type should be 'genome' or 'promotor'!");
		}
		
		while ((line = BR.readLine()) != null){
			String[] splitLine2=line.split("\t");
			String chrom=splitLine2[0];
			peakValue=peakValueHash.get(chrom);
			if(typeIsGenome){
				int startPeakPos=Integer.parseInt(splitLine2[1]);
				int endPeakPos=Integer.parseInt(splitLine2[2])-1;
				float aktPeak=Float.parseFloat(splitLine2[6]);
				for(int p=startPeakPos;p<=endPeakPos;p++){
					peakValue[p]=aktPeak;
				}
			}else{
				peakValue[Integer.parseInt(splitLine2[1])-1]=Float.parseFloat(splitLine2[2]);
			}
			peakValueHash.put(chrom,peakValue);
		}
		BR.close();
		for(String chrom :chromLengthHash.keySet() ){
			peakBool=peakBoolHash.get(chrom);
			peakValue=peakValueHash.get(chrom);
			chromLength=chromLengthHash.get(chrom);
			int lastPos=chromLengthHash.get(chrom)-1;
			for(int i=0;i<chromLength;i++){
				int windowStart=((i-before<0)?0:(i-before));
				int windowend=((i+after>lastPos)?lastPos:(i+after));
				for(int s=windowStart;s<=windowend;s++){
					if(peakValue[s]>0.0f){
						peakBool[0][i]=true;
						break;
					}
				}	
			}
			
			for(int i=0;i<chromLength;i++){
				int windowStart=((i-after<0)?0:(i-after));
				int windowend=((i+before>lastPos)?lastPos:(i+before));
				for(int s=windowStart;s<=windowend;s++){
					if(peakValue[s]>0.0f){
						peakBool[1][i]=true;
						break;
					}
				}	
			}
			peakBoolHash.put(chrom, peakBool);
			narrowpeakProfiles.put(chrom, new Narrowpeakprofil(chrom, chromLengthHash.get(chrom), peakBoolHash.get(chrom), before,after));
		}
	}
	
	public Narrowpeakprofil getNarrowpeakprofil(String chrom){
		return this.narrowpeakProfiles.get(chrom);
	}
}
