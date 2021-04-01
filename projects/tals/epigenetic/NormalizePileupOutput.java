package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import de.jstacs.io.FileManager;
import de.jstacs.parameters.FileParameter;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;

public class NormalizePileupOutput implements JstacsTool{
	
	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new NormalizePileupOutput());
		
		cli.run(args);
	}
	
	public NormalizePileupOutput() {

	}

	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter pileupFile = new FileParameter("pileup-output-File","Pileup output file.","tsv.gz,tsv,txt",true);
		return new ToolParameterSet(this.getShortName(),pileupFile);
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol,
			ProgressUpdater progress, int threads) throws Exception {
		progress.setLast(1.0);
		progress.setCurrent(0.0);
		String pileupFile = parameters.getParameterAt(0).getValue().toString();
		BufferedReader BR=null;
		if(pileupFile.endsWith("gz")){
			BR=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(pileupFile)))));
		}else{
			BR=new BufferedReader(new InputStreamReader(new FileInputStream(new File(pileupFile))));
		}
		
		File out = File.createTempFile("pileup.normalized", ".temp.tsv.gz", new File("."));
		out.deleteOnExit();
		
		GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(out));
		PrintStream os_ps=new PrintStream(os);
		
		//HashMap<String, HashMap<Integer,String>> tempPileup=new HashMap<>();
		
		String line="";
		String[] splitLine;
	//	HashMap <Integer,String> temp=null;
		HashMap <Integer,Integer> tempCov=null;
		String chrom_before="";
		String chrom="";
		int pos=0;
		
		boolean first=true;
		while ((line = BR.readLine()) != null){
			if(line.matches("Chr([\\d]+)\\t([\\d]+)\\t([\\d]+)")){
//				Chr1    1015    1
//				Chr1    1016    3
				splitLine=line.split("\t");
				chrom=splitLine[0];
				pos=Integer.parseInt(splitLine[1]);
				if(chrom.equals(chrom_before)){
					tempCov.put(pos, Integer.parseInt(splitLine[2]));
				}else{//neues chromosom beginnt
					if(first){
						first=false;
					}else{
						//System.out.println(line);
						normalizeChrom(chrom_before,os_ps,tempCov);
					}
					
					tempCov=new HashMap<>();
					tempCov.put(pos, Integer.parseInt(splitLine[2]));
				}
				chrom_before=chrom;
			}
		}
		normalizeChrom(chrom_before,os_ps,tempCov);
		BR.close();


	
		os.close();
		
		TextResult tr = new TextResult("Normalized pileup file", "Normalized pileup file", new FileParameter.FileRepresentation(out.getAbsolutePath()), "tsv.gz", getToolName(), null, true);
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );

	}
	
	private void normalizeChrom(String chrom,PrintStream os_ps,HashMap <Integer,Integer> tempCov){
		double window=10000.0;
		int half=(int)(window/2);
		Integer[] tempPileupCovKeySetArray= new Integer[tempCov.keySet().size()];
		tempCov.keySet().toArray(tempPileupCovKeySetArray);
		Arrays.sort(tempPileupCovKeySetArray);
		Integer lastPos=tempPileupCovKeySetArray[tempPileupCovKeySetArray.length-1];
		Double[] originalCov=new Double[lastPos];
		Double[] normalizeCov=new Double[lastPos];
		Arrays.fill(originalCov, 0.0);
		Arrays.fill(normalizeCov, 0.0);
		
		for(int i=0;i<originalCov.length;i++){
			if(tempCov.containsKey(i)){
				originalCov[i]=tempCov.get(i).doubleValue();
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

			if(tempCov.containsKey(i)){
				if(normalizeCov[i]>0.0){
					os_ps.print(chrom+"\t"+i+"\t"+normalizeCov[i]+"\n");
				}
			}
		}
	}

	@Override
	public String getToolName() {
		return "NormalizePileupOutput";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "normpileup";
	}

	@Override
	public String getDescription() {
		return "Normalizes pileup output";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( NormalizePileupOutput.class.getClassLoader().getResourceAsStream( "projects/tals/epigenetic/toolHelpFiles/NormalizePileupOutput.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void clear() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String[] getReferences() {
		// TODO Auto-generated method stub
		return null;
	}
}
