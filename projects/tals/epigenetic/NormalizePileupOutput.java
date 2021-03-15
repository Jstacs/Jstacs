package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

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
		FileParameter pileupFile = new FileParameter("pileup-output-File","Pileup output file.","tsv.gz,tsv",true);
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
		
		HashMap<String, HashMap<Integer,String>> tempPileup=new HashMap<>();
		HashMap<String, HashMap<Integer,Integer>> tempPileupCov=new HashMap<>();
		
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
					os_ps.print(chrom+"\t"+i+"\t"+normalizeCov[i]+"\n");
				}
			}
		}
		os.close();
		
		TextResult tr = new TextResult("Normalized pileup file", "Normalized pileup file", new FileParameter.FileRepresentation(out.getAbsolutePath()), "tsv.gz", getToolName(), null, true);
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );

	}

	@Override
	public String getToolName() {
		return "NormalizePielupOutput";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "NormalizePielupOutput";
	}

	@Override
	public String getDescription() {
		return "Normalizes pileup output";
	}

	@Override
	public String getHelpText() {
		// TODO Auto-generated method stub
		return null;
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
