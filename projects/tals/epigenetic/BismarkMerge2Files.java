package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
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

public class BismarkMerge2Files implements JstacsTool{
	
	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new BismarkMerge2Files());
		
		cli.run(args);
	}
	
	public BismarkMerge2Files() {

	}

	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter bismarkFile1 = new FileParameter("bismark-File-1","Methylationinformation in bismark format file 1","cov.gz,cov",true);
		FileParameter bismarkFile2 = new FileParameter("bismark-File-2","Methylationinformation in bismark format file 2","cov.gz,cov",true);
		return new ToolParameterSet(this.getShortName(),bismarkFile1,bismarkFile2);
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol,
			ProgressUpdater progress, int threads) throws Exception {
		progress.setLast(1.0);
		progress.setCurrent(0.0);
		String bismarkFile1 = parameters.getParameterAt(0).getValue().toString();
		BufferedReader BR=null;
		if(bismarkFile1.endsWith("gz")){
			BR=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(bismarkFile1)))));
		}else{
			BR=new BufferedReader(new InputStreamReader(new FileInputStream(new File(bismarkFile1))));
		}
		
		String bismarkFile2 = parameters.getParameterAt(1).getValue().toString();
		BufferedReader BR2=null;
		if(bismarkFile2.endsWith("gz")){
			BR2=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(bismarkFile2)))));
		}else{
			BR2=new BufferedReader(new InputStreamReader(new FileInputStream(new File(bismarkFile2))));
		}
		
		File out = File.createTempFile("merged.bismark", ".temp.gz", new File("."));
		out.deleteOnExit();
		
		GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(out));
		PrintStream os_ps=new PrintStream(os);
		
		HashMap<String, HashMap<Integer,String>> tempBismark=new HashMap<>();
		
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
		
		for(String chrom : tempBismark.keySet()){
			HashMap <Integer,String> temp=tempBismark.get(chrom);
			List<Integer> sortedByKey = new ArrayList<>(temp.keySet());
			Collections.sort(sortedByKey);
			for (Integer i : sortedByKey) {
				os_ps.print(temp.get(i)+"\n");
			}			
		}
		os.close();
		
		TextResult tr = new TextResult("Merged Bismark files", "Merged Bismark files", new FileParameter.FileRepresentation(out.getAbsolutePath()), "cov.gz", getToolName(), null, true);
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );

	}

	@Override
	public String getToolName() {
		return "BismarkFileMerger";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "bismarkfilemerger";
	}

	@Override
	public String getDescription() {
		return "Merges 2 Files in bismark format";
	}
	
	@Override
	public String getHelpText() {
		String helpText="Merges files generated by Bismark methylation extractor with "
				+ "parameters ”–bedGraph –CX -p”. In plants, cytosine methylation "
				+ "occurs in the following three contexts: ’CpG’, ’CpHpG’ and ’CpHpH’ "
				+ "(H = ’A’, ’C’ or ’T’). With option ”–CX”, the output contains the methylation of "
				+ "cytosines in all three contexts. The output contains a coverage file, which "
				+ "contains the columns: <chromosome> <start position> <end position> <methylation "
				+ "percentage> <count methylated> <count unmethylated>";
		return helpText;
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
