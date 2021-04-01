package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
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


public class BismarkConvertToPromoterSearch implements JstacsTool{
	
	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new BismarkConvertToPromoterSearch());
		
		cli.run(args);
	}
	
	public BismarkConvertToPromoterSearch() {

	}

	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter bismarkFile = new FileParameter("bismark-File","Methylationinformation in bismark format","cov.gz,cov",true);
		FileParameter promotorFasta = new FileParameter("promoter fasta file","Promoter fastA file","fa,fasta",true);
		return new ToolParameterSet(this.getShortName(),bismarkFile,promotorFasta);
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol,
			ProgressUpdater progress, int threads) throws Exception {
		progress.setLast(1.0);
		progress.setCurrent(0.0);
		
		String bismarkFile = parameters.getParameterAt(0).getValue().toString();
		BufferedReader BR=null;
		if(bismarkFile.endsWith("gz")){
			BR=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(bismarkFile)))));
		}else{
			BR=new BufferedReader(new InputStreamReader(new FileInputStream(new File(bismarkFile))));
		}
		
		String promotorFasta = parameters.getParameterAt(1).getValue().toString();
		BufferedReader FA=new BufferedReader(new FileReader(promotorFasta));
		
		File out = File.createTempFile("bimark.promoter", ".temp.cov.gz", new File("."));
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
				//System.out.println(splitLine[0]);
				
				temp=new HashMap<>();
			}
			temp.put(Integer.parseInt(splitLine[1]), line);
			
			tempBismark.put(splitLine[0], temp);
		}
		BR.close();
		
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
				gene=line.substring(1).trim();
				int idx = gene.indexOf(" ");
				if (idx > 0) {
					gene = gene.substring(0, idx);
				}
				
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
								os_ps.print(gene+"\t"+startBismark+"\t"+startBismark+"\t"+splitBismark[3]+"\t"+splitBismark[4]+"\t"+splitBismark[5]+"\n");
							}
							
						}
					}
					if(strand==false){
						pos=endPos-k;
						
						if(tempBismark.containsKey(chrom)){
							if(tempBismark.get(chrom).containsKey(pos)){
								splitBismark=tempBismark.get(chrom).get(pos).split("\t");
								startBismark=k+1;
								os_ps.print(gene+"\t"+startBismark+"\t"+startBismark+"\t"+splitBismark[3]+"\t"+splitBismark[4]+"\t"+splitBismark[5]+"\n");
							}
						}
					}
					k++;
				}
			}
		}
		FA.close();
		os.close();
		
		TextResult tr = new TextResult("Bismark promoter file", "Bismark file in promoter region", new FileParameter.FileRepresentation(out.getAbsolutePath()), "cov.gz", getToolName(), null, true);
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );

	}

	@Override
	public String getToolName() {
		return "BismarkConvertToPromoter";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "bis2prom";
	}

	@Override
	public String getDescription() {
		return "Creates Bismark file in promoter region";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( BismarkConvertToPromoterSearch.class.getClassLoader().getResourceAsStream( "projects/tals/epigenetic/toolHelpFiles/BismarkConvertToPromoterSearch.txt" ) ).toString();
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
