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

public class PileupConvertToPromoterSearch implements JstacsTool{
	
	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new PileupConvertToPromoterSearch());
		
		cli.run(args);
	}
	
	public PileupConvertToPromoterSearch() {

	}

	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter pileupFile = new FileParameter("Normalized pileup output file","Normalized pileup output file.","tsv.gz,tsv",true);
		FileParameter promoterFasta = new FileParameter("Promoter fasta file","Promoter fastA file","fa,fasta",true);
		return new ToolParameterSet(this.getShortName(),pileupFile,promoterFasta);
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
				String promotorFasta = parameters.getParameterAt(1).getValue().toString();
				BufferedReader FA=new BufferedReader(new FileReader(promotorFasta));
				
				File out = File.createTempFile("promoter.pileup", ".temp.tsv.gz", new File("."));
				out.deleteOnExit();
				
				GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(out));
				PrintStream os_ps=new PrintStream(os);
				
				HashMap<String, HashMap<Integer,String>> tempPileup=new HashMap<>();
				
				String line="";
				String[] splitLine;
				while ((line = BR.readLine()) != null){
					if(line.matches(".*\t(\\d*)\t.*")){
//					Chr1    1015    1
//					Chr1    1016    3
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
//				boolean strand;
				while ((line = FA.readLine()) != null){
					if(line.startsWith(">")){
						splitHeader=line.split(" ");
						splitArea=splitHeader[2].split(":");
//						if(splitArea[2].equals("+")){
//							strand=true;
//						}else{
//							strand=false;
//						}
						splitPos=splitArea[1].split("-");
						gene=line.substring(1).trim();
						int idx = gene.indexOf(" ");
						if (idx > 0) {
							gene = gene.substring(0, idx);
						}
						
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
										String output=gene+"\t"+startBAM+"\t"+splitBAM[2];
										os_ps.print(output+"\n");
										
									}
									
								}
						}
					}
				}
				FA.close();
				os.close();	
								
				TextResult tr = new TextResult("Pileup promoter file", "Pileup promoter file", new FileParameter.FileRepresentation(out.getAbsolutePath()), "tsv.gz", getToolName(), null, true);
				return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
	}

	@Override
	public String getToolName() {
		return "PileupConvertToPromoter";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "pile2prom";
	}

	@Override
	public String getDescription() {
		return "Creates Pileup file in promoter region";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( PileupConvertToPromoterSearch.class.getClassLoader().getResourceAsStream( "projects/tals/epigenetic/toolHelpFiles/PileupConvertToPromoterSearch.txt" ) ).toString();
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
