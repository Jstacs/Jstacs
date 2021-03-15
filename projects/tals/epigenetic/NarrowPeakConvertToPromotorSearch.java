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

public class NarrowPeakConvertToPromotorSearch implements JstacsTool{
	
	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new NarrowPeakConvertToPromotorSearch());
		
		cli.run(args);
	}
	
	public NarrowPeakConvertToPromotorSearch() {

	}

	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter narrowPeakFile = new FileParameter("bismark-File-1","Methylationinformation in bismark format file 1","cov.gz,cov",true);
		FileParameter promotorFasta = new FileParameter("promotor fasta file","Promotor fastA file","fa,fasta",true);
		return new ToolParameterSet(this.getShortName(),narrowPeakFile,promotorFasta);
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol,
			ProgressUpdater progress, int threads) throws Exception {
		progress.setLast(1.0);
		progress.setCurrent(0.0);
		
		String narrowPeakFile = parameters.getParameterAt(0).getValue().toString();
		BufferedReader BR=null;
		if(narrowPeakFile.endsWith("gz")){
			BR=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(narrowPeakFile)))));
		}else{
			BR=new BufferedReader(new InputStreamReader(new FileInputStream(new File(narrowPeakFile))));
		}
		
		String promotorFasta = parameters.getParameterAt(1).getValue().toString();
		BufferedReader FA=new BufferedReader(new FileReader(promotorFasta));
		
		File outF = File.createTempFile("promotor.peaks.narrowPeak", ".temp.gz", new File("."));
		outF.deleteOnExit();
		
		GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(outF));
		PrintStream os_ps=new PrintStream(os);
		
		HashMap<String, HashMap<Integer,NarrowPeak>> NarrowPeakHash=new HashMap<>();

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
				gene=line.substring(1).trim();
				int idx = gene.indexOf(" ");
				if (idx > 0) {
					gene = gene.substring(0, idx);
				}
				
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
									os_ps.print(out+"\n");
								}else if((startPosPeak<=startPos)&(endPosPeak<=endPos)){
									out=gene+"\t"+0+"\t"+(endPosPeak-startPos)+"\t"+"."+"\t"+aktNarrowPeak.getPeakScore()+"\t"+"."+"\t"+aktNarrowPeak.getPeakValue();
									os_ps.print(out+"\n");
								}else if((startPosPeak>=startPos)&(endPosPeak>=endPos)){
									out=gene+"\t"+(startPosPeak-startPos)+"\t"+(endPos-startPos)+"\t"+"."+"\t"+aktNarrowPeak.getPeakScore()+"\t"+"."+"\t"+aktNarrowPeak.getPeakValue();
									os_ps.print(out+"\n");
								}else if((startPosPeak<=startPos)&(endPosPeak>=endPos)){
									out=gene+"\t"+0+"\t"+(endPos-startPos)+"\t"+"."+"\t"+aktNarrowPeak.getPeakScore()+"\t"+"."+"\t"+aktNarrowPeak.getPeakValue();
									os_ps.print(out+"\n");
								}
							}
					}
					
				}
			
			}
		}
		FA.close();
		os.close();
		
		TextResult tr = new TextResult("Narrow peak promotor file", "Narrow peak promotor file", new FileParameter.FileRepresentation(outF.getAbsolutePath()), "narrowPeak.gz", getToolName(), null, true);
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );

	}

	@Override
	public String getToolName() {
		return "NarrowPeakConvertToPromotor";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "NarrowPeakConvertToPromotor";
	}

	@Override
	public String getDescription() {
		return "Creates NarrowPeak file in promotor region";
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
