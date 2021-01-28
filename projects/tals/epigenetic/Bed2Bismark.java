package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Date;
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

public class Bed2Bismark implements JstacsTool {

	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new Bed2Bismark());
		
		cli.run(args);

	}
	
	public Bed2Bismark() {

	}
	
	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter bedMethylFile = new FileParameter("bedMethyl-file","Methylationinformation in bedMethyl format","bed.gz,bed",true);
		return new ToolParameterSet(this.getShortName(),bedMethylFile);
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol,
			ProgressUpdater progress, int threads) throws Exception {
		progress.setLast(1.0);
		progress.setCurrent(0.0);
		String bedMethylFile = parameters.getParameterAt(0).getValue().toString();
		BufferedReader BR=null;
		if(bedMethylFile.endsWith("gz")){
			BR=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(bedMethylFile)))));
		}else{
			BR=new BufferedReader(new InputStreamReader(new FileInputStream(new File(bedMethylFile))));
		}
		
		File out = File.createTempFile("bismark", ".temp.gz", new File("."));
		out.deleteOnExit();

		GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(out));
		PrintStream os_ps=new PrintStream(os);
		String line="";
		String[] splitLine;
		String sep="\t";
		
		while ((line = BR.readLine()) != null){
				splitLine=line.split("\t");
				double methylation_level=Double.parseDouble(splitLine[10])/100.0;
				os_ps.print(splitLine[0]+sep+splitLine[1]+sep+splitLine[2]+sep+splitLine[10]+sep+(int)Math.round(Double.parseDouble(splitLine[9])*methylation_level)+sep+(int)Math.round(Double.parseDouble(splitLine[9])*(1.0-methylation_level))+"\n");
		}
		BR.close();
		os.close();
		
		TextResult tr = new TextResult("Converted Bismark file", "Bismark file converted from bed methyl file", new FileParameter.FileRepresentation(out.getAbsolutePath()), "cov.gz", getToolName(), null, true);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
	}

	@Override
	public String getToolName() {
		return "BedMethylConverter";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "bedmethylconverter";
	}

	@Override
	public String getDescription() {
		return "converts bedMethyl files to bismark";
	}


	@Override
	public String getHelpText() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
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
