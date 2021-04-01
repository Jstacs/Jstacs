package projects.tals.epigenetic;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Date;
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
import projects.encodedream.ObjectStream;
import projects.encodedream.Pileup;
import projects.encodedream.Pileup.Pile;

public class PileupTool implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		
		FileParameter fp = new FileParameter("BAM file", "Mapped reads from DNase-seq or ATAC-seq experiment", "bam", true);
		
		ToolParameterSet tps = new ToolParameterSet(getToolName(), fp);
		return tps;
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		String bam = ((FileParameter)parameters.getParameterAt(0)).getFileContents().getFilename();
		ObjectStream<Pile> ps = new ObjectStream<>(10000);
		
		File out = File.createTempFile("pileup", ".temp.gz", new File("."));
		out.deleteOnExit();

		GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(out));
		
		PrintStream print = new PrintStream(os);
		
		new Thread( ()->{
			try {
				Pileup.pileup(bam, ps, false, true, true);
				ps.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}).start();
		
		ps.print(print);
		os.close();
		
		TextResult tr = new TextResult("Pileup", "Pileup of 5' ends of reads", new FileParameter.FileRepresentation(out.getAbsolutePath()), "tsv.gz", getToolName(), null, true);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
	}

	@Override
	public String getToolName() {
		return "Chromatin pileup";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "pileup";
	}

	@Override
	public String getDescription() {
		return "computes coverage pileup from BAM";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( PileupTool.class.getClassLoader().getResourceAsStream( "projects/tals/epigenetic/toolHelpFiles/PileupTool.txt" ) ).toString();
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
		return null;
	}

	@Override
	public void clear() {

	}

	@Override
	public String[] getReferences() {
		return null;
	}

}
