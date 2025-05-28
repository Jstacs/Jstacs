package projects.gemorna;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.FileExistsValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import projects.gemoma.Analyzer;
import projects.gemoma.Analyzer.Transcript;

public class PredictCDSFromGFF implements JstacsTool{

	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new PredictCDSFromGFF());
		cli.run(args);
	}
	
	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<Parameter>();
		
		pars.add(new FileParameter("Genome", "Genome sequence as FastA", "fa,fna.fasta", true));
		
		pars.add(new FileParameter("predicted annotation", "\"GFF or GTF file containing the predicted annotation\"", "gff,gff3,gff.gz,gff3.gz,gtf,gtf.gz", true, new FileExistsValidator()));
		
		try {
			pars.add(new SimpleParameter(DataType.INT,"Minimum protein length","Minimum length of protein in AA",true,70));
		} catch (DatatypeNotValidException e) {
			e.printStackTrace();
		} catch (IllegalValueException e) {
			e.printStackTrace();
		}
		
		return new ToolParameterSet(this.getToolName(), pars);
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		String genome = ((FileParameter)parameters.getParameterForName("Genome")).getFileContents().getFilename();
		
		FileParameter fp = (FileParameter) parameters.getParameterForName("predicted annotation");
		String truth = fp.getValue();
		
		int minProteinLength = (int) parameters.getParameterForName("Minimum protein length").getValue(); 
		
		Genome.init(genome);
		
		HashMap<String,int[]> attributesTruth = new HashMap<String,int[]>();
		attributesTruth.put("ALL",new int[1]);
		ArrayList<String> attTruth = new ArrayList<String>();
		HashMap<String,HashMap<String,Transcript>> res = new Analyzer().readGFF( "exon", truth, protocol, attributesTruth, attTruth );
		
		File out = File.createTempFile("predictions", ".temp", new File("."));
		out.deleteOnExit();
		
		PrintWriter wr = new PrintWriter(out);
		
		for(String key : res.keySet()) {
			HashMap<String,Transcript> curr = res.get(key);
			for(String key2 : curr.keySet()) {
				Transcript t = curr.get(key2);
				SplicingGraph sg = new SplicingGraph(t);
				projects.gemorna.SplicingGraph.Transcript t2 = sg.createTranscript(t);
				try {
					t2.addStrandAndCDS(minProteinLength,false);
				}catch(StringIndexOutOfBoundsException e) {
					t2.setCDSStartEnd(-1, -1);
					protocol.appendWarning("Range of "+t.getID()+" outside of chromosome\n");
				}
				wr.println(t2);
			}
		}
		
		wr.close();
		
		TextResult tr = new TextResult("Transcripts with CDS", "Transcripts in GFF format", new FileParameter.FileRepresentation(out.getAbsolutePath()), "gff3", getToolName(), null, true);

		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	@Override
	public String getToolName() {
		return "Predict CDS from GFF";
	}

	@Override
	public String getToolVersion() {
		return "1.2";
	}

	@Override
	public String getShortName() {
		return "predictCDS";
	}

	@Override
	public String getDescription() {
		return "";
	}

	@Override
	public String getHelpText() {
		return "";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
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
