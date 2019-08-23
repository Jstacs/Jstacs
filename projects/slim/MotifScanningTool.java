package projects.slim;

import java.io.File;
import java.io.FileOutputStream;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.io.FileManager;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.SafeOutputStream;
import projects.dimont.ThresholdedStrandChIPper;

public class MotifScanningTool implements JstacsTool {

	public static void main( String[] args ) throws Exception {
		
		CLI cli = new CLI(new MotifScanningTool());
		
		cli.run(args);

	}
	
	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> pars = new LinkedList<>();
		pars.add(new FileParameter("Input sequences", "Input sequences in FastA format", "fasta,fa,fas", true));
		
		pars.add(new FileParameter("Model", "Model XML", "xml", true));
		
		return new ToolParameterSet(getShortName(), pars.toArray(new Parameter[0]));
	}


	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		GenDisMixClassifier cl = new GenDisMixClassifier(FileManager.readFile((String) parameters.getParameterAt(1).getValue()));
		
		AlphabetContainer con = cl.getAlphabetContainer();
		
		
		String filename = (String)(parameters.getParameterAt( 0 ).getValue());
		DataSet data= new DataSet( con, new SparseStringExtractor( filename, '>', new SimpleSequenceAnnotationParser()) );
		
		
		
		ThresholdedStrandChIPper model = (ThresholdedStrandChIPper) cl.getDifferentiableSequenceScore(0);
		
		DifferentiableStatisticalModel motif = model.getFunction(0);
		
		
		File out = File.createTempFile("dimontscan", "_dgs.temp", new File("."));
		out.deleteOnExit(); 
		
		SafeOutputStream sos = SafeOutputStream.getSafeOutputStream(new FileOutputStream(out));
		
		
		
		for(int i=0;i<data.getNumberOfElements();i++){
			Sequence seq = data.getElementAt(i);
			String id = (String) seq.getSequenceAnnotationByType("unparsed comment line", 0).getResultAt(0).getValue();
			
			DoubleList temp = new DoubleList();
			
			double max = Double.NEGATIVE_INFINITY;
			int maxStart = 0;
			String strand = "+";
			
			for(int j=0;j<seq.getLength()-motif.getLength()+1;j++){
				double score = motif.getLogScoreFor(seq, j);
				temp.add(score);
				if(score > max){
					max = score;
					maxStart = j;
				}
			}
			
			Sequence rc = seq.reverseComplement();
			for(int j=0;j<rc.getLength()-motif.getLength()+1;j++){
				double score = motif.getLogScoreFor(rc, j);
				if(score > max){
					max = score;
					maxStart = j;
					strand = "-";
				}
				temp.add(score);
			}
			
			double sum = Normalisation.getLogSum(temp.toArray()) - Math.log(seq.getLength()*2);
			
			String outseq = null;
			if(strand.equals("+")){
				outseq = seq.toString(maxStart,maxStart+motif.getLength());
			}else{
				outseq = rc.toString(maxStart,maxStart+motif.getLength());
				
				maxStart = rc.getLength()-motif.getLength()-maxStart;
			}
			
			sos.writeln((i+1)+"\t"+(maxStart+1)+"\t"+strand+"\t"+max+"\t"+sum+"\t"+outseq+"\t"+id);
			
		}
		
		return new ToolResult("predictions", "", null, new ResultSet( new TextResult("predictions", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "txt", getToolName(), null, true)), parameters, getToolName(), new Date());
	}

	@Override
	public String getToolName() {
		return "MotifScan";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "scan";
	}

	@Override
	public String getDescription() {
		return "Scan input sequences for motif matches";
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
