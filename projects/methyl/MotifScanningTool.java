package projects.methyl;

import java.io.File;
import java.io.FileOutputStream;
import java.io.StringReader;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.Parameter;
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
		
		pars.add(new FileParameter("Model", "Model XML as output by Methyl SlimDimont", "xml", true));
		
		return new ToolParameterSet(getShortName(),pars);
	}


	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		progress.setLast(1.0);
		progress.setCurrent(0.0);
		
		GenDisMixClassifier cl = new GenDisMixClassifier( new StringBuffer( ((FileParameter)parameters.getParameterAt(1)).getFileContents().getContent()));
		
		AlphabetContainer con = cl.getAlphabetContainer();
		
		
		FileRepresentation fr = ((FileParameter) parameters.getParameterAt( 0 )).getFileContents();
		
		DataSet data= new DataSet( con, new SparseStringExtractor( new StringReader(fr.getContent()), '>', "", new SimpleSequenceAnnotationParser()) );
		
		
		
		ThresholdedStrandChIPper model = (ThresholdedStrandChIPper) cl.getDifferentiableSequenceScore(0);
		
		DifferentiableStatisticalModel motif = model.getFunction(0);
		
		
		File out = File.createTempFile("dimontscan", "_dgs.temp");
		out.deleteOnExit(); 
		
		SafeOutputStream sos = SafeOutputStream.getSafeOutputStream(new FileOutputStream(out));
		
		progress.setCurrent(0.2);
		
		double tot = data.getNumberOfElements();
		
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
				temp.add(score);
				if(score > max){
					max = score;
					maxStart = j;
					strand = "-";
				}
			}
			
			double sum = Normalisation.getLogSum(temp.toArray()) - Math.log(seq.getLength()*2);
			
			String outseq = null;
			if(strand.equals("+")){
				outseq = seq.toString(maxStart,maxStart+motif.getLength());
			}else{
				outseq = rc.toString(maxStart,maxStart+motif.getLength());
				
				maxStart = rc.getLength()-motif.getLength()-maxStart;
			}
			
			
			//double fgScore = cl.getScore(seq, 0);
			//double llr = fgScore - cl.getScore(seq, 1);
			
			sos.writeln((i+1)+"\t"+(maxStart+1)+"\t"+strand+"\t"+max+"\t"+sum+"\t"+outseq+"\t"+id/*+"\t"+llr+"\t"+fgScore*/);
			progress.setCurrent(0.2 + (i/tot)*0.8);
		}
		progress.setCurrent(1.0);
		return new ToolResult("Result of "+getToolName(), "", null, new ResultSet( new TextResult("Predictions", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), true, "tsv", getToolName(), null, true)), parameters, getToolName(), new Date());
	}

	@Override
	public String getToolName() {
		return "Sequence Scoring";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "score";
	}

	@Override
	public String getDescription() {
		return "Scan input sequences for motif matches";
	}

	@Override
	public String getHelpText() {
		return "**"+getToolName()+"** scans a set of input sequences (e.g., sequences under ChIP-seq peaks) for a given motif model (provided as XML as output by \"Methyl SlimDimont\" and provides per sequence information of i)"
				+ " the start position and strand of the best motif match, ii) the corresponding maximum score, iii) the log-sum occupancy score, iv) the matching sequence, and v) the ID (FastaA header) of the sequence.\n\n"
				+ "The purpose of this tool mainly is to determine per-sequence scores for classification, for instance, distinguishing bound from unbound sequences.\n\n" +
				"If you experience problems using "+getToolName()+", please contact_ us.\n" + 
				"\n" + 
				".. _contact: mailto:grau@informatik.uni-halle.de";
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
