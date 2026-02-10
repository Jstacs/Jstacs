package projects.sigma;

import java.io.FileReader;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.HMMFactory;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.Pair;

public class ScanHMMer implements JstacsTool {

	public static void main(String[] args) throws Exception {
		
		CLI cli = new CLI(new boolean[]{true,false,false,false}, new SigmaHMM(), new Predict(), new GenomicScan(), new ScanHMMer());
		cli.run(args);

	}

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<>();
		
		pars.add(new FileParameter("Input sequences", "", "fasta,fas,fa", true));
		pars.add(new FileParameter("Model", "", "hmm", true));
		
		try {
			pars.add(new SimpleParameter(DataType.INT, "Offset", "", true,42));
			pars.add(new SimpleParameter(DataType.INT, "Model length", "", true,48));

		} catch (ParameterException e) {
			e.printStackTrace();
		}
		
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		String seqFile = (String) parameters.getParameterAt(0).getValue();
		String modelFile = (String) parameters.getParameterAt(1).getValue();
		
		int off = (int) parameters.getParameterAt(2).getValue();
		int len = (int) parameters.getParameterAt(3).getValue();
		
		DataSet ds = new DNADataSet(seqFile, '>', new SimpleSequenceAnnotationParser());
		
		
		StringBuffer cons = new StringBuffer();
		
		Pair<AbstractHMM, HomogeneousMMDiffSM> hmmer = HMMFactory.parseProfileHMMFromHMMer( new FileReader(modelFile),cons , null, null );
		
		AbstractHMM hmm = hmmer.getFirstElement();
		HomogeneousDiffSM bg = hmmer.getSecondElement();
		//System.out.println(bg);
		
		LinkedList<ResultSet> coll = new LinkedList<>();
		
		
		for(int i=0;i<ds.getNumberOfElements();i++) {
			Sequence seq = ds.getElementAt(i);
			//Sequence sub = seq.getSubSequence(off, len);
			//System.out.println(sub);
			double hmmprob = hmm.getLogProbFor(seq,off,off+len-1);
			//double hmmprob2 = hmm.getLogProbFor(sub);
			//System.out.println(hmmprob+" "+hmmprob2);
			double bgprob = bg.getLogProbFor(seq,off,off+len-1);
			//System.out.println(seq.getAnnotation()[0].getResultAt(0).getValue()+"\t"+hmmprob+"\t"+bgprob+"\t"+(hmmprob-bgprob));
			ResultSet rs = new ResultSet(new Result[]{
					new NumericalResult("Index", "", i+1),
					new NumericalResult("Score", "", hmmprob),
					new NumericalResult("Bg", "", bgprob),
					new NumericalResult("LLR", "", hmmprob-bgprob),
					new CategoricalResult("annotation", "", (String) seq.getAnnotation()[0].getResultAt(0).getValue())
			});
			coll.add(rs);
		}
		
		ListResult lr = new ListResult("Predictions", "", null, coll);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(lr), parameters, getToolName(), new Date(System.currentTimeMillis()));
	}

	@Override
	public String getToolName() {
		return  "Scan HMMer";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "hmmer";
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
