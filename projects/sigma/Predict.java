package projects.sigma;

import java.io.StringReader;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.DifferentiableHigherOrderHMM;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;

public class Predict implements JstacsTool {

	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new boolean[]{false}, new Predict());
		cli.run(args);
	}
	
	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> pars = new LinkedList<>();
		
		pars.add(new FileParameter("Input sequences", "", "fasta,fas,fa", true));
		pars.add(new FileParameter("Model", "", "xml", true));
		
		try {
			pars.add(new SimpleParameter(DataType.INT, "First number of components", "", true,new NumberValidator<Integer>(0,20),3));
			pars.add(new SimpleParameter(DataType.INT, "First length", "", true,new NumberValidator<Integer>(3,20),10));
			pars.add(new SimpleParameter(DataType.INT, "Second number of components", "", true,new NumberValidator<Integer>(0,20),3));
			pars.add(new SimpleParameter(DataType.INT, "Second length", "", true,new NumberValidator<Integer>(3,20),10));

		} catch (ParameterException e) {
			e.printStackTrace();
		}
		
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
		
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		DataSet data = new DataSet(DNAAlphabetContainer.SINGLETON, new SparseStringExtractor(new StringReader(((FileParameter)parameters.getParameterAt(0)).getFileContents().getContent()), '>', "", new SimpleSequenceAnnotationParser()));

		DifferentiableHigherOrderHMM hmm = new DifferentiableHigherOrderHMM(new StringBuffer(((FileParameter)parameters.getParameterAt(1)).getFileContents().getContent()));
		
		int nComponentsFirst = (int) parameters.getParameterAt(2).getValue();
		int nFirst = (int) parameters.getParameterAt(3).getValue();
		int nComponentsSecond = (int) parameters.getParameterAt(4).getValue();
		int nSecond = (int) parameters.getParameterAt(5).getValue();
		
		
		LinkedList<Result> ress = SigmaHMM.evaluateJoint(this.getShortName(), hmm, nComponentsFirst, nComponentsSecond, nFirst, nSecond, data);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(ress), parameters, getToolName(), new Date(System.currentTimeMillis()));
		
	}
	
	

	@Override
	public String getToolName() {
		return "Prediction";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "predict";
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
		// TODO Auto-generated method stub
		
	}

	@Override
	public String[] getReferences() {
		// TODO Auto-generated method stub
		return null;
	}

}
