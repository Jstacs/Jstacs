package projects.sigma;

import java.io.StringReader;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.annotation.SequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.DataSetResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.DifferentiableHigherOrderHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.HigherOrderHMM;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import projects.sigma.SigmaIter.IntermediateResult;

public class PredictIter implements JstacsTool {

	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new boolean[]{false}, new PredictIter());
		cli.run(args);
	}
	
	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> pars = new LinkedList<>();
		
		pars.add(new FileParameter("Input sequences", "", "fasta,fas,fa", true));
		//pars.add(new FileParameter("Model", "", "xml", true));
		try {
			pars.add(new SimpleParameter(DataType.INT, "Length", "", true,new NumberValidator<Integer>(3,20),10));

		} catch (ParameterException e) {
			e.printStackTrace();
		}
		
		try {
			pars.add(new ParameterSetContainer(new ExpandableParameterSet(new SimpleParameterSet(
					new FileParameter("Model", "", "xml", true)
					),"Models","")));
		} catch (CloneNotSupportedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		
		
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
		
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		SequenceAnnotationParser parser = new SimpleSequenceAnnotationParser();
		
		DataSet data = new DataSet(DNAAlphabetContainer.SINGLETON, new SparseStringExtractor(new StringReader(((FileParameter)parameters.getParameterAt(0)).getFileContents().getContent()), '>', "", parser ));
		
		DataSet orig = data;
		
		int length = (int) parameters.getParameterAt(1).getValue();
		
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterAt(2).getValue();
		
		int nModels = eps.getNumberOfParameters();
		HigherOrderHMM[] models = new HigherOrderHMM[nModels];
		for(int i=0;i<models.length;i++) {
			HigherOrderHMM hmm = new HigherOrderHMM(new StringBuffer(((FileParameter) ((ParameterSetContainer)eps.getParameterAt(i)).getValue().getParameterAt(0) ).getFileContents().getContent()));
			models[i] = hmm;
		}
		
		int off = 1;
		LinkedList<Result> ress = new LinkedList<>();
		
		int[] idx = new int[data.getNumberOfElements()];
		for(int i=0;i<idx.length;i++){
			idx[i] = i;
		}
		
		ResultSet[] rs = new ResultSet[data.getNumberOfElements()];
		
		
		for(int m=0;m<models.length;m++) {
		
			//System.out.println(m);
			
			HigherOrderHMM hmm = models[m];
			
			IntermediateResult ir = SigmaIter.evaluateModel(data, hmm, idx, length, off, rs);
			
			//System.out.println(ir.in.size());
						
			DataSet inSet = new DataSet("",ir.in);
			
			
			DataSetResult dsr = new DataSetResult("Sequences for model "+off, "", inSet);
			dsr.setParser(parser);
			
			ress.add(dsr);
			
			if(ir.out.size() == 0) {
				break;
			}
			
			data = new DataSet("",ir.out);
			idx = ir.idx.toArray();
			
			
			off++;
		}
		
		for(int i=0;i<rs.length;i++) {
			if(rs[i] == null) {
				rs[i] = new ResultSet(new Result[]{
						new NumericalResult("Index", "", i),
						new NumericalResult("start","",-1),
						new CategoricalResult("seq", "", "NA"),
						new NumericalResult("Score", "", Double.NEGATIVE_INFINITY),
						new NumericalResult("component", "", -1),
						new NumericalResult("LLR", "", Double.NEGATIVE_INFINITY),
						new CategoricalResult("annotation", "", (String) orig.getElementAt(i).getAnnotation()[0].getResultAt(0).getValue())
				});
			}
		}
		

		ress.add(new ListResult("Predictions", "", null, rs));

		
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
