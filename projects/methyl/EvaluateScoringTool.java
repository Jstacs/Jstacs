package projects.methyl;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.classifiers.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.classifiers.performanceMeasures.PRCurve;
import de.jstacs.classifiers.performanceMeasures.ROCCurve;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.DoubleList;

public class EvaluateScoringTool implements JstacsTool {


	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<Parameter>();
		FileParameter fp =  new FileParameter("Positives", "Output of \"Sequence Scoring\" for the positive test sequences.", "tsv,tabular", true);
		pars.add(fp);
		FileParameter fp2 =  new FileParameter("Negatives", "Output of \"Sequence Scoring\" for the negative test sequences.", "tsv,tabular", true);
		pars.add(fp2);
		try {
			SimpleParameter sp = new SimpleParameter(DataType.BOOLEAN, "Curves", "Also compute and draw ROC and PR curves", true, false);
			pars.add(sp);
			
			SimpleParameter sp2 = new SimpleParameter(DataType.BOOLEAN, "Use sum-occupancy", "Use log-sum occupancy score instead of maximum score", true, false);
			pars.add(sp2);
		} catch (DatatypeNotValidException | IllegalValueException e) {
			throw new RuntimeException(e);
		}
		
		return new ToolParameterSet(getToolName(), pars);
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		FileParameter.FileRepresentation posFile = ((FileParameter)parameters.getParameterAt(0)).getFileContents();
		FileParameter.FileRepresentation negFile = ((FileParameter)parameters.getParameterAt(1)).getFileContents();
		
		boolean curves = (boolean) parameters.getParameterAt(2).getValue();
		
		boolean sumOcc = (boolean) parameters.getParameterAt(3).getValue();
		
		int col = 3;
		if(sumOcc) {
			col = 4;
		}
		
		double[] pos = parse(posFile,col);
		double[] neg = parse(negFile,col);
		
		
		ROCCurve roc = new ROCCurve();
		ResultSet rocRes = roc.compute(pos, neg);
		
		PRCurve pr = new PRCurve();
		ResultSet prRes = pr.compute(pos,neg);
				
		LinkedList<Result> ress = new LinkedList<Result>();
		
		ListResult lr = new ListResult("Areas under curve", "", null, 
				new ResultSet(new Result[] {new CategoricalResult("Curve", "", "ROC"),
						new NumericalResult("Value", "", (double)rocRes.getResultAt(0).getValue())}),
				new ResultSet(new Result[] {new CategoricalResult("Curve", "", "PR"),
						new NumericalResult("Value", "", (double)prRes.getResultAt(1).getValue())}));
		ress.add(lr);
		
		if(curves) {
		
			DoubleTableResult dtrROC = (DoubleTableResult) rocRes.getResultAt(1);
			ListResult lrROC = getListResult(dtrROC,"FPR","Sn", "ROC points");
			ress.add(lrROC);
			
			PlotGeneratorResult pgrROC = getPlot(dtrROC,"FPR", "Sn", "ROC curve");
			ress.add(pgrROC);
			
			DoubleTableResult dtrPR = (DoubleTableResult) prRes.getResultAt(2);
			ListResult lrPR = getListResult(dtrPR,"Recall","Precision", "PR points");
			ress.add(lrPR);
			
			PlotGeneratorResult pgrPR = getPlot(dtrPR,"Recall","Precision", "PR curve");
			ress.add(pgrPR);
			
		
		}
		return new ToolResult("Result of "+getToolName(), "", null, new ResultSet(ress), parameters, getToolName(), new Date());
	}

	private PlotGeneratorResult getPlot(DoubleTableResult dtr, String xlab, String ylab, String name) {
		CurvePlotter plotter = new CurvePlotter(dtr, xlab, ylab);
		
		return new PlotGeneratorResult(name, "", plotter, true);
	}

	private ListResult getListResult(DoubleTableResult dtr, String x, String y, String name) {
		
		LinkedList<ResultSet> li = new LinkedList<ResultSet>();
		
		
		for(int i=0;i<dtr.getNumberOfLines();i++) {
			double[] line = dtr.getLine(i);
			li.add(new ResultSet(new Result[] {new NumericalResult(x, "", line[0]),new NumericalResult(y, "", line[1])}));
		}
		
		ListResult lr = new ListResult(name, "", null, li);
		
		return lr;
	}

	private double[] parse(FileRepresentation file, int col) throws NumberFormatException, IOException {
		
		BufferedReader br = new BufferedReader(new StringReader(file.getContent()));
		
		String line = null;
		
		DoubleList vals = new DoubleList();
		
		while( (line = br.readLine()) != null ) {
			String[] parts = line.split("\t");
			double val = Double.parseDouble(parts[col]);
			vals.add(val);
		}
		vals.sort();
		return vals.toArray();
		
	}

	@Override
	public String getToolName() {
		return "Evaluate Scoring";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "eval";
	}

	@Override
	public String getDescription() {
		return "Evaluate a scoring of sequences";
	}

	@Override
	public String getHelpText() {
		return "**"+getToolName()+"** computes the area under the ROC curve and under the precision recall curve "
				+ "based on the scoring of a positive and a negative set of sequences. Optionally, also the curves may be drawn.\n\n" +
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
