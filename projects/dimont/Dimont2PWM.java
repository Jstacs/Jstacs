package projects.dimont;

import java.io.File;
import java.io.PrintWriter;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.MarkovModelDiffSM;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.SeqLogoPlotter.SeqLogoPlotGenerator;

public class Dimont2PWM implements JstacsTool {

	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new Dimont2PWM());
		
		cli.run(args);
	}

	@Override
	public ToolParameterSet getToolParameters() {
		return new ToolParameterSet(getToolName(), new FileParameter("Dimont XML", "", "xml", true));
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		String name = ((FileParameter)parameters.getParameterAt(0)).getFileContents().getFilename();
		GenDisMixClassifier cl = new GenDisMixClassifier( new StringBuffer(((FileParameter)parameters.getParameterAt(0)).getFileContents().getContent()) );

				
		double[][] modelPwm = ((MarkovModelDiffSM)((ThresholdedStrandChIPper)((GenDisMixClassifier)cl).getDifferentiableSequenceScore(0)).getMotifModel()).getPWM();

		StringBuffer sb = new StringBuffer();
		sb.append(">"+name+"\n");
		for(int i=0;i<modelPwm.length;i++) {
			for(int j=0;j<modelPwm[i].length;j++) {
				if(j>0) {
					sb.append("\t");
				}
				sb.append(modelPwm[i][j]);
			}
			sb.append("\n");
		}
		
		File outfile = File.createTempFile("pwm", ".temp.gz", new File("."));
		outfile.deleteOnExit();
		
		PrintWriter wr = new PrintWriter(outfile);
		wr.print(sb);
		wr.close();
		
		
		LinkedList<Result> results = new LinkedList<>();
		
		
		TextResult tr = new TextResult("PWM", "", new FileParameter.FileRepresentation(outfile.getAbsolutePath()), "pwm" , getToolName(), null, true);
		
		
		results.add(tr);
		
		results.add(new PlotGeneratorResult( "Sequence logo", "Sequence logo of motif ", 
				new SeqLogoPlotGenerator(modelPwm, 1000), true));
		
		results.add(new PlotGeneratorResult( "Sequence logo (rc)", "Sequence logo of the reverse complement of motif ", 
				new SeqLogoPlotGenerator(PFMComparator.getReverseComplement( DNAAlphabet.SINGLETON, modelPwm ), 1000), true));
		
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(results), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	@Override
	public String getToolName() {
		return "Dimont2PWM";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "dimont2pwm";
	}

	@Override
	public String getDescription() {
		return "Extract the internal Dimont model into a PWM";
	}

	@Override
	public String getHelpText() {
		return "Only works for PWM models within Dimont";
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
