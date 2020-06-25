package projects.methyl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.util.Date;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.DataType;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.tools.DataColumnParameter;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.PFMComparator;
import projects.dimont.ThresholdedStrandChIPper;

public class MethylationSensitivity implements JstacsTool {

	public static void main(String[] args) throws NonParsableException, IOException, CloneNotSupportedException, IllegalArgumentException, WrongAlphabetException {
		
		
		
	}

	private static double getScore(DifferentiableStatisticalModel motif, String mg) throws IllegalArgumentException, WrongAlphabetException {
		Sequence seq = Sequence.create(motif.getAlphabetContainer(), mg);
		return motif.getLogScoreFor(seq);
	}

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> parameters = new LinkedList<Parameter>();
		
		parameters.add(new FileParameter("Model", "The XML file containing the model", "xml", true));
		FileParameter fp = new FileParameter("Predictions","The file containing the predictions from the training run","tsv",true); 
		parameters.add(fp);
		
		try {
			parameters.add(new DataColumnParameter(fp.getName(), "Sequence column", "The column of the predictions file containing the sequences in adjusted strand orientation", true,8));
			parameters.add(new SimpleParameter(DataType.BOOLEAN, "Verbose", "Output MpH sensitivity profile for every input sequence", true, false));
		} catch (DatatypeNotValidException | IllegalValueException e) {
			throw new RuntimeException(e);
		}
		
		return new ToolParameterSet(getToolName(), parameters);
		
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		FileRepresentation clFile = ((FileParameter) parameters.getParameterAt(0)).getFileContents();
		FileRepresentation predFile = ((FileParameter) parameters.getParameterAt(1)).getFileContents();
		int seqcol = (int) parameters.getParameterAt(2).getValue() - 1;
		boolean verbose = (boolean) parameters.getParameterAt(3).getValue();
		
		GenDisMixClassifier cl = new GenDisMixClassifier(new StringBuffer(clFile.getContent()));
		
		
		ThresholdedStrandChIPper model = (ThresholdedStrandChIPper) cl.getDifferentiableSequenceScore(0);
		
		DifferentiableStatisticalModel motif = model.getFunction(0);

		BufferedReader read = null;
		if((new File(predFile.getFilename())).exists()) {
			read = new BufferedReader(new FileReader(predFile.getFilename()));
		}else {
			read = new BufferedReader(new StringReader(predFile.getContent()));
		}
		
		
		Pattern pat = Pattern.compile("CG");
		
		double[] cpgs = null;
		double[] mgs = null;
		double[] chs = null;
		double[] mhs = null;
		
		PrintWriter wr = null;
		File tempFile = null;
		if(verbose) {
			tempFile = File.createTempFile("mhprof", ".tsv");
			tempFile.deleteOnExit();
			
			wr = new PrintWriter(tempFile);
			wr.print("sequence");
			for(int i=0;i<motif.getLength();i++) {
				wr.print("\tP"+(i+1));
			}
			wr.println();
			
		}
		
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		
		String str = null;
		int n=0;
		while( (str = read.readLine()) != null ){
		
			if(!str.startsWith("#")) {
			
				String[] parts = str.split("\t");
				String seqstr = parts[ seqcol ];
				seqs.add(Sequence.create(motif.getAlphabetContainer(), seqstr));
				
				if(cpgs == null){
					cpgs = new double[seqstr.length()];
					mgs = new double[seqstr.length()];
					chs = new double[seqstr.length()];
					mhs = new double[seqstr.length()];
				}
				
				double[] mhl = new double[seqstr.length()];
				
				seqstr = seqstr.replaceAll("M", "C").replaceAll("H", "G");
				
				double base = getScore(motif,seqstr);
				
				Matcher m = pat.matcher(seqstr);
				
				while(m.find()){
					int start = m.start();
					int end = m.end();
					String mg = seqstr.substring(0, start)+"MG"+seqstr.substring(end);
					String ch = seqstr.substring(0, start)+"CH"+seqstr.substring(end);
					String mh = seqstr.substring(0, start)+"MH"+seqstr.substring(end);
					
					double mgVal = getScore(motif,mg);
					double chVal = getScore(motif,ch);
					double mhVal = getScore(motif,mh);
					
					cpgs[start]++;
					mgs[start] += mgVal-base;
					chs[start] += chVal-base;
					mhs[start] += mhVal-base;
					mhl[start] = mhVal-base;
					
				}
				n++;
				
				if(verbose) {
					wr.print(seqstr);
					for(int i=0;i<mhl.length;i++) {
						wr.print("\t");
						wr.print(mhl[i]);
					}
					wr.println();
				}
			}
		}
		read.close();
		
		LinkedList<Result> resll = new LinkedList<>();
		
		if(verbose) {
			wr.close();
			
			TextResult fr = new TextResult("MH profile", "", new FileParameter.FileRepresentation(tempFile.getAbsolutePath()),true,"tsv",this.getToolName(),null,true);
			
			resll.add(fr);
			
		}
		
		StringBuffer sb = new StringBuffer();
		
		sb.append("Position\tCpG\tMG\tCH\tMH\n");
		for(int i=0;i<cpgs.length-1;i++){
			cpgs[i] /= n;
			mhs[i] /= n;
			sb.append((i+1)+"\t"+(cpgs[i])+"\t"+(mgs[i]/n)+"\t"+(chs[i]/n)+"\t"+(mhs[i])+"\n");
		}
		
		TextResult tr = new TextResult("Average methylation sensitivity", "Methylation sensitivity averaged over all input sequences per content", new FileParameter.FileRepresentation("", sb.toString()), "tsv", this.getToolName(), null, true);
		resll.add(tr);
		
		DataSet data = new DataSet("", seqs);
		
		double[][] pwm = PFMComparator.getPWM(data, 0, data.getNumberOfElements());
		
		
		PlotGeneratorResult pgr = new PlotGeneratorResult("Methylation sensitivity", "Plot of a sequence logo and corresponding MpH sensitivity profile", new MLogoPlotter.MLogoPlotGenerator(pwm, cpgs, mhs, 1000), true);
		
		resll.add(pgr);
		
		
		return new ToolResult("Result of "+getToolName(), "", null, new ResultSet(resll), parameters, getToolName(), new Date());
		
		
	}

	@Override
	public String getToolName() {
		return "Methylation Sensitivity";
	}

	@Override
	public String getToolVersion() {
		return "0.2";
	}

	@Override
	public String getShortName() {
		return "msens";
	}

	@Override
	public String getDescription() {
		return "";
	}

	@Override
	public String getHelpText() {
		return "**"+getToolName()+"** determines average methylation sensitivity profiles for CpG dinucleotides converted to MpG, CpH, and MpH. As input, it needs a model XML as generated by \"Methyl SlimDimont\", "
				+ "and a prediction file as output from the corresponding training run.\n\n"
				+ "Optionally, Methylation Sensitivity also generates per-sequence methylation sensitivity profiles for the MpH context.\n\n" +
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
