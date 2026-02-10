package projects.sigma;

import java.io.StringReader;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.LinkedList;
import java.util.Locale;

import de.jstacs.DataType;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.DataSetResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.HigherOrderHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.SilentEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete.DiscreteEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.ViterbiParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ProtocolOutputStream;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter;

public class SigmaIter implements JstacsTool {

	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new boolean[]{true,}, new SigmaIter());
		cli.run(args);

	}

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<>();
		
		pars.add(new FileParameter("Input sequences", "", "fasta,fas,fa", true));
		
		try {
			pars.add(new SimpleParameter(DataType.INT, "Start", "", true,50));
			pars.add(new SimpleParameter(DataType.INT, "End", "", true,70));
			pars.add(new SimpleParameter(DataType.INT, "Length", "", true,new NumberValidator<Integer>(3,20),10));
			pars.add(new SimpleParameter(DataType.INT, "Number of Starts", "", true,500));
		} catch (ParameterException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		SequenceAnnotationParser parser = new SimpleSequenceAnnotationParser();
		
		DataSet data = new DataSet(DNAAlphabetContainer.SINGLETON, new SparseStringExtractor(new StringReader(((FileParameter)parameters.getParameterAt(0)).getFileContents().getContent()), '>', "", parser ));
		
		DataSet orig = data;
		
		int start = (int) parameters.getParameterAt(1).getValue();
		int end = (int) parameters.getParameterAt(2).getValue();
		int length = (int) parameters.getParameterAt(3).getValue();
		int nStarts = (int) parameters.getParameterAt(4).getValue();
		
		int off = 1;
		LinkedList<Result> ress = new LinkedList<>();
		
		int[] idx = new int[data.getNumberOfElements()];
		for(int i=0;i<idx.length;i++){
			idx[i] = i;
		}
		
		ResultSet[] rs = new ResultSet[data.getNumberOfElements()];
		
		
		while(data.getNumberOfElements() > 10) {
		
			HigherOrderHMM hmm = trainModel(data,start,end,length,16.0,nStarts,protocol,threads);
		
			ress.add(new StorableResult("HMM"+(off), "", hmm));
			
			IntermediateResult ir = evaluateModel(data, hmm, idx, length, off, rs);
			
			
			if(ir.in.size() == 0) {
				break;
			}
			
			addSeqLogos(ir.models.getFirstElement(), ress, off);
			
			
			
			
			
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
	
	public static class IntermediateResult{
		
		IntList idx; 
		LinkedList<Sequence> in; 
		LinkedList<Sequence> out;
		Pair<double[][],double[]> models;
		
		public IntermediateResult(IntList idx, LinkedList<Sequence> in, LinkedList<Sequence> out,
				Pair<double[][], double[]> models) {
			super();
			this.idx = idx;
			this.in = in;
			this.out = out;
			this.models = models;
		}
		
		
		
	}
	
	static IntermediateResult evaluateModel(DataSet data, HigherOrderHMM hmm, int[] idx, int length, int off, ResultSet[] rs) throws Exception {
		Pair<double[][],double[]> models = getModels(hmm,length);
		
		LinkedList<Sequence> in = new LinkedList<>();
		IntList idx2 = new IntList();
		
		LinkedList<Sequence> out = new LinkedList<>();
	
		
		for(int i=0;i<data.getNumberOfElements();i++) {
			
			Sequence seq = data.getElementAt(i);
			
			Pair<IntList,Double> pair2 = hmm.getViterbiPathFor(seq);
			String[] path = hmm.decodePath(pair2.getFirstElement());
			//System.out.println(Arrays.toString(path));
			int seqStart = -1;
			int seqEnd = -1;
			for(int j=0;j<path.length;j++) {
				if(path[j].startsWith("M0")) {
					seqStart = j;
				}
				if(path[j].startsWith("E")) {
					seqEnd = j;
					break;
				}
			}
			
			if(seqStart>-1) {
			
				Sequence sub = seq.getSubSequence(seqStart, seqEnd-seqStart);
				
				double llr = getLLR(sub, models.getFirstElement(), models.getSecondElement());
				//System.out.println(i+" "+idx[i]+" "+seqStart+" "+llr);
				if(llr > Math.log(2)) {
					in.add(seq);
				}else {
					out.add(seq);
					idx2.add(idx[i]);
				}
				
				rs[idx[i]] = new ResultSet(new Result[]{
						new NumericalResult("Index", "", idx[i]),
						new NumericalResult("start","",seqStart),
						new CategoricalResult("seq", "", sub.toString()),
						new NumericalResult("Score", "", pair2.getSecondElement()),
						new NumericalResult("component", "", off),
						new NumericalResult("LLR", "", llr),
						new CategoricalResult("annotation", "", (String) seq.getAnnotation()[0].getResultAt(0).getValue())
				});
			
			}else {
				System.out.println("else");
				out.add(seq);
				idx2.add(idx[i]);
				
				rs[idx[i]] = new ResultSet(new Result[]{
						new NumericalResult("Index", "", idx[i]),
						new NumericalResult("start","",seqStart),
						new CategoricalResult("seq", "", "NA"),
						new NumericalResult("Score", "", Double.NEGATIVE_INFINITY),
						new NumericalResult("component", "", -1),
						new NumericalResult("LLR", "", Double.NEGATIVE_INFINITY),
						new CategoricalResult("annotation", "", (String) seq.getAnnotation()[0].getResultAt(0).getValue())
				});
			}
			
		}
		
		return new IntermediateResult(idx2,in,out,models);
	}
	
	private static double getLLR(Sequence sub, double[][] pwm, double[] bgMod) {
		double score = 0;
		for(int i=0;i<pwm.length;i++){
			score += pwm[i][sub.discreteVal(i)] - bgMod[sub.discreteVal(i)];
		}
		return score;
	}
	
	private static void addSeqLogos(double[][] logPWM, LinkedList<Result> ress, int off) throws CloneNotSupportedException {
		
		logPWM = ArrayHandler.clone(logPWM);
		
		for(int j=0;j<logPWM.length;j++){//component
			Normalisation.logSumNormalisation(logPWM[j]);
	
		}
		
		ress.add(new PlotGeneratorResult("SeqLogo_"+off, "", new SeqLogoPlotter.SeqLogoPlotGenerator(logPWM, 100), true));
	}
	
	private static Pair<double[][],double[]> getModels(HigherOrderHMM hmm, int length) throws CloneNotSupportedException{
		Emission[] emission = (Emission[]) hmm.getEmissions();
		
		double[][] logPWM = new double[length][];
		
		for(int i=0;i<length;i++) {
			DifferentiableEmission temp = ((DiscreteEmission)emission[i+1]).clone();
			temp.setParameterOffset(0);
			double[] pars = new double[ temp.getNumberOfParameters() ];
			temp.fillCurrentParameter(pars);
			logPWM[i] = pars.clone();
		}
		DifferentiableEmission temp = ((DiscreteEmission)emission[0]).clone();
		temp.setParameterOffset(0);
		double[] hom = new double[ temp.getNumberOfParameters() ];
		temp.fillCurrentParameter(hom);
		
		return new Pair<double[][],double[]>(logPWM,hom);
	}
	
	private HigherOrderHMM trainModel(DataSet data, int start, int end, int len, double ess, int nStarts, Protocol protocol, int threads) throws Exception {
		
		int seqLen = data.getElementLength();
		
		DiscreteEmission insert = new DiscreteEmission(DNAAlphabetContainer.SINGLETON, ess*(seqLen-len) + ess*len/2.0);

		DifferentiableEmission[] emission = new DifferentiableEmission[1 + len + 1];
		emission[0] = insert;
		for(int i=0;i<len;i++){
			emission[i+1] = new DiscreteEmission(DNAAlphabetContainer.SINGLETON, ess/2.0);
		}
		emission[emission.length-1] = new SilentEmission();
		
		
		int[] emissionIdx = new int[start-1 + 1 + 2*len + 1 + seqLen-end-1 + 1];
		String[] name = new String[emissionIdx.length];

		ArrayList<TransitionElement> tes = new ArrayList<>();
		int k=0;
		
		tes.add(new TransitionElement(null, new int[]{0}, new double[]{ess}));
		
		for(int i=0;i<start-1;i++,k++) {
			emissionIdx[k] = 0;
			name[k] = "O"+k;
			tes.add(new TransitionElement(new int[]{k}, new int[]{k+1}, new double[]{ess}));

		}
		
		double fac = ((end-start+1)-len)/2.0;
		
		emissionIdx[k] = 0;
		name[k] = "B";
		tes.add(new TransitionElement(new int[]{k}, new int[]{k,k+1,k+1+len}, new double[]{ess*fac,ess/2.0,ess/2.0}));
		k++;
		
		for(int i=0;i<len;i++,k++) {
			emissionIdx[k] = i+1;
			name[k] = "M"+i;
			if(i<len-1) {
				tes.add(new TransitionElement(new int[]{k}, new int[]{k+1}, new double[]{ess/2.0}));
			}else {
				tes.add(new TransitionElement(new int[]{k}, new int[]{k+1+len}, new double[]{ess/2.0}));
			}
		}
		
		for(int i=0;i<len;i++,k++) {
			emissionIdx[k] = 0;
			name[k] = "A"+i;
			tes.add(new TransitionElement(new int[]{k}, new int[]{k+1}, new double[]{ess/2.0}));
		}
		
		emissionIdx[k] = 0;
		name[k] = "E";
		tes.add(new TransitionElement(new int[]{k}, new int[]{k,k+1}, new double[]{ess*fac,ess}));
		k++;
		
		for(int i=end+1;i<seqLen;i++,k++) {
			emissionIdx[k] = 0;
			name[k] = "O"+k;

			tes.add(new TransitionElement(new int[]{k}, new int[]{k+1}, new double[]{ess}));
		}
		
		
		name[k] = "F";
		emissionIdx[k] = emission.length-1;
		//tes.add(new TransitionElement(new int[]{k}, new int[]{k,k+1}, new double[]{ess*0.9,ess*0.1}));
		
		boolean[] forward = new boolean[name.length];
		Arrays.fill(forward, true);
		
		ViterbiParameterSet trainingParameterSet = new ViterbiParameterSet(nStarts, new SmallDifferenceOfFunctionEvaluationsCondition(1E-6), threads);
		
		HigherOrderHMM hmm = new HigherOrderHMM(trainingParameterSet, name, emissionIdx, forward, emission, tes.toArray(new TransitionElement[0]));
		
		hmm.setOutputStream(new ProtocolOutputStream(protocol, true));
		
		//System.out.println(hmm.getGraphvizRepresentation(new DecimalFormat("##.###",new DecimalFormatSymbols(Locale.US))));
		
		hmm.train(data);
		
		return hmm;
		
	}

	@Override
	public String getToolName() {
		return "iter";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "iter";
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
