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
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.results.TextResult;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.DifferentiableHigherOrderHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.HigherOrderHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.SilentEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete.DiscreteEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.BaumWelchParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.MaxHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.NumericalHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.ViterbiParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.NumericalHMMTrainingParameterSet.TrainingType;
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
import de.jstacs.utils.graphics.PDFAdaptor;

public class SigmaHMM implements JstacsTool {

	public enum Training{
		VITERBI,
		BAUMWELCH,
		GRADIENT
	}
	
	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new boolean[]{true,false,false,false}, new SigmaHMM(), new Predict(), new GenomicScan(), new ScanHMMer());
		cli.run(args);
		
	}

	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> pars = new LinkedList<>();
		
		pars.add(new FileParameter("Input sequences", "", "fasta,fas,fa", true));
		
		try {
			//pars.add(new SimpleParameter(DataType.INT, "Number of components", "", true,new NumberValidator<Integer>(1,20),3));
			pars.add(new SelectionParameter(DataType.PARAMETERSET, new String[]{"joint","iterative"}, new ParameterSet[]{
					new SimpleParameterSet(new SimpleParameter(DataType.INT, "Number of components", "", true,new NumberValidator<Integer>(0,20),3)),
					new SimpleParameterSet(new SimpleParameter(DataType.INT, "Minimum number of sequences", "", true,new NumberValidator<Integer>(1,Integer.MAX_VALUE),10),
							new SimpleParameter(DataType.DOUBLE, "Threshold", "", true,0.0))
					
			}, "Strategy", "Training strategy", true));
			pars.add(new SimpleParameter(DataType.INT, "First length", "", true,new NumberValidator<Integer>(3,20),10));
			pars.add(new SimpleParameter(DataType.INT, "Second length", "", true,new NumberValidator<Integer>(3,20),10));
			
			pars.add(new SelectionParameter(DataType.PARAMETERSET, new String[]{"minimum","range"}, new ParameterSet[]{
					new SimpleParameterSet(
							new SimpleParameter(DataType.INT, "Minimum distance", "", true,new NumberValidator<Integer>(1,Integer.MAX_VALUE),10)
					),
					new SimpleParameterSet(
							new SimpleParameter(DataType.INT, "Minimum distance", "", true,new NumberValidator<Integer>(1,Integer.MAX_VALUE),10),
							new SimpleParameter(DataType.INT, "Maximum distance", "", true,new NumberValidator<Integer>(1,Integer.MAX_VALUE),10)
					)
			}, "Distance", "Distance between boxes", true));
			
			pars.add(new SimpleParameter(DataType.INT, "Offset", "", true,new NumberValidator<Integer>(1,Integer.MAX_VALUE),40));
			pars.add(new SimpleParameter(DataType.DOUBLE, "ESS", "", true,new NumberValidator<Double>(0.0,Double.MAX_VALUE),4.0));
			pars.add(new EnumParameter(Training.class, "", true));
			pars.add(new SimpleParameter(DataType.INT, "Starts", "", true,new NumberValidator<Integer>(1,Integer.MAX_VALUE),40));
			//pars.add(new SimpleParameter(DataType.BOOLEAN, "Allow switch", "Allow switches between first and second component", true, false));
			pars.add(new SelectionParameter(DataType.PARAMETERSET, new String[]{"true","false"}, new ParameterSet[]{
					new SimpleParameterSet(new SimpleParameter(DataType.INT, "Number of components 2", "", true,new NumberValidator<Integer>(0,20),3)),
					new SimpleParameterSet()
					
			}, "Allow switch", "Allow switches between first and second component", true));
			
		} catch (ParameterException e) {
			e.printStackTrace();
		}
		
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
		
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		
		DataSet data = new DataSet(DNAAlphabetContainer.SINGLETON, new SparseStringExtractor(new StringReader(((FileParameter)parameters.getParameterAt(0)).getFileContents().getContent()), '>', "", new SimpleSequenceAnnotationParser()));
		
		int nStarts = (int) parameters.getParameterAt(8).getValue();
		
		
		int nFirst = (int) parameters.getParameterAt(2).getValue();
		int nSecond = (int) parameters.getParameterAt(3).getValue();
		int seqlen = data.getElementLength();
		
		int minDist = 0;
		int maxDist = Integer.MAX_VALUE;
		if( ((SelectionParameter)parameters.getParameterAt(4)).getSelected()==0 ){
			minDist = (int) ((ParameterSet)parameters.getParameterAt(4).getValue()).getParameterAt(0).getValue();
		}else{
			minDist = (int) ((ParameterSet)parameters.getParameterAt(4).getValue()).getParameterAt(0).getValue();
			maxDist = (int) ((ParameterSet)parameters.getParameterAt(4).getValue()).getParameterAt(1).getValue();
		}

		int minOffset = (int) parameters.getParameterAt(5).getValue();
		
		
		double ess = (double) parameters.getParameterAt(6).getValue();
		
		Enum<Training> val = ((EnumParameter)parameters.getParameterAt(7)).getValue();
		
		
		boolean allowSwitch = ((SelectionParameter)parameters.getParameterAt(9)).getSelected()==0;
		
		
		
		LinkedList<Result> ress = null;
		
		if(((SelectionParameter)parameters.getParameterAt(1)).getSelected()==0){
			int nComponents = (int) ((ParameterSet) parameters.getParameterAt(1).getValue()).getParameterAt(0).getValue();
			int nComponents2 = nComponents;
			if(allowSwitch){
				nComponents2 = (int) ((ParameterSet)((SelectionParameter)parameters.getParameterAt(9)).getValue()).getParameterAt(0).getValue();
			}
			
			HigherOrderHMM hmm = trainModel(data, nComponents, nComponents2, nStarts, seqlen, nFirst, nSecond, minOffset, minDist, maxDist, ess, val, allowSwitch, threads, protocol);
		
			ress = evaluateJoint(this.getShortName(),hmm, nComponents, nComponents2, nFirst, nSecond, data);
		}else{
			
			int minNumber = (int) ((ParameterSet) parameters.getParameterAt(1).getValue()).getParameterAt(0).getValue();
			double threshold = (double) ((ParameterSet) parameters.getParameterAt(1).getValue()).getParameterAt(1).getValue();
			
			ress = trainIteratively(data, minNumber, nStarts, seqlen, nFirst, nSecond, minOffset, minDist, maxDist, ess, val, allowSwitch, threads, protocol, threshold);
			
		}
		
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(ress), parameters, getToolName(), new Date(System.currentTimeMillis()));
	}
	
	
	private LinkedList<Result> trainIteratively(DataSet data, int minNumber, int nStarts, int seqlen, int nFirst, int nSecond,
			int minOffset, int minDist, int maxDist, double ess, Enum<Training> train, boolean allowSwitch, int threads, Protocol protocol, double threshold) throws Exception{
		
		LinkedList<Result> ress = new LinkedList<>();
		
		DataSet current = data;
		int[] idx = new int[current.getNumberOfElements()];
		for(int i=0;i<idx.length;i++){
			idx[i] = i;
		}
		ResultSet[] rs = new ResultSet[current.getNumberOfElements()];
		
		int comp = 0;
		
		while(current.getNumberOfElements()>minNumber){
			
			HigherOrderHMM hmm = trainModel(current, 1, 1, nStarts, seqlen, nFirst, nSecond, minOffset, minDist, maxDist, ess, train, allowSwitch, threads, protocol);
			
			ress.add(new StorableResult("HMM"+(comp+1), "", hmm));
			
			DifferentiableEmission[] emission = (DifferentiableEmission[]) hmm.getEmissions();
			
			Pair<double[][][][],double[]> mods = getModels(1,1, nFirst, nSecond, emission);
			double[][][][] logPWMs = mods.getFirstElement();
			double[] bgMod = mods.getSecondElement();
			
			addSeqLogos(logPWMs, ress, comp);
			
			LinkedList<int[]> use = new LinkedList<>();
			IntList idx2 = new IntList();
			
			for(int i=0;i<current.getNumberOfElements();i++){
				
				Sequence seq = current.getElementAt(i);
				
				Pair<IntList,Double> pair2 = hmm.getViterbiPathFor(seq);
				String[] path = hmm.decodePath(pair2.getFirstElement());
				
				Pair<int[], Sequence[]> pair = getViterbiSeqs(path, seq);
				
				double llr = getLLR(pair.getSecondElement()[0], pair.getSecondElement()[1], logPWMs[0][0],logPWMs[0][1], bgMod);
				
				int startF = pair.getFirstElement()[2];
				int startS = pair.getFirstElement()[3];
				
				Sequence sub1 = pair.getSecondElement()[0];
				Sequence sub2 = pair.getSecondElement()[1];
				
				
				int myComp = 0;
				
				if(llr > threshold){
					
					myComp = comp+1;
					
				}else{
					use.add(new int[]{i,i+1});
					idx2.add(idx[i]);
				}
				
				rs[idx[i]] = new ResultSet(new Result[]{
						new NumericalResult("Index", "", i+1),
						new NumericalResult("start1","",startF),
						new CategoricalResult("seq1", "", sub1.toString()),
						new NumericalResult("start2","",startS),
						new CategoricalResult("seq2", "", sub2.toString()),
						new NumericalResult("Score", "", pair2.getSecondElement()),
						new NumericalResult("component", "", (myComp)),
						new NumericalResult("LLR", "", llr),
						new CategoricalResult("annotation", "", (String) seq.getAnnotation()[0].getResultAt(0).getValue())
				});
				
			}
			
			idx = idx2.toArray();
			
			if(use.size() > 0){
				current = current.getPartialDataSet(use.toArray(new int[0][]));
			}else{
				break;
			}
			
			comp++;
			
		}
		
		ress.add(new ListResult("Predictions", "", null, rs));
		
		return ress;
		
	}
	
	public static LinkedList<Result> evaluateJoint(String shortName, HigherOrderHMM hmm, int nComponents1, int nComponents2, int nFirst, int nSecond, DataSet data) throws Exception{
		
				
		LinkedList<Result> ress = new LinkedList<>();
		
		ress.add(new TextResult("Model", "", 
				new FileParameter.FileRepresentation("", hmm.getGraphvizRepresentation(new DecimalFormat())), 
				"dot", shortName, null, true));
		
		
		Emission[] emission = (Emission[]) hmm.getEmissions();
		DifferentiableEmission temp = null;
		
		
		Pair<double[][][][],double[]> mods = getModels(nComponents1, nComponents2, nFirst, nSecond, emission);
		
		double[][][][] logPWMs = mods.getFirstElement();
		double[] bgMod = mods.getSecondElement();
		
		
		addSeqLogos(logPWMs, ress, 0);
		
		
		LinkedList<ResultSet> coll = new LinkedList<>();
		
		
		for(int i=0;i<data.getNumberOfElements();i++){
			Sequence seq = data.getElementAt(i);
			
			Pair<IntList,Double> pair2 = hmm.getViterbiPathFor(seq);
			String[] path = hmm.decodePath(pair2.getFirstElement());
			
			Pair<int[], Sequence[]> pair = getViterbiSeqs(path, seq);
			
			int comp1 = pair.getFirstElement()[0];
			int comp2 = pair.getFirstElement()[1];
			int startF = pair.getFirstElement()[2];
			int startS = pair.getFirstElement()[3];
			
			Sequence sub1 = pair.getSecondElement()[0];
			Sequence sub2 = pair.getSecondElement()[1];
			
			ResultSet rs = new ResultSet(new Result[]{
					new NumericalResult("Index", "", i+1),
					new NumericalResult("start1","",startF),
					new CategoricalResult("seq1", "", sub1.toString()),
					new NumericalResult("start2","",startS),
					new CategoricalResult("seq2", "", sub2.toString()),
					new NumericalResult("Score", "", pair2.getSecondElement()),
					new CategoricalResult("component", "", comp1+","+comp2),
					new NumericalResult("LLR", "", getLLR(sub1,sub2,comp1>0 ? logPWMs[0][comp1-1] : null,comp2>0 ? logPWMs[1][comp2-1] : null,bgMod)),
					new CategoricalResult("annotation", "", (String) seq.getAnnotation()[0].getResultAt(0).getValue())
			});
			coll.add(rs);
		}
		
		ListResult lr = new ListResult("Predictions", "", null, coll);
		ress.add(lr);
		
		ress.add(new StorableResult("HMM", "", hmm));
		
		return ress;
	}
	
	private static void addSeqLogos(double[][][][] logPWMs, LinkedList<Result> ress, int off) throws CloneNotSupportedException {
		
		logPWMs = ArrayHandler.clone(logPWMs);
		
		for(int i=0;i<logPWMs.length;i++){//first/second
			for(int j=0;j<logPWMs[i].length;j++){//component
				int height = 100;
				
				for(int k=0;k<logPWMs[i][j].length;k++){//position
					Normalisation.logSumNormalisation(logPWMs[i][j][k]);
				}
	
				ress.add(new PlotGeneratorResult("SeqLogo_"+(off+j+1)+"_"+(i+1), "", new SeqLogoPlotter.SeqLogoPlotGenerator(logPWMs[i][j], height), true));
			}
		}
	}

	public static Pair<double[][][][],double[]> getModels(int nComponents1, int nComponents2, int nFirst, int nSecond, Emission[] emission) throws CloneNotSupportedException{
		
		double[][][][] logPWMs = new double[2][][][];
		double[] bgMod = null;
		
		
		logPWMs[0] = new double[nComponents1][][];
		logPWMs[1] = new double[nComponents2][][];
		
		for(int j=0;j<nComponents1;j++){

			//double[][] pwm1 = new double[nFirst][];
			
			logPWMs[0][j] = new double[nFirst][];
			
			for(int i=0;i<nFirst;i++){
				DifferentiableEmission temp = ((DiscreteEmission)emission[i+1 + j*(nFirst)]).clone();
				temp.setParameterOffset(0);
				double[] pars = new double[ temp.getNumberOfParameters() ];
				temp.fillCurrentParameter(pars);
				logPWMs[0][j][i] = pars.clone();
				//Normalisation.logSumNormalisation(pars);
				//pwm1[i] = pars;
			}
					
			
		}
		
		for(int j=0;j<nComponents2;j++){

			//double[][] pwm2 = new double[nSecond][];
			
			logPWMs[1][j] = new double[nSecond][];
			

			for(int i=0;i<nSecond;i++){
				DifferentiableEmission temp = ((DiscreteEmission)emission[i+1+ nFirst*nComponents1 + j*(nSecond)]).clone();
				temp.setParameterOffset(0);
				double[] pars = new double[ temp.getNumberOfParameters() ];
				temp.fillCurrentParameter(pars);
				logPWMs[1][j][i] = pars.clone();
				//Normalisation.logSumNormalisation(pars);
				//pwm2[i] = pars;
			}

			
			
		}
		
		DifferentiableEmission temp =((DiscreteEmission)emission[0]).clone();
		temp.setParameterOffset(0);
		bgMod = new double[temp.getNumberOfParameters()];
		temp.fillCurrentParameter(bgMod);
		
		return new Pair<double[][][][],double[]>(logPWMs,bgMod);
	}
	
	public static Pair<int[], Sequence[]> getViterbiSeqs(String[] path, Sequence seq) throws Exception{
		
		int startF = -1;
		int startS = -1;
		int endF = -1;
		int endS = -1;
		int comp1 = -1;
		int comp2 = -1;
		for(int j=0;j<path.length;j++){
			if(path[j].startsWith("F")){
				if(startF < 0){
					startF = j;
					comp1 = Integer.parseInt( path[j].replaceAll("\\-.*", "").substring(1) );
				}
				endF = j;
			}
			if(path[j].startsWith("S")){
				if(startS < 0){
					startS = j;
					comp2 = Integer.parseInt( path[j].replaceAll("\\-.*", "").substring(1) );
				}
				endS = j;
			}
		}
		
		Sequence sub1 = seq.getSubSequence(startF, endF-startF+1);
		Sequence sub2 = seq.getSubSequence(startS, endS-startS+1);
		
		
		return new Pair<int[],Sequence[]>(new int[]{comp1,comp2,startF,startS},new Sequence[]{sub1,sub2});
		
	}
	
	private static double getLLR(Sequence sub1, Sequence sub2, double[][] pwm1, double[][] pwm2, double[] bgMod) {
		double score = 0;
		if(pwm1 != null){
			for(int i=0;i<pwm1.length;i++){
				score += pwm1[i][sub1.discreteVal(i)] - bgMod[sub1.discreteVal(i)];
			}
		}
		if(pwm2 != null){
			for(int i=0;i<pwm2.length;i++){
				score += pwm2[i][sub2.discreteVal(i)] - bgMod[sub2.discreteVal(i)];
			}
		}
		return score;
	}

	private HigherOrderHMM trainModel(DataSet data, int nComponents1, int nComponents2, int nStarts, int seqlen, int nFirst, 
			int nSecond, int minOffset, int minDist, int maxDist, double ess, Enum<Training> train, boolean allowSwitch, int threads, Protocol protocol) throws Exception{
		int maxOffset = seqlen - nFirst - nSecond - minDist;
		int maxDiff = seqlen - nFirst - nSecond - minOffset;
		
		DiscreteEmission insert = new DiscreteEmission(DNAAlphabetContainer.SINGLETON, ess*(seqlen - nFirst - nSecond));
		
		DifferentiableEmission[] emission = new DifferentiableEmission[nFirst*nComponents1 + nSecond*nComponents2 + 1 + 1];
		emission[0] = insert;
		for(int i=0;i<nFirst*nComponents1+nSecond*nComponents2;i++){
			if(i<nFirst*nComponents1){
				emission[i+1] = new DiscreteEmission(DNAAlphabetContainer.SINGLETON, ess/(double)nComponents1);
			}else{
				emission[i+1] = new DiscreteEmission(DNAAlphabetContainer.SINGLETON, ess/(double)nComponents2);
			}
		}
		emission[emission.length-1] = new SilentEmission();
		
		
		ArrayList<TransitionElement> tes = new ArrayList<>();
		
		int[] emissionIdx = null;
		if(maxDist == Integer.MAX_VALUE){
			emissionIdx = new int[ minOffset+ (nFirst + minDist)*(nComponents1+1) + nSecond*(nComponents2+1) + 1 + 1 ];
		}else{
			emissionIdx = new int[ minOffset+ (nFirst + maxDist)*(nComponents1+1) + nSecond*(nComponents2+1) + 1 + 1 ];
		}
			
		String[] name = new String[emissionIdx.length];
		int k=0;
		for(int i=0;i<minOffset;i++,k++){
			name[k] = "O"+i;
			emissionIdx[k] = 0;
			
			if(i==0){
				tes.add(new TransitionElement(null, new int[]{k}, new double[]{ess}));
			}else{
				tes.add(new TransitionElement(new int[]{k-1}, new int[]{k}, new double[]{ess}));
			}
		}
		
		
		int[] states = new int[nComponents1+1+1];
		double[] hypers = new double[nComponents1+1+1];
		states[0] = k-1;
		hypers[0] = ess*(maxOffset - minOffset)/2.0;
		for(int i=0;i<nComponents1+1;i++){
			if(maxDist == Integer.MAX_VALUE){
				states[i+1] = k + i*(nFirst+minDist);
			}else{
				states[i+1] = k + i*(nFirst+maxDist);
			}
			hypers[i+1] = ess/(nComponents1+1);
		}
		
		tes.add(new TransitionElement(new int[]{k-1}, states, hypers));
		
		int end = k + (nComponents1+1)*(nFirst+minDist)+(nComponents2+1)*nSecond;
		if(maxDist != Integer.MAX_VALUE){
			end = k + (nComponents1+1)*(nFirst+maxDist)+(nComponents2+1)*nSecond;
		}
		
		for(int j=0;j<nComponents1+1;j++){

			for(int i=0;i<nFirst;i++,k++){
				name[k] = "F"+j+"-"+i;
				if(j==0){
					emissionIdx[k] = 0;
				}else{
					emissionIdx[k] = 1 + (j-1)*(nFirst) + i;
				}
				
				if(i>0){
					tes.add(new TransitionElement(new int[]{k-1}, new int[]{k}, new double[]{ess}));
				}
			}
			if(maxDist == Integer.MAX_VALUE){
				for(int i=0;i<minDist;i++,k++){
					name[k] = "G"+j+"-"+i;
					emissionIdx[k] = 0;

					tes.add(new TransitionElement(new int[]{k-1}, new int[]{k}, new double[]{ess}));
				}
			}else{
				int i=0;
				for(;i<minDist-1;i++,k++){
					name[k] = "G"+j+"-"+i;
					emissionIdx[k] = 0;

					tes.add(new TransitionElement(new int[]{k-1}, new int[]{k}, new double[]{ess}));
				}
				
				for(;i<maxDist;i++,k++){
					name[k] = "G"+j+"-"+i;
					emissionIdx[k] = 0;
					if(i==minDist-1) {
						int[] states2 = new int[maxDist-i];
						double[] hyperParameters = new double[maxDist-i];
						for(int l=0;l<states2.length;l++){
							states2[l] = k+l;
							hyperParameters[l] = ess/hyperParameters.length;
						}
						
						tes.add(new TransitionElement(new int[]{k-1}, states2, hyperParameters));
					}else {
						tes.add(new TransitionElement(new int[]{k-1}, new int[] {k}, new double[] {ess}));
					}
				}
				if(!allowSwitch) {
					int next = k+nFirst*nComponents1 + maxDist*(nComponents1-j) ;
					tes.add(new TransitionElement(new int[] {k-1}, new int[] {next}, new double[] {ess}));
				}
			}
			
		}
		
		
		for(int j=0;j<nComponents1+1;j++){
			if(allowSwitch){
				if(maxDist == Integer.MAX_VALUE){
					int prev = minOffset+(j+1)*(nFirst+minDist)-1;
					int[] to = new int[nComponents2+2];
					double[] hypers2 = new double[nComponents2+2];
					to[0] = prev;
					hypers2[0] = ess*(maxDiff - minDist)/2.0;
					for(int l=0;l<nComponents2+1;l++){
						to[l+1] = (l)*nSecond + minOffset + (nFirst+minDist)*(nComponents1+1);
						hypers2[l+1] = (l==j ? ess/2.0 : ess/2.0/(nComponents2));
					}
					//System.out.println(Arrays.toString(to));
					tes.add(new TransitionElement(new int[]{prev}, to, hypers2));
				}else{
					int prev = minOffset+(j+1)*(nFirst+maxDist)-1;
					int[] to = new int[nComponents2+1];
					double[] hypers2 = new double[nComponents2+1];
					for(int l=0;l<nComponents2+1;l++){
						to[l] = (l)*nSecond + minOffset + (nFirst+maxDist)*(nComponents1+1);
						hypers2[l] = (l==j ? ess/2.0 : ess/2.0/(nComponents2));
					}
					//System.out.println(Arrays.toString(to));
					tes.add(new TransitionElement(new int[]{prev}, to, hypers2));
				}
			}else{
				if(maxDist == Integer.MAX_VALUE){
					int prev = minOffset+(j+1)*(nFirst+minDist)-1;
					tes.add(new TransitionElement(new int[]{prev}, new int[]{prev,k}, new double[]{ess*(maxDiff - minDist)/2.0,ess}));
				}/*else{
					
					int prev = minOffset+(j+1)*(nFirst+maxDist)-1;
					System.out.println(prev+" "+k);
					tes.add(new TransitionElement(new int[]{prev}, new int[]{k}, new double[]{ess*(maxDiff - maxDist)/2.0}));
				}*/
			}
		}
		
		for(int j=0;j<nComponents2+1;j++){
			for(int i=0;i<nSecond;i++,k++){
				name[k] = "S"+j+"-"+i;
				if(j==0){
					emissionIdx[k] = 0;
				}else{
					emissionIdx[k] = 1 + nFirst*nComponents1 + (j-1)*(nSecond) + i;
				}
				/*if(i==0){
					if(allowSwitch){
						int prev = minOffset+(j+1)*(nFirst+minDiff)-1;
						int[] to = new int[nComponents2+2];
						double[] hypers2 = new double[nComponents2+2];
						to[0] = prev;
						hypers2[0] = ess*(maxDiff - minDiff)/2.0;
						for(int l=0;l<nComponents2+1;l++){
							to[l+1] = (l)*nSecond + minOffset + (nFirst+minDiff)*(nComponents1+1);
							hypers2[l+1] = (l==j ? ess/2.0 : ess/2.0/(nComponents2));
						}
						//System.out.println(Arrays.toString(to));
						tes.add(new TransitionElement(new int[]{prev}, to, hypers2));
					}else{
						int prev = minOffset+(j+1)*(nFirst+minDiff)-1;
						tes.add(new TransitionElement(new int[]{prev}, new int[]{prev,k}, new double[]{ess*(maxDiff - minDiff)/2.0,ess}));
					}
				}else{*/
				if(i>0){
					tes.add(new TransitionElement(new int[]{k-1}, new int[]{k}, new double[]{ess}));
				}
			}
			
			tes.add(new TransitionElement(new int[]{k-1}, new int[]{end}, new double[]{ess}));

		}
		
		
		
		name[k] = "U";
		emissionIdx[k] = 0;
		
		
		name[k+1] = "E";
		emissionIdx[k+1] = emission.length-1;
		tes.add(new TransitionElement(new int[]{k}, new int[]{k,k+1}, new double[]{ess*0.9,ess*0.1}));
		
		boolean[] forward = new boolean[name.length];
		Arrays.fill(forward, true);
		
		//System.out.println(emissionIdx.length+" "+Arrays.toString(emissionIdx));
		for(int i=0;i<name.length;i++){
			System.out.println(i+" "+name[i]+" "+emissionIdx[i]+" "+emission[emissionIdx[i]].getClass());	
		}
		
		for(int i=0;i<tes.size();i++){
			System.out.println(tes.get(i));
		}
		
		MaxHMMTrainingParameterSet trainingParameterSet = null;
		
		if(train == Training.VITERBI){
			trainingParameterSet = new ViterbiParameterSet(nStarts, new SmallDifferenceOfFunctionEvaluationsCondition(1E-6), threads);
		}else if(train == Training.BAUMWELCH){
			trainingParameterSet = new BaumWelchParameterSet(nStarts, new SmallDifferenceOfFunctionEvaluationsCondition(1E-6), threads);
		}else{
			trainingParameterSet = new NumericalHMMTrainingParameterSet(nStarts, new SmallDifferenceOfFunctionEvaluationsCondition(1E-9), threads, Optimizer.CONJUGATE_GRADIENTS_PRP, 1E-9, 1E-6,TrainingType.VITERBI,true);
		}
		
		
		DifferentiableHigherOrderHMM hmm = new DifferentiableHigherOrderHMM(trainingParameterSet, name, emissionIdx, forward, emission, ess, tes.toArray(new TransitionElement[0]));
		
		hmm.setOutputStream(new ProtocolOutputStream(protocol, true));
		
		System.out.println(hmm.getGraphvizRepresentation(new DecimalFormat("##.###",new DecimalFormatSymbols(Locale.US))));
		
		hmm.train(data);
		
		return hmm;
	}
	
	

	@Override
	public String getToolName() {
		return "Sigma motifs";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "sigma";
	}

	@Override
	public String getDescription() {
		return "find sigma sites";
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
