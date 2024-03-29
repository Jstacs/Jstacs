/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package projects.methyl;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.Random;

import javax.naming.OperationNotSupportedException;

import de.jstacs.DataType;
import de.jstacs.Storable;
import de.jstacs.algorithms.optimization.ConstantStartDistance;
import de.jstacs.algorithms.optimization.NegativeDifferentiableFunction;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.StartDistanceForecaster;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.algorithms.optimization.termination.CombinedCondition;
import de.jstacs.algorithms.optimization.termination.IterationCondition;
import de.jstacs.algorithms.optimization.termination.MultipleIterationsCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LogGenDisMixFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.CompositeLogPrior;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.ComplementableDiscreteAlphabet;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.alphabets.GenericComplementableDiscreteAlphabet;
import de.jstacs.data.sequences.ArbitraryFloatSequence;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile;
import de.jstacs.motifDiscovery.MutableMotifDiscoverer;
import de.jstacs.motifDiscovery.MutableMotifDiscovererToolbox;
import de.jstacs.motifDiscovery.MutableMotifDiscovererToolbox.InitMethodForDiffSM;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ImageResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;
import de.jstacs.results.StorableResult;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.MarkovModelDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.UniformHomogeneousDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.localMixture.LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder;
import de.jstacs.sequenceScores.statisticalModels.differentiable.localMixture.LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder.PriorType;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.SeqLogoPlotter;
import de.jstacs.utils.SeqLogoPlotter.SeqLogoPlotGenerator;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.ToolBox.TiedRanks;
//import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.MarkovModelDiffSM;
import projects.dimont.AbstractSingleMotifChIPper;
import projects.dimont.HeuristicOneDataSetLogGenDisMixFunction;
import projects.dimont.Interpolation;
import projects.dimont.ThresholdedStrandChIPper;
import projects.encodedream.tools.MotifScores;
import projects.quickscan.QuickBindingSitePredictionTool;


public class MethylSlimDimontTool implements JstacsTool { 
	
	private static final double ALPHA = 1E-3;
	private static Random r = new Random(1123);
	
	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> parameters = new LinkedList<Parameter>();
		try{
		parameters.add(new SimpleParameter(DataType.STRING, "Alphabet", "Characters of the alphabet as a string of unseparated characters, "
				+ "first listing the symbols in forward orientation and then their complement in the same order. For instance, a methylation-aware alphabet would be specified as ACGTMH,TGCAHM and a standard DNA alphabet as ACGT,TGCA", true,"ACGTMH,TGCAHM"));
		
		parameters.add(new FileParameter("Input file", "The file name of the file containing the input sequences in annotated FastA format as generated by the Data Extractor tool", "fasta,fa,fas", true));
		
		FileParameter bgFile = new FileParameter( "Background file", "The file name of the file containing background sequences in annotated FastA format.", "fasta,fa,fas",false );
		
		SelectionParameter sp = new SelectionParameter(DataType.PARAMETERSET,new String[] {"background file","shuffled input","none"},
				new ParameterSet[]{new SimpleParameterSet(bgFile),new SimpleParameterSet(),new SimpleParameterSet()},
				"Background sample","Background sample containing negative examples, may be di-nucleotide shuffled input sequences",true);
		sp.setDefault("shuffled input");
		parameters.add( sp );
		
		parameters.add( new SimpleParameter( DataType.STRING, "Position tag", "The tag for the position information in the FastA-annotation of the input file", true,"peak" ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Value tag", "The tag for the value information in the FastA-annotation of the input file", true, "signal" ) );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Standard deviation", "The standard deviation of the position distribution centered at the position specified by the position tag", true, new NumberValidator<Double>( 1.0, 1E4 ), 75.0 ) );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Weighting factor", "The value for weighting the data, between 0 and 1", true, new NumberValidator<Double>(0.0, 1.0),0.2) );
		
		
		parameters.add( new SimpleParameter( DataType.INT, "Starts", "The number of pre-optimization runs.", true, new NumberValidator<Integer>(1,100), 20 ) );
		
		
		parameters.add( new SimpleParameter( DataType.INT, "Initial motif width", "The motif width that is used initially, may be adjusted during optimization.", true, new NumberValidator<Integer>(1,50), 20 ) );
				
		SelectionParameter modelType = new SelectionParameter(DataType.PARAMETERSET, new String[] {"LSlim model","Markov model"}, 
				new ParameterSet[] {
						new SimpleParameterSet( new SimpleParameter( DataType.INT, "Maximum distance", "The maximum distance considered in the LSlim model", true, new NumberValidator<Integer>(1,Integer.MAX_VALUE), 5 )),
						new SimpleParameterSet( new SimpleParameter( DataType.INT, "Order", "The order of the Markov model", true, new NumberValidator<Integer>(0,5), 0 ))
				}, "Model type", "The type of the motif model; a PWM model corresponds to a Markov model of order 0.", true);
		
		parameters.add( modelType );
		
		parameters.add( new SimpleParameter( DataType.INT, "Markov order of background model", "The Markov order of the model for the background sequence and the background sequence, -1 defines uniform distribution.", true, new NumberValidator<Integer>(-1,5), -1 ) );

		
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Equivalent sample size", "Reflects the strength of the prior on the model parameters.", true, new NumberValidator<Double>(0d, Double.POSITIVE_INFINITY), 4d ) );
		
		parameters.add( new SimpleParameter( DataType.BOOLEAN, "Delete BSs from profile", "A switch for deleting binding site positions of discovered motifs from the profile before searching for futher motifs.", true, true ) );
		
		parameters.add( new SimpleParameter( DataType.BOOLEAN, "Adjust for shifts", "Adjust for shifts of the motif.", true, true ) );
		
		}catch(Exception e){
			throw new RuntimeException();
		}
		
		//parameters.add( new SimpleParameter( DataType.INT, "Compute threads", "The number of threads that are use to evaluate the objective function and its gradient.", false, new NumberValidator<Integer>(1,128) ) );

		return new ToolParameterSet(getShortName(),parameters);
		
	}

	
	public AlphabetContainer getAlphabet(SimpleParameter par) throws IllegalArgumentException, DoubleSymbolException{
		String chars = (String)par.getValue();
		String[] fr = chars.split(",");
		String[] syms = fr[0].split("");
		String[] symsr = fr[1].split("");
		if(syms.length != symsr.length){
			throw new IllegalArgumentException("Not the same number of reverse complementary symbols as original symbols: "+Arrays.toString(syms)+" <-> "+Arrays.toString(symsr));
		}
		int[] compl = new int[syms.length];
		for(int i=0;i<symsr.length;i++){
			for(int j=0;j<syms.length;j++){
				if(symsr[i].equals(syms[j])){
					compl[i] = j;
				}
			}
		}
		
		return new AlphabetContainer(new GenericComplementableDiscreteAlphabet(true, syms, compl));
		
	}
	
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		
		
		progress.setLast( 1.0 );
		progress.setCurrent( 0.0 );
		
		
		AlphabetContainer con = getAlphabet((SimpleParameter)parameters.getParameterAt(0));
		
		DataSet fgData = new DataSet( con, new SparseStringExtractor( 
				new StringReader(
						((FileParameter)parameters.getParameterAt(1)).getFileContents().getContent())
				, '>', "", new SplitSequenceAnnotationParser(":",";") ) );
		
		DataSet bgData = null;

		/*if(parameters.getParameterAt(2).hasDefaultOrIsSet()){
			if("shuffled".equals( ((FileParameter)parameters.getParameterAt(2)).getFileContents().getFilename()) ){
				bgData = shuffle(fgData);
			}else{
				bgData = new DataSet( con, new SparseStringExtractor( 
						new StringReader(
								((FileParameter)parameters.getParameterAt(2)).getFileContents().getContent())
						, '>', "", new SplitSequenceAnnotationParser(":",";") ) );
			}
		}*/
		int bgSample = ((SelectionParameter)parameters.getParameterAt(2)).getSelected();
		if(bgSample == 0) {
			bgData = new DataSet( con, new SparseStringExtractor( 
					new StringReader(
							((FileParameter)((ParameterSet)((SelectionParameter)parameters.getParameterAt(2)).getValue()).getParameterAt(0)).getFileContents().getContent())
					, '>', "", new SplitSequenceAnnotationParser(":",";") ) );
		}else if(bgSample == 1){
			bgData = shuffle(fgData);
		}else {
			bgData = null;
		}
		
		
		String position = parameters.getParameterAt(3).getValue().toString();
		String value = parameters.getParameterAt(4).getValue().toString();
		double sd = (Double) parameters.getParameterAt(5).getValue();
		double wf = (double) parameters.getParameterAt(6).getValue();
		int restarts = (Integer) parameters.getParameterAt(7).getValue();
		int motifLength = (Integer) parameters.getParameterAt(8).getValue();
		int modelType = (Integer) ((SelectionParameter)parameters.getParameterAt(9)).getSelected();
		int ordDist =  (int) ((ParameterSet)((SelectionParameter)parameters.getParameterAt(9)).getValue()).getParameterAt(0).getValue();
		int bgOrder = (Integer) parameters.getParameterAt(10).getValue();
		double ess = (Double) parameters.getParameterAt(11).getValue();
		boolean delete = (Boolean) parameters.getParameterAt(12).getValue();
		boolean modify = (Boolean) parameters.getParameterAt(13).getValue();
		
		
		double filterThreshold = 0.3;
		double filterThresholdEnd = 0.3;
		
	
		boolean free = false;
		
		DataSet[] data = { fgData };


		Sequence[] annotated = new Sequence[data[0].getNumberOfElements() + (bgData==null?0:bgData.getNumberOfElements())];
		double[][] weights = new double[2][];
		double[] raw =new double[data[0].getNumberOfElements()];
		
		//read annotation
		double[] mean = new double[annotated.length];
		Arrays.fill( mean, Double.NaN );
		for( int j = 0; j < raw.length; j++ ) {
			annotated[j] = data[0].getElementAt(j);
			SequenceAnnotation[] seqAn = annotated[j].getAnnotation();
			mean[j] = Double.NaN;
			for( int i = 0; i < seqAn.length; i++ ) {
				if( seqAn[i].getType().equals(value) ) {
					raw[j] = Double.parseDouble( seqAn[i].getIdentifier() );
				} else if( seqAn[i].getType().equals(position) ) {
					mean[j] = Double.parseDouble( seqAn[i].getIdentifier() );
				}
			}
		}
		if( bgData != null ) {
			for( int j = 0; j < bgData.getNumberOfElements(); j++ ) {
				Sequence s = bgData.getElementAt(j);
				annotated[raw.length+j] = s;
				mean[raw.length+j] = s.getLength()/2d;
			}
		}
		
		

		double[] w = Interpolation.getWeight( data[0], raw, wf, Interpolation.RANK_LOG );
		if( bgData == null ) {
			weights[0] = w;
		} else {
			weights[0] = new double[w.length + bgData.getNumberOfElements()];
			System.arraycopy(w, 0, weights[0], 0, w.length);
		}
		weights[1] = Interpolation.getBgWeight( weights[0] );
		
		if(bgData != null) {
			double fac2 = ((double)fgData.getNumberOfElements())/((double)bgData.getNumberOfElements());
			for(int i=w.length;i<weights.length;i++){
				weights[1][i] *= fac2;
			}
		}
		
		boolean[][] allowed = new boolean[annotated.length][];//TODO Jan? BitSets?
		for(int i=0;i<annotated.length;i++){
			allowed[i] = new boolean[annotated[i].getLength()];
			Arrays.fill( allowed[i], true );
		}
		
		//create initial annotation
		data[0] = annotate(annotated,weights,mean,sd,allowed, fgData.getNumberOfElements());
		
		//complete vs small fg data/weights
		DataSet completeData = data[0];
		double[][] completeWeight = { weights[0], weights[1] };
		
		//create models
		DifferentiableStatisticalModel motif = null;
		if(modelType == 1){
				motif = new MarkovModelDiffSM( con, motifLength, ess, true, ordDist, null );
		}else{
			motif = new LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder( con, motifLength, 1, ordDist, ess, 0.9, PriorType.BDeu );
			//motif = new LimitedSparseLocalInhomogeneousMixtureDiffSM( con, motifLength, 2, -fgOrder, ess, 0.9, priorType );
			//motif = new SparseLocalInhomogeneousMixtureDiffSM( con, motifLength, -fgOrder, ess, false, 0.9, PriorType.Complex_Mixture );
		}
		DifferentiableStatisticalModel fg =	new ThresholdedStrandChIPper( 1, 0.5, motif );
		fg.initializeFunctionRandomly( false );

		double fac = (1-wf)/wf;
		
		DifferentiableStatisticalModel bg = getBgSF( con, bgOrder, ess*fac, data[0].getAverageElementLength() );
		bg.initializeFunction(0, false, data, weights );

		DifferentiableStatisticalModel[] score = { fg, bg };
		double[] beta = LearningPrinciple.getBeta( LearningPrinciple.MSP );
		
		//prepare for optimization
		LearningPrinciple initKey = (beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX]>0) ? LearningPrinciple.MCL : LearningPrinciple.ML;	
			
		HeuristicOneDataSetLogGenDisMixFunction initObjective = new HeuristicOneDataSetLogGenDisMixFunction( threads, score, completeData, completeWeight.clone(), new CompositeLogPrior(), LearningPrinciple.getBeta( initKey ), true, free );
		HeuristicOneDataSetLogGenDisMixFunction objective = new HeuristicOneDataSetLogGenDisMixFunction( threads, score, data[0], weights, new CompositeLogPrior(), beta, true, free );
		NegativeDifferentiableFunction neg = new NegativeDifferentiableFunction( objective );
		
		double eps = 1E-4;
		AbstractTerminationCondition stop = new CombinedCondition( 2,
				new MultipleIterationsCondition( 5, new SmallDifferenceOfFunctionEvaluationsCondition(eps) ),
				new IterationCondition(100) 
		);
		StartDistanceForecaster start = new ConstantStartDistance(1);
		double[] params = null;
		double[][] pwm = null;
		double[] entropy = new double[motifLength], kl = new double[motifLength];
		SignificantMotifOccurrencesFinder smof;
		byte algo = Optimizer.CONJUGATE_GRADIENTS_PRP;
		GenDisMixClassifierParameterSet genDisMixParams = new GenDisMixClassifierParameterSet( con, 0, algo, eps, eps*1E-1, 1, free, KindOfParameter.PLUGIN, true, threads);
		objective.reset( score );
		DataSet smallData = null;
		double[][] smallWeight = new double[2][];
		double[] p = null;
		
		//create small data set
		Pair<DataSet,double[][]> small = getSmallDataSets( completeWeight, annotated, 0.3, 1000 );
		smallData = small.getFirstElement();
		smallWeight = small.getSecondElement();
			
		initObjective.setDataAndWeights( new DataSet[]{smallData}, smallWeight );
		//initObjective.setDataAndWeights( data, weights );
		initObjective.reset( score );
		
		double percentKmers = 1.0;
		
		//get (good) start paramaters
		ComparableElement<double[], Double>[] array2, sortedPars = new ComparableElement[restarts];
		int numKmers = 0;
		if( percentKmers > 0.0 ) {		
			int k=7;
			ComparableElement<String, Double>[] array = getKmereSequenceStatistic(Math.max( 50, (int)Math.ceil( percentKmers*restarts )), k, smallData, smallWeight[0]);
			array2 = new ComparableElement[array.length];
			int a = 4;
			double d = 0.1/(a-1), h;
			d = (1d-a*d)/(a*d);
			h = d * motif.getESS();
			for( int z = 0; z < array.length; z++ ) {
				//((AbstractSingleMotifChIPper)score[0]).initializeMotif( 0, new DataSet("", Sequence.create(con, array[z].getElement()) ), new double[]{h} );
				Pair<DataSet,double[]> temp = getInitData( smallData, smallWeight, Sequence.create(con, array[z].getElement()) , motifLength );
				((AbstractSingleMotifChIPper)score[0]).initializeMotif( 0, temp.getFirstElement(), temp.getSecondElement() );
				p = objective.getParameters(KindOfParameter.PLUGIN);
			
				initObjective.reset(score);
				initObjective.resetHeuristics();
				double val = initObjective.evaluateFunction(p);//TODO
				//double val = array[z].getWeight();
				//array[z] = new ComparableElement<String, Double>( array[z].getElement(), val );
				
				array2[z] = new ComparableElement<double[], Double>( p, val );
			}
			Arrays.sort( array2 );
			numKmers = Math.min( array2.length, (int)Math.ceil(restarts*percentKmers) );
			for(int z=0;z<numKmers;z++){
				sortedPars[z] = array2[array2.length-1-z];
			}
			
		}
		
		
		if( numKmers != sortedPars.length ) {
			InitMethodForDiffSM[] initMeth = {InitMethodForDiffSM.PLUG_IN, InitMethodForDiffSM.NOTHING};
			ComparableElement<double[], Double>[] sortedPars2 = MutableMotifDiscovererToolbox.getSortedInitialParameters( score, initMeth, initObjective, Math.max( 100, restarts ), SafeOutputStream.getSafeOutputStream( null ), 0 );
			for(int z=0;z<sortedPars.length-numKmers;z++){
				sortedPars[numKmers+z] = sortedPars2[sortedPars2.length-1-z];
			}
			Arrays.sort( sortedPars );
		}
		
		progress.setCurrent(0.05);
		
		//preoptimization
		AbstractTerminationCondition stop2 = new CombinedCondition( 2,
				new MultipleIterationsCondition( 5, new SmallDifferenceOfFunctionEvaluationsCondition(eps) ),
				new IterationCondition(25) 
		);
		ComparableElement<double[], Double>[] preOpt = new ComparableElement[restarts];
		for( int r = 0; r < restarts; r++ ) {
			data[0] = smallData;
			weights = smallWeight;
			objective.setDataAndWeights( data, weights );
			objective.resetHeuristics();
			
			protocol.append("-----------------------------------------\npre-optimization " + r+"\n");
			
			start.reset();
			
			p = sortedPars[sortedPars.length-1-r].getElement();			
			objective.setParams(p);
			//System.out.println(score[0]);
			Optimizer.optimize( algo, neg, p, stop, eps*1E-1, start, null );
			
			data[0] = completeData;
			weights = completeWeight;
			objective.setDataAndWeights( data, weights );
			preOpt[r] = new ComparableElement<double[], Double>( p, objective.evaluateFunction(p) );
			
			//out.writeln("consensus: "+getConsensus( DNAAlphabetContainer.SINGLETON, (((MarkovModelDiffSM)((AbstractSingleMotifChIPper) score[0]).getFunction( 0 ) ).getPWM())));
			protocol.append("model: "+((AbstractSingleMotifChIPper) score[0]).getFunction( 0 )+"\n");
			protocol.append("score: "+preOpt[r].getWeight()+"\n" );
			//System.out.println(score[0]);
			((AbstractSingleMotifChIPper)score[0]).resetPositions();
			progress.setCurrent( 0.05 + ((r+1.0)/(double)restarts)*0.2 );
		}
		protocol.append("-----------------------------------------\n");
		progress.setCurrent(0.25);
		
		//filter out redundant motifs
		Arrays.sort( preOpt );
		ArrayList<ComparableElement<double[], Double>> list = filter2( (AbstractSingleMotifChIPper) score[0], preOpt, smallData, filterThreshold ,motifLength, protocol );

		//final optimization
		MutableMotifDiscoverer[] best = new MutableMotifDiscoverer[list.size()];		
		Storable[] storables = new Storable[best.length];
		double[] opts = new double[best.length];
		Pair<double[][][],int[][]>[] pairs = new Pair[best.length];		
		for( int m = 0; m < best.length; m++ ) {
			protocol.append("+++++++++++++++++++++++++++++++++++++++++++++++++++" +
					"\n\nfinalOpt " + m + " -----------------------------------------\n");
			
			best[m] = (MutableMotifDiscoverer)score[0];
			for( int j = 0; j < best[m].getNumberOfMotifs(); j++ ) {
				if( best[m].getMotifLength(j) != motifLength ) {
					best[m].modifyMotif( j, 0, motifLength-best[m].getMotifLength(j) );
				}
				//System.out.println( "motif " + j + ": " + best[m].getMotifLength(j));
			}
			((AbstractSingleMotifChIPper) score[0]).reset();
			((AbstractSingleMotifChIPper) score[0]).resetPositions();
			objective.reset( score );
			objective.resetHeuristics();
			
			start.reset();
			p = list.get(m).getElement();
			objective.setParams(p);
			//System.out.println("before "+getConsensus( DNAAlphabetContainer.SINGLETON, (((MarkovModelDiffSM)((AbstractSingleMotifChIPper) score[0]).getFunction( 0 ) ).getPWM())));
			
			//optimize on the complete data set
			data[0] = annotate( annotated, weights, mean, sd, allowed, fgData.getNumberOfElements() );
			objective.setDataAndWeights( data, weights );
			Optimizer.optimize( algo, neg, p, stop, eps*1E-1, start, SafeOutputStream.getSafeOutputStream( null ) );

			progress.setCurrent(0.25 + 0.7*((m+0.5)/(double)best.length));
			
			//heuristic for motif length
			double[] sds = new double[1];
			heuristic((MutableMotifDiscoverer) score[0], completeData, completeWeight, objective,mean,sds, protocol, modify);
			
			
			//final optimization
			((AbstractSingleMotifChIPper)score[0]).reset();
			
			double newSd = Math.sqrt( sds[0]*sd);
			if(Double.isNaN( newSd ) || Double.isInfinite( newSd ) || newSd <= 0){
				protocol.append("Did not adjust sd to "+newSd+" using "+sds[0]+" and "+sd+"\n");
				newSd = sd;
			}
			
			data[0] = annotate( annotated, weights, mean, newSd, allowed, fgData.getNumberOfElements() );
			objective.setDataAndWeights( data, weights );
			((AbstractSingleMotifChIPper) score[0]).resetPositions();
			objective.reset( score ); // reset heuristics: toBeUsed..., ...
			
			p = objective.getParameters(KindOfParameter.LAST);
			objective.setParams(p);
			Optimizer.optimize( algo, neg, p, stop, eps*1E-1, start, SafeOutputStream.getSafeOutputStream( null ) );
			
			protocol.append("model: "+((AbstractSingleMotifChIPper) score[0]).getFunction( 0 )+"\n");
			//out.writeln("consensus: "+getConsensus( DNAAlphabetContainer.SINGLETON, (((MarkovModelDiffSM)((AbstractSingleMotifChIPper) score[0]).getFunction( 0 ) ).getPWM())));
			
			//save information (xml, model, predictions)
			best[m] = (MutableMotifDiscoverer) score[0].clone();
			opts[m] = objective.evaluateFunction(p);
			
			GenDisMixClassifier cl = new GenDisMixClassifier( genDisMixParams, new CompositeLogPrior(), opts[m], LearningPrinciple.getBeta(LearningPrinciple.MSP), score );
			cl.setClassWeights( false, objective.getClassParams(p) );
			
			storables[m] = cl;
						
			smof = new SignificantMotifOccurrencesFinder(best[m],completeData,completeWeight[1],ALPHA);			
			Pair<double[][][],int[][]> pair = smof.getPWMAndPositions( 0, completeData, completeWeight[0], 0, 0 );
			pairs[m] = pair;
			
			//"delete" positions from the profile 
			if( delete && m+1 < best.length ) {
				delete(pair.getSecondElement(),allowed,motifLength);
			}
			
			progress.setCurrent(0.25 + 0.7*((m+1.0)/(double)best.length));
		}
		
		//rank and filter final motifs
		int[] o = ToolBox.rank( opts, TiedRanks.IN_ORDER );
		int[] index = new int[o.length];
		for( int i = 0; i<o.length; i++ ) {
			index[o[i]] = i;
		}
		boolean[] use = postFilter( best, index, smallData, filterThresholdEnd, motifLength );
		
		LinkedList<ResultSetResult> results = new LinkedList();
		
		//write final xmls and sequence logos
		for(int m=0,n=0;m<best.length;m++){
			if(use[m]){
				LinkedList<Result> result = new LinkedList<Result>();
				
				result.add(new StorableResult( "SlimDimont "+(n+1), "The SlimDimont classifier", storables[index[m]] ) );
				
				LinkedList<Sequence> bs = new LinkedList<Sequence>();
				DoubleList bsWeights = new DoubleList();
				
				//TODO Jens->Jan: fgData vs. completeData: pairs[index[m]]
				result.add(getListResult(fgData, completeWeight[0],pairs[index[m]], ((ThresholdedStrandChIPper)((GenDisMixClassifier)storables[index[m]]).getDifferentiableSequenceScore( 0 )).getMotifLength( 0 ), n, ((ThresholdedStrandChIPper)((GenDisMixClassifier)storables[index[m]]).getDifferentiableSequenceScore( 0 )).getMotifModel(), bs, bsWeights, value ));
				
				pwm = pairs[index[m]].getFirstElement()[0];
				
				if(!Double.isNaN( pwm[0][0] )){
					try{
						result.add(new PlotGeneratorResult("Motif "+(n+1), "Sequence logo of motif "+(n+1), new SeqLogoPlotGenerator(pwm,1000), true));
						result.add(new PlotGeneratorResult("Motif "+(n+1)+" (rc)", "Sequence logo of the reverse complement of motif "+(n+1), new SeqLogoPlotGenerator( PFMComparator.getReverseComplement((ComplementableDiscreteAlphabet) con.getAlphabetAt(0), pwm) ,750), true));
						/*int height = SeqLogoPlotter.getHeight( 750, pwm );
						result.add(new ImageResult( "Motif "+(n+1), "Sequence logo of motif "+(n+1), SeqLogoPlotter.plotLogoToBufferedImage( height, pwm ) ));
						result.add(new ImageResult( "Motif "+(n+1)+" (rc)", "Sequence logo of the reverse complement of motif "+(n+1), SeqLogoPlotter.plotLogoToBufferedImage( height, PFMComparator.getReverseComplement( DNAAlphabet.SINGLETON, pwm ) ) ));
						*///result.add( new ImageResult( "Dependency logo "+(n+1), "Dependency logo of motif "+(n+1), SeqLogoPlotter.plotDefaultDependencyLogoToBufferedImage( new DataSet("",bs), bsWeights.toArray(), 600 ) ) );
					}catch(Exception e){
						e.printStackTrace();
					}
				}
				ResultSetResult rsr = new ResultSetResult("Motif "+(n+1), "The Dimont results for motif "+(n+1), null, new ResultSet(result));
				results.add(rsr);
				n++;
			}
		}
		initObjective.stopThreads();
		objective.stopThreads();
		progress.setCurrent(1.0);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(results), parameters, getToolName(), new Date(System.currentTimeMillis()));

		
	}
	
	
	private DataSet shuffle(DataSet fgData) throws Exception {
		Sequence[] seqs = new Sequence[fgData.getNumberOfElements()];
		for(int i=0;i<seqs.length;i++){
			seqs[i] = shuffle(fgData.getElementAt(i),2);
		}
		
		return new DataSet("",seqs);
	}
	
	public static Sequence shuffle( Sequence original, int k /*, boolean randomRotation*/ ) throws Exception {
		//Random r = new Random();
		//r.setSeed(1331);//TODO added for consistency, 12/19/2017
		
		int n = original.getLength();
		int[] shuffle = new int[n], help = shuffle.clone();
		for( int i = 0; i < n; i++ ) {
			shuffle[i] = original.discreteVal(i);
			//System.out.print( ( i % 10 == 0 ) ? "." : " " );
		}
		//System.out.println();
				 
		int anz = 0;
		for( int i = 0; i < shuffle.length; i++ ) {
			int a = r.nextInt(n-4*(k+1));
			int b = a+1+r.nextInt(n-3*(k+1)-a);
			int c = b+1+r.nextInt(n-2*(k+1)-b);
			int d = c+1+r.nextInt(n-k+1-c);
			
			boolean matches;
			int j = 0;
			//System.out.println();
			while( j < k-1 && shuffle[a+j] == shuffle[c+j] ) {
				//System.out.println( j + "\t" + (a+j) + "\t" + shuffle[a+j] + "\t" + (c+j) + "\t" + shuffle[c+j] );
				j++;
			}
			matches = j == k-1;
			j = 0;
			while( matches && j < k-1 && shuffle[b+j] == shuffle[d+j] ) {
				//System.out.println( j + "\t" + (b+j) + "\t" + shuffle[b+j] + "\t" + (d+j) + "\t" + shuffle[d+j] );
				j++;
			}
			matches &= j == k-1;
			if( matches ) {
				anz++;
				//System.out.println( i + "\t" + a + "\t" + b + "\t" + c + "\t" + d + "\t" + anz );
				Arrays.fill(help, 0);
				System.arraycopy( shuffle, 0, help, 0, a );
				System.arraycopy( shuffle, c, help, a, d-c );
				System.arraycopy( shuffle, b, help, a+d-c, c-b );
				System.arraycopy( shuffle, a, help, a+d-b, b-a );
				System.arraycopy( shuffle, d, help, d, n-d );
				System.arraycopy( help, 0, shuffle, 0, n );
			}
		}
		
		Sequence seq = new IntSequence( original.getAlphabetContainer(), shuffle );
		
		return seq;
	}


	public static ListResult getListResult( DataSet data, double[] weights, Pair<double[][][], int[][]> pair, int motifLength, int motifIndex, DifferentiableStatisticalModel model, LinkedList<Sequence> bs, DoubleList bsWeights, String valueKey ) throws Exception {
		
		SplitSequenceAnnotationParser pars = new SplitSequenceAnnotationParser( ":", ";" );
		
		LinkedList<ResultSet> set = new LinkedList<ResultSet>();
		int[][] pos = pair.getSecondElement();
		double[][] pvals = pair.getFirstElement()[1];
		for(int i=0;i<data.getNumberOfElements();i++){//TODO Jens->Jan: pos.length vs. data.getNumberOfElements()
			for(int j=0;j<pos[i].length;j++){
				int curr = pos[i][j];
				boolean rc = false;
				if(curr < 0){
					curr = -curr-1;
					rc = true;
				}
				Sequence sub = data.getElementAt( i ).getSubSequence( curr, motifLength );
				Sequence sub2 = sub;
				if(rc){
					sub2 = sub.reverseComplement();
				}
				double score = model.getLogScoreFor( sub2 );
				
				ResultSet rs = new ResultSet( new Result[]{ 
				                                           new NumericalResult( "Sequence index", "The index of the sequence", i+1 ),
				                                           new NumericalResult( "Position", "The starting position of the motif within the sequence", curr+1 ),
				                                           new CategoricalResult( "Strand", "The strand of the predicted BS", rc ? "-" : "+" ),
				                                           new NumericalResult( "p-value", "The p-value of the predicted BS", pvals[i][j] ),
				                                           new NumericalResult( "-log10(p-value)", "The negative logarithm of the p-value of the predicted BS", -Math.log10( pvals[i][j] ) ),
				                                           new NumericalResult( "Score", "The model score of the predicted BS", score ),
				                                           new CategoricalResult( "Binding site", "The binding site as in the sequence", sub.toString() ),
				                                           new CategoricalResult( "Adjusted binding site", "The binding site in predicted orientation", sub2.toString() ),
				                                           new CategoricalResult( "Signal", "The signal of the sequence annotation", data.getElementAt( i ).getSequenceAnnotationByType( valueKey, 0 ).getIdentifier() ),
				                                           new CategoricalResult( "Sequence annotation", "The annotation of the original sequence", pars.parseAnnotationToComment( ' ', data.getElementAt( i ).getAnnotation() ).substring( 1 ) )
				} );
				set.add( rs );
				bs.add(sub2);
				bsWeights.add(-Math.log10( pvals[i][j] ));
			}
		}
		
		ListResult lr = new ListResult( "Predictions for motif "+(motifIndex+1), "", null, set.toArray( new ResultSet[0] ) );
		return lr;
	}

	private static void delete( int[][] positions, boolean[][] allowed, int length ) {
		for(int i=0;i<positions.length;i++){
			for(int j=0;j<positions[i].length;j++){
				int pos = positions[i][j];
				if(pos < 0){
					pos = -pos - 1;
				}
				Arrays.fill( allowed[i], 
						Math.max( 0, pos-length/2), 
						Math.min( allowed[i].length, pos + length/2 ),
						false );
			}
		}		
	}
	
	private static Pair<DataSet,double[][]> getSmallDataSets(double[][] completeWeight, Sequence[] data, double percent, int maxN) throws EmptyDataSetException, WrongAlphabetException{
		int[] ofg = ToolBox.order( completeWeight[0], false );
		int[] obg = ToolBox.order( completeWeight[1], false );
		boolean[] used = new boolean[completeWeight[0].length];
		
		double[] sums = {ToolBox.sum( completeWeight[0] ), ToolBox.sum( completeWeight[1] )};
		
		
		double[] currSums = new double[2];
		int idxFg = ofg.length;
		int idxBg = obg.length;
		int totN = 0;
		int nfg = 0;
		
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		DoubleList w = new DoubleList();
		
		while(totN < maxN && currSums[0] + completeWeight[0][ofg[idxFg-1]] < sums[0]*percent){
			//System.out.println((idxFg-1)+" "+completeWeight[0][ofg[idxFg-1]]);
			currSums[0] += completeWeight[0][ofg[idxFg-1]];
			if(!used[ofg[idxFg-1]]){
				seqs.add( data[ofg[idxFg-1]] );
				w.add( completeWeight[0][ ofg[idxFg-1] ] );
				used[ofg[idxFg-1]] = true;
				totN++;
				nfg++;
			}
			idxFg--;
			double tempPerc = currSums[0]/sums[0];
			while(totN < maxN && currSums[1] + completeWeight[1][obg[idxBg-1]] < sums[1]*tempPerc){
				currSums[1] += completeWeight[1][obg[idxBg-1]];
				if(!used[obg[idxBg-1]]){
					seqs.add( data[obg[idxBg-1]] );
					w.add( completeWeight[0][ obg[idxBg-1] ] );
					used[obg[idxBg-1]] = true;
					totN++;
				}
				idxBg--;
			}
			//System.out.println(totN+" "+(ofg.length-idxFg)+" "+(obg.length-idxBg)+" "+tempPerc);
		}
		double[] rw = w.toArray();
		return new Pair<DataSet, double[][]>( new DataSet( "", seqs.toArray( new Sequence[0] ) ), new double[][]{rw, Interpolation.getBgWeight( rw )} );
	}

	static DataSet annotate( Sequence[] annotated, double[][] weights, double[] mean, double sd, boolean[][] allowed, int fg ) throws WrongAlphabetException, WrongSequenceTypeException, EmptyDataSetException {
		AlphabetContainer ref = new AlphabetContainer( new ContinuousAlphabet() );
		float[][] histogram = new float[annotated.length][];
		for( int j = 0; j < weights[0].length; j++ ) {
			histogram[j] = new float[annotated[j].getLength()];
			if( j >= fg ) {
				float sum = 0;
				for( int i = 0; i < histogram[j].length; i++ ) {
					if(allowed[j][i]){
						histogram[j][i] = 1;
						sum++;
					}
				}
				for(int i=0;i<histogram[j].length;i++){
					histogram[j][i] /= sum;
				}
			} else {
				float max = 0, sum = 0;
				for( int i = 0; i < histogram[j].length; i++ ) {
					if(allowed[j][i]){
						histogram[j][i] = (float) ((i - mean[j]) / sd);
						histogram[j][i] = (float) Math.exp( -0.5 * histogram[j][i]*histogram[j][i] );
						sum += histogram[j][i];
					}
				}
				for( int i = 0; i < histogram[j].length; i++ ) {
					histogram[j][i] /= sum;
					if( histogram[j][i] > max ) {
						max = histogram[j][i];
					}
				}
			
				float[] histBg = histogram[j].clone();
				sum = 0;
				for(int i=0;i<histBg.length;i++){
					if(allowed[j][i]){
//Jens->Jan: if real bg exists then mixture Gau� and uniform?	
						if(fg<weights[0].length){
							histBg[i] = 1;
						}else{
							histBg[i] = max - histBg[i];
						}
						sum += histBg[i];
					}
				}
				for(int i=0;i<histBg.length;i++){
					histBg[i] /= sum;
				}
				
				for( int i = 0; i < histogram[j].length; i++ ) {
					histogram[j][i] = (float)(weights[0][j] * histogram[j][i] + weights[1][j] * histBg[i]);
				}
			}
			annotated[j] = annotated[j].annotate( false, new ReferenceSequenceAnnotation( "reads", new ArbitraryFloatSequence( ref, histogram[j] ) ) );
		}
		return new DataSet( "", annotated );
	}

	private static Pair<DataSet,double[]> getInitData(DataSet completeData, double[][] completeWeight, Sequence kmer, int length) throws Exception {
		
		//System.out.println("    "+kmer);
		
		int off = (length-kmer.getLength())/2;//TODO even!
		
		int evCorr = (length-kmer.getLength()) % 2 == 0 ? 0 : 1;
		
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		DoubleList ws = new DoubleList();
		
		for(int i=0;i<completeData.getNumberOfElements();i++){
			Sequence seq = completeData.getElementAt( i );
			//System.out.println(seq.getLength());
			for(int j=off;j<seq.getLength()-off-kmer.getLength();j++){
				int dist = getHammingDistance(kmer,seq,j);
				if(dist < 3 && dist > -4 && (dist>=0 || j-off-evCorr>=0)){
					if(dist >=0){
						seqs.add( seq.getSubSequence( j-off, length ) );
					}else{
						seqs.add( seq.getSubSequence( j-off-evCorr, length ).reverseComplement() );
					}

					//System.out.println(seqs.getLast()+" "+dist);
					if(dist < 0){
						dist = -(dist+1);
					}
					ws.add( completeWeight[0][i]/Math.pow( 4, dist ) );
				}
			}			
		}
		//System.out.println(seqs);
		return new Pair<DataSet, double[]>( new DataSet( "", seqs ), ws.toArray() );
	}
	
	private static int getHammingDistance(Sequence kmer, Sequence seq, int start) throws WrongAlphabetException, OperationNotSupportedException{
	//	System.out.println(seq.getLength()+" "+start);
		int distfwd = kmer.getHammingDistance( seq.getSubSequence( start, kmer.getLength() ) );
		int distrc = kmer.getHammingDistance( seq.getSubSequence( start, kmer.getLength() ).reverseComplement() );
		if(distfwd <= distrc){
			return distfwd;
		}else{
			return -distrc-1;
		}
	}
	
	private static boolean heuristic( MutableMotifDiscoverer best, DataSet completeData, double[][] completeWeight, LogGenDisMixFunction objective, double[] mean, double[] sds, Protocol protocol, boolean modify ) throws Exception {
		SignificantMotifOccurrencesFinder smof = new SignificantMotifOccurrencesFinder(best,completeData,completeWeight[1],ALPHA);//XXX Jan
		boolean modified = false;
		double log2 = Math.log( 2 );
		for(int im=0;im<best.getNumberOfMotifs();im++){
			int add = 5;
			if(completeData.getAverageElementLength() - best.getMotifLength( im ) < 20){
				add = (int)((completeData.getAverageElementLength() - best.getMotifLength( im ))/4.0);
			}
			//add = 0;//TODO
			//System.out.println("add: "+add);
			LinkedList<Sequence> bs = new LinkedList<Sequence>();
			DoubleList bsWeights = new DoubleList();
			DoubleList bsScores = new DoubleList();
			Pair<double[][],double[]> pair = smof.getPWMAndPosDist( im, completeData, completeWeight[0], mean, add, add, bs, bsWeights, bsScores );
			double[][] pwm = pair.getFirstElement();
			sds[im] = pair.getSecondElement()[0];
			double[] entropy = new double[pwm.length], kl = new double[pwm.length];
			double[] bgDistr = getCounts(completeData,completeWeight[1]);//PFMComparator.getCounts(completeData);
			PFMComparator.normalize( bgDistr );

			for( int i = 0; i < pwm.length; i++ ) {
				kl[i] = 0;
				entropy[i] = Math.log( pwm[i].length )/log2;
				for( int j = 0; j < pwm[i].length; j++ ) {
					if( pwm[i][j] > 0 ) {
						entropy[i] += pwm[i][j] * Math.log( pwm[i][j] )/log2;
						//kl[i] += (pwm[i][j]-bgDistr[j]) * Math.log( pwm[i][j]/bgDistr[j] );
						kl[i] += pwm[i][j] * Math.log( pwm[i][j]/bgDistr[j] );
					}
				}
			}
			if(bs.size() == 0){
				return false;
			}
			//System.out.println("kl: "+Arrays.toString( kl ));
			DataSet bsDs = new DataSet( "", bs );
			//double[][] mi2 = TwoPointEvaluater.getMIInBits( bsDs, bsWeights.toArray() );
			/*for(int i=0;i<mi2.length;i++){
				System.out.println(i+" "+Arrays.toString( mi2[i] ));
			}*/
			double[][] mi = getMaxMI(bsDs, bsWeights.toArray(), bsScores.toArray());
			/*for(int i=0;i<mi.length;i++){
				System.out.println(i+" "+Arrays.toString( mi[i] ));
			}*/
			int maxOrder = 0;
			
			ThresholdedStrandChIPper chipper = (ThresholdedStrandChIPper)best;
			DifferentiableStatisticalModel motifModel = chipper.getMotifModel();
			if(motifModel instanceof MarkovModelDiffSM){
				maxOrder = ((MarkovModelDiffSM)motifModel).getOrder();
			}else if(motifModel instanceof LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder){
				maxOrder = ((LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder)motifModel).getDistance() + ((LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder)motifModel).getOrder() -1;
			}
			
		//	double[] maxMi = new double[mi.length];
			for(int i=0;i<mi.length;i++){
				for(int j=Math.max( 0, i-maxOrder );j<Math.min( mi[i].length, i+maxOrder+1 );j++){
					if(i != j){
						if(mi[i][j] > kl[i]){
							kl[i] = mi[i][j];
						}
					}
				}
				/*for(int j=0;j<mi[i].length;j++){
					if(i != j){
						if(mi[i][j] > maxMi[i]){
							maxMi[i] = mi[i][j];
						}
					}
				}*/
			}
			//System.out.println("kl: "+Arrays.toString( kl ));
			//System.out.println(bsDs.getNumberOfElements()+": "+bsDs.getAverageElementLength()+", "+bsDs.getElementLength());
			//System.out.println(Arrays.toString( maxMi ));
			int left = -add, right = add;
			double thresh = 0.2;//full: 0.0, normal: 0.1, Dimont: 0.2
			
			if(!modify){
				return false;//TODO
				
			}
			
			while( left+add < kl.length && kl[left+add] < thresh /*&& maxMi[left+add] < thresh*/ ){
				left++;
			}
			
			while(right-add > -kl.length && kl[kl.length-1+right-add] < thresh /*&& maxMi[kl.length-1+right-add] < thresh*/){
				right--;
			}
			
			protocol.append( "left: "+ left + ", right: " + right+"\n" );
			//if(!modify){//TODO start
				if(left <= 0 && right < 0 ){
					if(left <= right){
						left = right;
					}else{
						left = right = (left+right)/2;
					}
					//right = left = Math.max( left, right );
				}else if(right >= 0 && left > 0){
					if(right >= left){
						right = left;
					}else{
						left = right = (left+right)/2;
					}
					//left = right = Math.min( right, left );
				}else{
					left = right = 0;
				}
//			}//TODO end
			
			
				protocol.append( "left: "+ left + ", right: " + right+"\n" );
			if( left+add == kl.length || right-add == -kl.length ) {
				protocol.append( "tried to remove the complete motif: no modifications\n");
			} else {
				double normOld = ((DifferentiableStatisticalModel) best).getLogNormalizationConstant();
				if( left != 0 || right!= 0 ) {
					if( best.modifyMotif( im, left, right ) ){
						double w = normOld - ((DifferentiableStatisticalModel) best).getLogNormalizationConstant();
						//System.out.println(w);
						objective.addTermToClassParameter( 0, w );
					}
					protocol.append("modified motif\n");
					modified = true;
				} else {
					protocol.append("no modifications for the motif\n");
				}
			}
		}
		return modified;
	}
	
	private static double[][] getMaxMI( DataSet bsDs, double[] weights, double[] scores ) {
		ComparableElement<Pair<Sequence,Double>, Double>[] els = new ComparableElement[bsDs.getNumberOfElements()];
		for(int i=0;i<bsDs.getNumberOfElements();i++){
			els[i] = new ComparableElement<Pair<Sequence,Double>, Double>( new Pair<Sequence, Double>( bsDs.getElementAt( i ), weights[i] ), scores[i] );
		}
		Arrays.sort( els );
		
		int al = (int)bsDs.getAlphabetContainer().getAlphabetAt( 0 ).length();
		
		double[][] stat = new double[bsDs.getElementLength()][bsDs.getElementLength()];
		double[][] mis = new double[bsDs.getElementLength()][bsDs.getElementLength()];
		double[][][][] counts = new double[bsDs.getElementLength()][bsDs.getElementLength()][al][al];
		
		double log2 = Math.log( 2 );
		
		double sum = 0.0;
		
		for(int i=els.length-1;i>=0;i--){
			Sequence seq = els[i].getElement().getFirstElement();
			double weight = els[i].getElement().getSecondElement();
			for(int j=0;j<seq.getLength();j++){
				for(int k=0;k<seq.getLength();k++){
					counts[j][k][seq.discreteVal( j )][seq.discreteVal( k )] += weight;
				}
			}
			sum += weight;
			
			for(int j=0;j<counts.length;j++){
				for(int k=0;k<counts[j].length;k++){
					double mi = 0.0;
					for(int a=0;a<counts[j][k].length;a++){
						for(int b=0;b<counts[j][k][a].length;b++){
							if(counts[j][k][a][b] > 0){
								mi += counts[j][k][a][b] * Math.log( sum*counts[j][k][a][b]/(counts[j][j][a][a]*counts[k][k][b][b]) )/log2;
							}
						}
					}
					if(mi > stat[j][k]){
						stat[j][k] = mi;
						mis[j][k] = mi/sum;
					}
				}
			}
			
		}
		
		return mis;
		
	}

	private static double[] getCounts( DataSet completeData, double[] ds ) {
		double[] counts = new double[(int)completeData.getAlphabetContainer().getAlphabetLengthAt( 0 )];
		for(int i=0;i<completeData.getNumberOfElements();i++){
			Sequence seq = completeData.getElementAt( i );
			Sequence ref = ((ReferenceSequenceAnnotation)seq.getSequenceAnnotationByTypeAndIdentifier( "reference", "reads" )).getReferenceSequence();
			for(int j=0;j<seq.getLength();j++){
				counts[seq.discreteVal( j )] += ds[i]*ref.continuousVal( j );
			}
		}
		return counts;
	}


	// creates a new homogeneous model	
	private static DifferentiableStatisticalModel getBgSF( AlphabetContainer con, int order, double ess, double length ) throws Exception {
		if( order >= 0 ) {
			return new HomogeneousMMDiffSM( con, order, ess, (int) Math.round( length ) );
		} else {
			return new UniformHomogeneousDiffSM( con, ess );
		}		
	}
	
	private static ArrayList<ComparableElement<double[], Double>> filter2( AbstractSingleMotifChIPper chipper, ComparableElement<double[], Double>[] pars, DataSet fg, double t, int length, Protocol protocol ) throws Exception {
		ArrayList<ComparableElement<double[], Double>> list = new ArrayList<ComparableElement<double[], Double>>(10);
		double[][] current;

		ArrayList<double[][]> profiles = new ArrayList<double[][]>();
				
		outerloop:
		for( int i = pars.length-1; i >= 0; i-- ) {
			chipper.setParameters( pars[i].getElement(), 2 );
			//current = ((MarkovModelDiffSM) chipper.getFunction(0)).getPWM();
			
			//Filter for Poly-N
			/*double max = Double.NEGATIVE_INFINITY;
			for(int j=0;j<current.length;j++){
				double m = ToolBox.max( current[j] );
				if(m > max){
					max = m;
				}
			}
			if(max < 0.4){
				continue outerloop;
			}/**/
			
			double[][] profile = new double[fg.getNumberOfElements()][];
			for(int k=0;k<profile.length;k++){
				profile[k] = chipper.getProfileOfScoresFor( 0, 0, fg.getElementAt( k ), 0, KindOfProfile.UNNORMALIZED_JOINT);
				//Normalisation.logSumNormalisation( profile[k] );
			}
			
			for(int j=0;j<profiles.size();j++){
				double corr = getCorrelation(profiles.get( j ),profile,length);
				
				
				if(corr > t){
					continue outerloop;
				}
			}
			
			profiles.add( profile );
			list.add( pars[i] );
			//out.writeln("added: "+getConsensus(chipper.getAlphabetContainer(), current));
		}
		
		if( list.size() == 0 ) {
			//add first
			int i = pars.length-1;
			double[][] profile = new double[fg.getNumberOfElements()][];
			for(int k=0;k<profile.length;k++){
				profile[k] = chipper.getProfileOfScoresFor( 0, 0, fg.getElementAt( k ), 0, KindOfProfile.UNNORMALIZED_JOINT);
				//Normalisation.logSumNormalisation( profile[k] );
			}
			profiles.add( profile );
			list.add( pars[i] );
			chipper.setParameters( pars[i].getElement(), 2 );
			//current = ((MarkovModelDiffSM) chipper.getFunction(0)).getPWM();
			//out.writeln("added: "+getConsensus(chipper.getAlphabetContainer(), current));
		}
		protocol.append("number of motifs: " + list.size()+"\n" );
		return list;
	}
	
	private static boolean[] postFilter( MutableMotifDiscoverer[] disc , int[] order, DataSet fg, double t, int length ) throws Exception {
		double[][] current;

		ArrayList<double[][]> profiles = new ArrayList<double[][]>();
		
		boolean[] use = new boolean[disc.length];
		
		outerloop:
		for( int i = 0; i <order.length; i++ ) {
		//	current = ((MarkovModelDiffSM) ((AbstractSingleMotifChIPper)disc[order[i]]).getFunction(0)).getPWM();
			double[][] profile = new double[fg.getNumberOfElements()][];
			for(int k=0;k<profile.length;k++){
				profile[k] = disc[order[i]].getProfileOfScoresFor( 0, 0, fg.getElementAt( k ), 0, KindOfProfile.UNNORMALIZED_JOINT);
			}
			
			for(int j=0;j<profiles.size();j++){
				double corr = getCorrelation(profiles.get( j ),profile,length);
				
				
				if(corr > t){
					continue outerloop;
				}
			}
			
			profiles.add( profile );
			use[i] = true;
		}
		return use;
	}
	
	private static double getCorrelation( double[][] ds, double[][] profile, int length ) throws Exception {
		
		double max = Double.NEGATIVE_INFINITY;
		
		for(int off=0;off<length;off++){
			double currCorr1 = 0;
			double currCorr2 = 0;
			for(int i=0;i<ds.length;i++){

				double p1 = ToolBox.pearsonCorrelation( ds[i], profile[i], 0, off );
				double p2 = ToolBox.pearsonCorrelation( ds[i], profile[i], off, 0 );
				
				currCorr1 += p1;
				currCorr2 += p2;
				
			}
			if(currCorr1 > max){
				max = currCorr1;
			}
			if(currCorr2 > max){
				max = currCorr2;
			}
		}
		return max/ds.length;
	}
	
	public static String getConsensus( AlphabetContainer con, double[][] pfm ) {
		String c = "";
		for( int s, m, p, l = 0; l < pfm.length; l++ ) {
			m = pfm[l][0] > pfm[l][1] ? 0 : 1;
			s = 1 - m;
			for( p = 2; p < pfm[l].length; p++ ) {
				if( pfm[l][m] < pfm[l][p] ) {
					s = m;
					m = p;
				} else if( pfm[l][s] < pfm[l][p] ){
					s = p;
				}
			}
			if( pfm[l][m] > 0.4 ) {
				if( pfm[l][m] - pfm[l][s] > 0.1 ) {
					c += con.getSymbol( l, m );
				} else {
					c += con.getSymbol( l, m ).toLowerCase();
				}
			} else {
				c += "N";
			}
		}
		return c;
	}
	
	public static ComparableElement<String, Double>[] getKmereSequenceStatistic( int numWanted, int k, DataSet data, double[] weights ) throws Exception {
		AlphabetContainer con = data.getAlphabetContainer();
		if( !con.isSimple() || !con.isDiscrete() ) {
			throw new WrongAlphabetException();
		}
		Hashtable<String, double[]> res = new Hashtable<String, double[]>();
		HashSet<String> used = new HashSet<String>();
		
		//run over all sequences
		int m, l, n;
		Sequence seq;
		String[] s = new String[2];
		for( n = 0; n < weights.length; n++ ) {
			seq = data.getElementAt( n );
			s[0] = seq.toString();
			s[1] = seq.reverseComplement().toString();
			m = seq.getLength()-k+1;
			
			
			//run over all k-mers
			used.clear();
			for( l = 0; l < m; l++ ) {
				String h0 = s[0].substring( l, l+k );
				String h1 = s[1].substring( s[0].length()-k-l, s[0].length()-l );
				h0 = h0.compareTo(h1) < 0 ? h0 : h1;
				if( !used.contains( h0 ) ) {
					used.add( h0 );
				}
				
			}
			
			Iterator<String> it = used.iterator();
			double[] h;
			while( it.hasNext() ) {
				s[0] = it.next();
				if( res.containsKey( s[0] ) ) {
					h = res.get( s[0] );
					h[0] += weights[n];
					h[1] += 1d-weights[n];
				} else {
					res.put( s[0], new double[]{weights[n],1d-weights[n]} );
				}
			}
		}
		
		
		double sumFg = ToolBox.sum( weights );
		
		ComparableElement<String, Double>[] array = new ComparableElement[res.size()];
		Iterator<Entry<String, double[]>> it = res.entrySet().iterator();
		Entry<String, double[]> e;
		for( int a = 0; a < array.length; a++ ) {
			e = it.next();
			double[] val = e.getValue();
			array[a] = new ComparableElement<String, Double>(e.getKey(), Math.log(val[0]+1)*(val[0]+1)/(val[1]+1) );
			
			//array[a] = new ComparableElement<String, Double>(e.getKey(), val[0]*corr + val[1]*wro + (sumFg-val[0])*wro + (sumBg-val[1])*corr );
		}
		Arrays.sort(array);
		
		if(numWanted > array.length){
			numWanted = array.length;
		}
		
		ComparableElement<String, Double>[] resArray = new ComparableElement[numWanted];
		Sequence[] prevs = new Sequence[numWanted];
		int j=resArray.length-1;
		outerloop:
		for(int i=array.length-1;i>=0;i--){
			String curr = array[i].getElement();
			
			Sequence currs = Sequence.create( con, curr);
			for(int a=resArray.length-1;a>j;a--){
				Sequence prev = prevs[a];
				if(getMinimumHammingDistance( currs, prev ) < 2){
					//System.out.println(currs+" == "+prev);
					continue outerloop;
				}
			}
			//System.out.println("added "+currs);
			resArray[j] = array[i];
			prevs[j] = currs;
			j--;
			if(j < 0){
				break;
			}
		}
		
		//avoid null in the array
		if( j >= 0 ) {
			ComparableElement<String, Double>[] help = new ComparableElement[numWanted-(j+1)];
			System.arraycopy(resArray, j+1, help, 0, help.length);
			resArray=help;
		}
		
		//resArray = new ComparableElement[1];
		//resArray[0] = new ComparableElement<String, Double>( "CACGTG", 0.0 );
		
		
		return resArray;
	}
	
	private static int getMinimumHammingDistance(Sequence curr, Sequence seq2) throws Exception {
		int min = Integer.MAX_VALUE;
		for(int i=0;i<=curr.getLength()/3;i++){
			Sequence sub1 = curr.getSubSequence( i );
			Sequence subRc = curr.reverseComplement().getSubSequence( i );
			Sequence sub2 = seq2.getSubSequence( 0, seq2.getLength()-i );
			int d1 = sub1.getHammingDistance( sub2 );
			int d2 = subRc.getHammingDistance( sub2 );
			if(d1 < min){
				min = d1;
			}
			if(d2 < min){
				min = d2;
			}
		}
		for(int i=1;i<=curr.getLength()/3;i++){
			Sequence sub1 = curr.getSubSequence(0, curr.getLength()- i );
			Sequence subRc = curr.reverseComplement().getSubSequence( 0, curr.getLength()- i );
			Sequence sub2 = seq2.getSubSequence( i );
			//	System.out.println(sub1.getLength()+" "+subRc.getLength()+" "+sub2.getLength());
			int d1 = sub1.getHammingDistance( sub2 );
			int d2 = subRc.getHammingDistance( sub2 );
			if(d1 < min){
				min = d1;
			}
			if(d2 < min){
				min = d2;
			}
		}
		return min;
	}


	@Override
	public String getToolName() {
		return "Methyl SlimDimont";
	}


	@Override
	public String getToolVersion() {
		return "0.1";
	}


	@Override
	public String getShortName() {
		return "slimdimont";
	}


	@Override
	public String getDescription() {
		return "";
	}


	@Override
	public String getHelpText() {
		return "**"+getToolName()+"** is a tool for de-novo motif discovery from DNA sequences including extended, e.g., methylation-aware alphabets.\n" + 
				"\n" + 
				"Input sequences must be supplied in an annotated FastA format as generated by the Data Extractor tool.\n" + 
				"Input sequences may also obtained from other sources. In this case, the annotation of each sequence needs to provide a value that reflects the confidence that this sequence is bound by the factor of interest.\n" + 
				"Such confidences may be peak statistics (e.g., number of fragments under a peak) for ChIP data or signal intensities for PBM data. In addition, you need to provide an anchor position within the sequence. \n" + 
				"In case of ChIP data, this anchor position could for instance be the peak summit.\n" + 
				"An annotated FastA file for ChIP-seq data comprising sequences of length 100 centered around the peak summit could look like::\n" + 
				"	\n" + 
				"	> peak: 50; signal: 515\n" + 
				"	ggccatgtgtatttttttaaatttccac...\n" + 
				"	> peak: 50; signal: 199\n" + 
				"	GGTCCCCTGGGAGGATGGGGACGTGCTG...\n" + 
				"	...\n" + 
				"\n" + 
				"where the anchor point is given as 50 for the first two sequences, and the confidence amounts to 515 and 199, respectively.\n" + 
				"The FastA comment may contain additional annotations of the format ``key1 : value1; key2: value2;...``.\n" + 
				"\n" + 
				"Accordingly, you would need to set the parameter \"Position tag\" to ``peak`` and the parameter \"Value tag\" to ``signal`` for the input file (default values).\n"
				+ "The parameter Alphabet specifies the symbols of the (extended) alphabet and their complementary symbols. Default is standard DNA alphabet.\n" + 
				"\n" + 
				"For the standard deviation of the position prior, the initial motif length and the number of pre-optimization runs, we provide default values that worked well in our studies on ChIP and PBM data. \n" + 
				"However, you may want adjust these parameters to meet your prior information.\n" + 
				"\n" + 
				"The parameter \"Markov order of the motif model\" sets the order of the inhomogeneous Markov model used for modeling the motif. If this parameter is set to ``0``, you obtain a position weight matrix (PWM) model. \n" + 
				"If it is set to ``1``, you obtain a weight array matrix (WAM) model. You can set the order of the motif model to at most ``3``.\n" + 
				"\n" + 
				"The parameter \"Markov order of the background model\" sets the order of the homogeneous Markov model used for modeling positions not covered by a motif. \n" + 
				"If this parameter is set to ``-1``, you obtain a uniform distribution, which worked well for ChIP data. For PBM data, orders of up to ``4`` resulted in an increased prediction performance in our case studies. The maximum allowed value is ``5``.\n" + 
				"\n" + 
				"The parameter \"Weighting factor\" defines the proportion of sequences that you expect to be bound by the targeted factor with high confidence. For ChIP data, the default value of ``0.2`` typically works well. \n" + 
				"For PBM data, containing a large number of unspecific probes, this parameter should be set to a lower value, e.g. ``0.01``.\n" + 
				"\n" + 
				"The \"Equivalent sample size\" reflects the strength of the influence of the prior on the model parameters, where higher values smooth out the parameters to a greater extent.\n" + 
				"\n" + 
				"The parameter \"Delete BSs from profile\" defines if BSs of already discovered motifs should be deleted, i.e., \"blanked out\", from the sequence before searching for futher motifs.\n" + 
				"\n" + 
				"You can also install this web-application within your local Galaxy server. Instructions can be found at the Dimont_ page of Jstacs. \n" + 
				"There you can also download a command line version of Dimont.\n" + 
				"\n" + 
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
		return new String[] {
				"@article{grau13a-general,\n" + 
				"	Author = {Grau, Jan and Posch, Stefan and Grosse, Ivo and Keilwagen, Jens},\n" + 
				"	Journal = {Nucleic Acids Research},\n" + 
				"	Number = {21},\n" + 
				"	Pages = {e197},\n" + 
				"	Title = {A general approach for discriminative de novo motif discovery from high-throughput data},\n" + 
				"	Volume = {41},\n" + 
				"	Year = {2013}}",
				"@article{keilwagen15varying,\n" + 
				"	Author = {Keilwagen, Jens and Grau, Jan},\n" + 
				"	Journal = {Nucleic Acids Research},\n" + 
				"	Month = {10},\n" + 
				"	Number = {18},\n" + 
				"	Pages = {e119-e119},\n" + 
				"	Title = {{Varying levels of complexity in transcription factor binding motifs}},\n" + 
				"	Volume = {43},\n" + 
				"	Year = {2015}}"
				
		};
	}
}