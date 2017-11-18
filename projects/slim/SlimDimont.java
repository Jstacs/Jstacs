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

package projects.slim;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;

import javax.imageio.ImageIO;
import javax.naming.OperationNotSupportedException;

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
import de.jstacs.algorithms.optimization.termination.TimeCondition;
import de.jstacs.classifiers.differentiableSequenceScoreBased.AbstractMultiThreadedOptimizableFunction;
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
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.ArbitraryFloatSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.SparseSequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.io.FileManager;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile;
import de.jstacs.motifDiscovery.MutableMotifDiscoverer;
import de.jstacs.motifDiscovery.MutableMotifDiscovererToolbox;
import de.jstacs.motifDiscovery.MutableMotifDiscovererToolbox.InitMethodForDiffSM;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder;
import de.jstacs.parameters.ParameterSetTagger;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ImageResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.MarkovModelDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.UniformHomogeneousDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.localMixture.LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder;
import de.jstacs.sequenceScores.statisticalModels.differentiable.localMixture.LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder.PriorType;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.SeqLogoPlotter;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.ToolBox.TiedRanks;
//import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.MarkovModelDiffSM;
import projects.dimont.AbstractSingleMotifChIPper;
import projects.dimont.HeuristicOneDataSetLogGenDisMixFunction;
import projects.dimont.Interpolation;
import projects.dimont.ThresholdedStrandChIPper;


public class SlimDimont { 
	
	private static final double ALPHA = 1E-3;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {

		//parse parameters
		ParameterSetTagger cParams = new ParameterSetTagger( SlimDimontParameterSet.PREFIX, new SlimDimontParameterSet() );
		cParams.fillParameters( "=", args );
		System.out.println( "parameters:" );
		System.out.println( cParams );
		System.out.println("_________________________________");
		if( !cParams.hasDefaultOrIsSet() ) {
			System.out.println( "Some of the required parameters are not specified." );
			System.exit( 1 );
		}

		
		
		String home = cParams.getValueFromTag( SlimDimontParameterSet.HOME, String.class );
		String fgData = home +File.separator+ cParams.getValueFromTag( SlimDimontParameterSet.DATA, String.class );
		String bg = cParams.getValueFromTag( SlimDimontParameterSet.BACKGROUND, String.class );
		String infix = cParams.getValueFromTag( SlimDimontParameterSet.INFIX, String.class );
		int motifLength = cParams.getValueFromTag( SlimDimontParameterSet.LENGTH, Integer.class );
		int restarts = cParams.getValueFromTag( SlimDimontParameterSet.STARTS, Integer.class );
		int fgOrder = cParams.getValueFromTag( SlimDimontParameterSet.MOTIF_ORDER, Integer.class );
		int bgOrder = cParams.getValueFromTag( SlimDimontParameterSet.BG_ORDER, Integer.class );
		String position = cParams.getValueFromTag( SlimDimontParameterSet.POSITION_TAG, String.class );
		String value = cParams.getValueFromTag( SlimDimontParameterSet.VALUE_TAG, String.class );
		double sd = (Double)cParams.getValueFromTag( SlimDimontParameterSet.SD, Double.class );
		String weightingFactor = cParams.getValueFromTag( SlimDimontParameterSet.WEIGHTING_FACTOR, String.class );
		double ess = cParams.getValueFromTag( SlimDimontParameterSet.ESS, Double.class );
		boolean delete = cParams.getValueFromTag( SlimDimontParameterSet.DELETE, Boolean.class );
		boolean modify = cParams.getValueFromTag( SlimDimontParameterSet.MODIFY, Boolean.class );
		
		int threads;
		if( cParams.isSet(SlimDimontParameterSet.THREADS) ) {
			threads = cParams.getValueFromTag( SlimDimontParameterSet.THREADS, Integer.class );
		} else {
			threads = AbstractMultiThreadedOptimizableFunction.getNumberOfAvailableProcessors();
		}
		
		SequenceAnnotationParser parser = new SplitSequenceAnnotationParser(":", ";");
		
		DataSet data = SparseSequence.getDataSet(DNAAlphabetContainer.SINGLETON, fgData, parser);
		if( data.getNumberOfElements() < 10 ) {
			System.out.println("WARNING: SlimDimont is not made for small data sets, i.e., a small number of sequences. Hence, it might return useless results.");
		}
		DataSet bgData = bg==null ? null : SparseSequence.getDataSet(DNAAlphabetContainer.SINGLETON, new SparseStringExtractor(home +File.separator+ bg,'>'));
		
		/*if(fgData.matches( ".*_part_[0-9].fa" )){
			String temp = fgData.replaceAll( "_part_[0-9].fa", ".fa" );
			System.out.println("Reading "+temp);
			DataSet all = SparseSequence.getDataSet(DNAAlphabetContainer.SINGLETON, temp, parser);
			System.out.println("test: "+data.getNumberOfElements());
			System.out.println("all: "+all.getNumberOfElements());
			data = DataSet.diff( all, data );
			System.out.println("diff: "+data.getNumberOfElements());
		}*/
		
		Result[][] res = run( data, bgData, motifLength,restarts,fgOrder,bgOrder,position,value,weightingFactor,ess,delete,threads, SafeOutputStream.getSafeOutputStream( System.out ),sd,modify);
		
		for(int i=0;i<res.length;i++){
			StorableResult sr = (StorableResult) res[i][0];
			FileManager.writeFile( new File(home+File.separator+infix+"-model-"+(i+1)+".xml"), sr.getValue() );
			
			System.out.println("+++++++++++++++++++++++++++++++++++++++++++++\nMotif model "+(i+1)+":");
			System.out.println( ( (GenDisMixClassifier)sr.getResultInstance() ).getDifferentiableSequenceScore( 0 ));
			
			ListResult lr = (ListResult) res[i][1];
			FileManager.writeFile( new File(home+File.separator+infix+"-predictions-"+(i+1)+".txt"), lr.toString() );
			
			if(res[i].length > 2){
				ImageResult ir = (ImageResult) res[i][2];
				ImageIO.write( ir.getValue() , "png", new File(home+File.separator+infix+"-logo-"+(i+1)+".png") );
				ir = (ImageResult) res[i][3];
				ImageIO.write( ir.getValue() , "png", new File(home+File.separator+infix+"-logo-rc-"+(i+1)+".png") );
				ir = (ImageResult) res[i][4];
				ImageIO.write( ir.getValue() , "png", new File(home+File.separator+infix+"-dependencylogo-"+(i+1)+".png") );
			}
			
		}
		
	}
	
	public static Result[][] run(DataSet fgData, DataSet bgData, int motifLength, int restarts, int fgOrder, int bgOrder, String position,
			String value, String weightingFactor, double ess, boolean delete, int threads, SafeOutputStream out, double sd,
			boolean modify) throws Exception {
	
		double filterThreshold = 0.3;
		double filterThresholdEnd = 0.3;
		
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		
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
		
		//create weights
		double wf;
		if( weightingFactor.endsWith("sd") ) {
			double h = Double.parseDouble( weightingFactor.substring(0,weightingFactor.length()-2) );
			double meanRaw = ToolBox.sum(raw) / raw.length;
			double sdRaw = 0;
			for( int i = 0; i < raw.length; i++ ) {
				sdRaw += (raw[i]-meanRaw) * (raw[i]-meanRaw);
			}
			sdRaw = Math.sqrt( sdRaw/raw.length );
			h = meanRaw + h*sdRaw;
			double anz = 0;
			for( int i = 0; i < raw.length; i++ ) {
				if( raw[i] >= h ) {
					anz++;
				}
			}
			anz=Math.max(50,anz);
			wf = anz/raw.length;
		} else {
			wf = Double.parseDouble( weightingFactor );
		}

		double[] w = Interpolation.getWeight( data[0], raw, wf, Interpolation.RANK_LOG );
		if( bgData == null ) {
			weights[0] = w;
		} else {
			weights[0] = new double[w.length + bgData.getNumberOfElements()];
			System.arraycopy(w, 0, weights[0], 0, w.length);
		}
		weights[1] = Interpolation.getBgWeight( weights[0] );
		
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
		if(fgOrder >= 0){
				motif = new MarkovModelDiffSM( con, motifLength, ess, true, fgOrder, null );
		}else{
			motif = new LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder( con, motifLength, 1, -fgOrder, ess, 0.9, PriorType.BDeu );
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
		
		//preoptimization
		AbstractTerminationCondition stop2 = new CombinedCondition( 3,
				new MultipleIterationsCondition( 5, new SmallDifferenceOfFunctionEvaluationsCondition(eps) ),
				new IterationCondition(50),
				new TimeCondition(900)
		);
		ComparableElement<double[], Double>[] preOpt = new ComparableElement[restarts];
		for( int r = 0; r < restarts; r++ ) {
			data[0] = smallData;
			weights = smallWeight;
			objective.setDataAndWeights( data, weights );
			objective.resetHeuristics();
			
			out.writeln("-----------------------------------------\npre-optimization " + r);
			
			start.reset();
			
			p = sortedPars[sortedPars.length-1-r].getElement();			
			objective.setParams(p);
			System.out.println(score[0]);
			long time = System.currentTimeMillis();
			Optimizer.optimize( algo, neg, p, stop2, eps, start, null );
			
			data[0] = completeData;
			weights = completeWeight;
			objective.setDataAndWeights( data, weights );
			preOpt[r] = new ComparableElement<double[], Double>( p, objective.evaluateFunction(p) );
			System.out.println("time: "+(System.currentTimeMillis()-time));
			
			//out.writeln("consensus: "+getConsensus( DNAAlphabetContainer.SINGLETON, (((MarkovModelDiffSM)((AbstractSingleMotifChIPper) score[0]).getFunction( 0 ) ).getPWM())));
			out.writeln("model: "+((AbstractSingleMotifChIPper) score[0]).getFunction( 0 ));
			out.writeln(  "score: "+preOpt[r].getWeight() );
			//System.out.println(score[0]);
			((AbstractSingleMotifChIPper)score[0]).resetPositions();
		}
		out.writeln("-----------------------------------------");
		
		//filter out redundant motifs
		Arrays.sort( preOpt );
		ArrayList<ComparableElement<double[], Double>> list = filter2( (AbstractSingleMotifChIPper) score[0], preOpt, smallData, filterThreshold ,motifLength, out );

		//final optimization
		MutableMotifDiscoverer[] best = new MutableMotifDiscoverer[list.size()];		
		Storable[] storables = new Storable[best.length];
		double[] opts = new double[best.length];
		Pair<double[][][],int[][]>[] pairs = new Pair[best.length];		
		for( int m = 0; m < best.length; m++ ) {
			out.writeln("+++++++++++++++++++++++++++++++++++++++++++++++++++" +
					"\n\nfinalOpt " + m + " -----------------------------------------");
			
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

			//heuristic for motif length
			double[] sds = new double[1];
			heuristic((MutableMotifDiscoverer) score[0], completeData, completeWeight, objective,mean,sds, out, modify);
			
			
			//final optimization
			((AbstractSingleMotifChIPper)score[0]).reset();
			
			double newSd = Math.sqrt( sds[0]*sd);
			if(Double.isNaN( newSd ) || Double.isInfinite( newSd ) || newSd <= 0){
				out.writeln("Did not adjust sd to "+newSd+" using "+sds[0]+" and "+sd);
				newSd = sd;
			}
			
			data[0] = annotate( annotated, weights, mean, newSd, allowed, fgData.getNumberOfElements() );
			objective.setDataAndWeights( data, weights );
			((AbstractSingleMotifChIPper) score[0]).resetPositions();
			objective.reset( score ); // reset heuristics: toBeUsed..., ...
			
			p = objective.getParameters(KindOfParameter.LAST);
			objective.setParams(p);
			Optimizer.optimize( algo, neg, p, stop, eps*1E-1, start, SafeOutputStream.getSafeOutputStream( null ) );
			
			out.writeln("model: "+((AbstractSingleMotifChIPper) score[0]).getFunction( 0 ));
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
		}
		
		//rank and filter final motifs
		int[] o = ToolBox.rank( opts, TiedRanks.IN_ORDER );
		int[] index = new int[o.length];
		for( int i = 0; i<o.length; i++ ) {
			index[o[i]] = i;
		}
		boolean[] use = postFilter( best, index, smallData, filterThresholdEnd, motifLength );
		
		LinkedList<Result[]> results = new LinkedList();
		
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
						int height = SeqLogoPlotter.getHeight( 750, pwm );
						result.add(new ImageResult( "Motif "+(n+1), "Sequence logo of motif "+(n+1), SeqLogoPlotter.plotLogoToBufferedImage( height, pwm ) ));
						result.add(new ImageResult( "Motif "+(n+1)+" (rc)", "Sequence logo of the reverse complement of motif "+(n+1), SeqLogoPlotter.plotLogoToBufferedImage( height, PFMComparator.getReverseComplement( DNAAlphabet.SINGLETON, pwm ) ) ));
						result.add( new ImageResult( "Dependency logo "+(n+1), "Dependency logo of motif "+(n+1), SeqLogoPlotter.plotDefaultDependencyLogoToBufferedImage( new DataSet("",bs), bsWeights.toArray(), 600 ) ) );
					}catch(Exception e){
						e.printStackTrace();
					}
				}
				results.add( result.toArray( new Result[0] ) );
				n++;
			}
		}
		initObjective.stopThreads();
		objective.stopThreads();
		
		return results.toArray( new Result[0][] );
		
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
	
	private static boolean heuristic( MutableMotifDiscoverer best, DataSet completeData, double[][] completeWeight, LogGenDisMixFunction objective, double[] mean, double[] sds, SafeOutputStream out, boolean modify ) throws Exception {
		SignificantMotifOccurrencesFinder smof = new SignificantMotifOccurrencesFinder(best,completeData,completeWeight[1],ALPHA);//XXX Jan
		boolean modified = false;
		double log2 = Math.log( 2 );
		for(int im=0;im<best.getNumberOfMotifs();im++){
			int add = 5;
			if(completeData.getAverageElementLength() - best.getMotifLength( im ) < 20){
				add = (int)((completeData.getAverageElementLength() - best.getMotifLength( im ))/4.0);
			}
			//add = 0;//TODO
			System.out.println("add: "+add);
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
			System.out.println("kl: "+Arrays.toString( kl ));
			DataSet bsDs = new DataSet( "", bs );
			//double[][] mi2 = TwoPointEvaluater.getMIInBits( bsDs, bsWeights.toArray() );
			/*for(int i=0;i<mi2.length;i++){
				System.out.println(i+" "+Arrays.toString( mi2[i] ));
			}*/
			double[][] mi = getMaxMI(bsDs, bsWeights.toArray(), bsScores.toArray());
			for(int i=0;i<mi.length;i++){
				System.out.println(i+" "+Arrays.toString( mi[i] ));
			}
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
			System.out.println("kl: "+Arrays.toString( kl ));
			System.out.println(bsDs.getNumberOfElements()+": "+bsDs.getAverageElementLength()+", "+bsDs.getElementLength());
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
			
			out.writeln( "left: "+ left + ", right: " + right );
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
			
			
			out.writeln( "left: "+ left + ", right: " + right );
			if( left+add == kl.length || right-add == -kl.length ) {
				out.writeln( "tried to remove the complete motif: no modifications");
			} else {
				double normOld = ((DifferentiableStatisticalModel) best).getLogNormalizationConstant();
				if( left != 0 || right!= 0 ) {
					if( best.modifyMotif( im, left, right ) ){
						double w = normOld - ((DifferentiableStatisticalModel) best).getLogNormalizationConstant();
						//System.out.println(w);
						objective.addTermToClassParameter( 0, w );
					}
					out.writeln("modified motif");
					modified = true;
				} else {
					out.writeln( "no modifications for the motif");
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
	
	private static ArrayList<ComparableElement<double[], Double>> filter2( AbstractSingleMotifChIPper chipper, ComparableElement<double[], Double>[] pars, DataSet fg, double t, int length, SafeOutputStream out ) throws Exception {
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
		out.writeln("number of motifs: " + list.size() );
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
			
			Sequence currs = Sequence.create( DNAAlphabetContainer.SINGLETON, curr);
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
}