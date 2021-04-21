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

package projects.dimont;

import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;

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
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.ArbitraryFloatSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.SparseSequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.io.FileManager;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile;
import de.jstacs.motifDiscovery.MutableMotifDiscoverer;
import de.jstacs.motifDiscovery.MutableMotifDiscovererToolbox;
import de.jstacs.motifDiscovery.MutableMotifDiscovererToolbox.InitMethodForDiffSM;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;
import de.jstacs.results.StorableResult;
import de.jstacs.results.TextResult;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.MarkovModelDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.UniformHomogeneousDiffSM;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.SeqLogoPlotter.SeqLogoPlotGenerator;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.ToolBox.TiedRanks;

public class DimontTool implements JstacsTool {

	private static final double ALPHA = 1E-3;
	
	
	public DimontTool() {
	}

	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> parameters = new LinkedList<Parameter>();
		
		
		parameters.add(new FileParameter("Input file", "The file name of the file containing the input sequences in annotated FastA format (see readme)", "fasta,fa,fas", true));
		try{
		parameters.add( new SimpleParameter( DataType.STRING, "Position tag", "The tag for the position information in the FastA-annotation of the input file", true,"peak" ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Value tag", "The tag for the value information in the FastA-annotation of the input file", true, "signal" ) );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Standard deviation", "The standard deviation of the position distribution centered at the position specified by the position tag", true, new NumberValidator<Double>( 1.0, 1E4 ), 75.0 ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Weighting factor", "The value for weighting the data; either a value between 0 and 1, or a description relative to the standard deviation (e.g. +4sd)", true, "" + 0.2 ) );
		
		
		parameters.add( new SimpleParameter( DataType.INT, "Starts", "The number of pre-optimization runs.", true, new NumberValidator<Integer>(1,100), 20 ) );
		
		
		parameters.add( new SimpleParameter( DataType.INT, "Initial motif width", "The motif width that is used initially, may be adjusted during optimization.", true, new NumberValidator<Integer>(1,50), 15 ) );
		
		parameters.add( new SimpleParameter( DataType.INT, "Markov order of motif model", "The Markov order of the model for the motif.", true, new NumberValidator<Integer>(0,3), 0 ) );
		parameters.add( new SimpleParameter( DataType.INT, "Markov order of background model", "The Markov order of the model for the background sequence and the background sequence, -1 defines uniform distribution.", true, new NumberValidator<Integer>(-1,5), -1 ) );
		
		
		
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Equivalent sample size", "Reflects the strength of the prior on the model parameters.", true, new NumberValidator<Double>(0d, Double.POSITIVE_INFINITY), 4d ) );
		
		parameters.add( new SimpleParameter( DataType.BOOLEAN, "Delete BSs from profile", "A switch for deleting binding site positions of discovered motifs from the profile before searching for futher motifs.", true, true ) );
		
		}catch(Exception e){
			throw new RuntimeException();
		}
		
		//parameters.add( new SimpleParameter( DataType.INT, "Compute threads", "The number of threads that are use to evaluate the objective function and its gradient.", false, new NumberValidator<Integer>(1,128) ) );

		return new ToolParameterSet(getShortName(),parameters.toArray(new Parameter[0]));
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		
		progress.setLast( 1 );
		progress.setCurrent( 0.0 );
		
		
		DataSet fgData = SparseSequence.getDataSet( DNAAlphabetContainer.SINGLETON, new SparseStringExtractor( 
				new StringReader(
						((FileParameter)parameters.getParameterAt(0)).getFileContents().getContent())
				, '>', "", new SplitSequenceAnnotationParser(":",";") ) );
		
		
		String position = parameters.getParameterAt(1).getValue().toString();
		String value = parameters.getParameterAt(2).getValue().toString();
		double sd = (Double) parameters.getParameterAt(3).getValue();
		String weightingFactor = parameters.getParameterAt(4).getValue().toString();
		int restarts = (Integer) parameters.getParameterAt(5).getValue();
		int motifLength = (Integer) parameters.getParameterAt(6).getValue();
		int fgOrder = (Integer) parameters.getParameterAt(7).getValue();
		int bgOrder = (Integer) parameters.getParameterAt(8).getValue();
		double ess = (Double) parameters.getParameterAt(9).getValue();
		boolean delete = (Boolean) parameters.getParameterAt(10).getValue();
		
		
		double filterThreshold = 0.3;
		double filterThresholdEnd = 0.3;
		
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		
		boolean free = false;
		
		DataSet[] data = { fgData };
		

		Sequence[] annotated = new Sequence[data[0].getNumberOfElements()];
		double[][] weights = new double[2][data[0].getNumberOfElements()];
		double[] raw = weights[0].clone();
		
		//read annotation
		double[] mean = new double[weights[0].length];
		Arrays.fill( mean, Double.NaN );
		for( int j = 0; j < weights[0].length; j++ ) {
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

		weights[0] = Interpolation.getWeight( data[0], raw, wf, Interpolation.RANK_LOG );
		weights[1] = Interpolation.getBgWeight( weights[0] );
		
		boolean[][] allowed = new boolean[annotated.length][];//TODO Jan? BitSets?
		for(int i=0;i<annotated.length;i++){
			allowed[i] = new boolean[annotated[i].getLength()];
			Arrays.fill( allowed[i], true );
		}
		
		//create initial annotation
		data[0] = annotate(annotated,weights,mean,sd,allowed);
		
		//complete vs small fg data/weights
		DataSet completeData = data[0];
		double[][] completeWeight = { weights[0], weights[1] };
		
		//create models
		DifferentiableStatisticalModel motif = //new MarkovRandomFieldDiffSM(con, 8, ess, "m2sx");
				new MarkovModelDiffSM( con, motifLength, ess, true, fgOrder, null );		
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
				((AbstractSingleMotifChIPper)score[0]).initializeMotif( 0, new DataSet("", Sequence.create(con, array[z].getElement()) ), new double[]{h} );
				p = objective.getParameters(KindOfParameter.PLUGIN);
			
				initObjective.reset(score);
				initObjective.resetHeuristics();
				double val = initObjective.evaluateFunction(p);
				
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
		AbstractTerminationCondition stop2 = new CombinedCondition( 2,
				new MultipleIterationsCondition( 5, new SmallDifferenceOfFunctionEvaluationsCondition(eps) ),
				new IterationCondition(25) 
		);
		ComparableElement<double[], Double>[] preOpt = new ComparableElement[restarts];
		for( int r = 0; r < restarts; r++ ) {
			
			progress.setCurrent(0.4/(double)restarts*(r));
			
			data[0] = smallData;
			weights = smallWeight;
			objective.setDataAndWeights( data, weights );
			objective.resetHeuristics();
			
			protocol.append("-----------------------------------------\npre-optimization " + r+"\n");
			
			start.reset();
			
			p = sortedPars[sortedPars.length-1-r].getElement();			
			objective.setParams(p);
			Optimizer.optimize( algo, neg, p, stop, eps*1E-1, start, null );
			
			data[0] = completeData;
			weights = completeWeight;
			objective.setDataAndWeights( data, weights );
			preOpt[r] = new ComparableElement<double[], Double>( p, objective.evaluateFunction(p) );
			
			protocol.append("consensus: "+getConsensus( DNAAlphabetContainer.SINGLETON, (((MarkovModelDiffSM)((AbstractSingleMotifChIPper) score[0]).getFunction( 0 ) ).getPWM()))+"\n");
			protocol.append(  "score: "+preOpt[r].getWeight()+"\n" );
			//System.out.println(score[0]);
			((AbstractSingleMotifChIPper)score[0]).resetPositions();
		}
		protocol.append("-----------------------------------------\n");
		progress.setCurrent(0.4);
		
		//filter out redundant motifs
		Arrays.sort( preOpt );
		ArrayList<ComparableElement<double[], Double>> list = filter2( (AbstractSingleMotifChIPper) score[0], preOpt, smallData, filterThreshold ,motifLength, protocol );

		//final optimization
		MutableMotifDiscoverer[] best = new MutableMotifDiscoverer[list.size()];		
		Storable[] storables = new Storable[best.length];
		double[] opts = new double[best.length];
		Pair<double[][][],int[][]>[] pairs = new Pair[best.length];		
		for( int m = 0; m < best.length; m++ ) {
			progress.setCurrent(0.4 + 0.6*( m*1.0/((double)best.length+1) ));
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
			data[0] = annotate( annotated, weights, mean, sd, allowed );
			objective.setDataAndWeights( data, weights );
			Optimizer.optimize( algo, neg, p, stop, eps*1E-1, start, SafeOutputStream.getSafeOutputStream( null ) );

			//heuristic for motif length
			double[] sds = new double[1];
			heuristic((MutableMotifDiscoverer) score[0], completeData, completeWeight, objective,mean,sds, protocol);
			
			
			//final optimization
			((AbstractSingleMotifChIPper)score[0]).reset();
			
			double newSd = Math.sqrt( sds[0]*sd);
			if(Double.isNaN( newSd ) || Double.isInfinite( newSd ) || newSd <= 0){
				protocol.append("Did not adjust sd to "+newSd+" using "+sds[0]+" and "+sd+"\n");
				newSd = sd;
			}
			
			data[0] = annotate( annotated, weights, mean, newSd, allowed );
			objective.setDataAndWeights( data, weights );
			((AbstractSingleMotifChIPper) score[0]).resetPositions();
			objective.reset( score ); // reset heuristics: toBeUsed..., ...
			
			p = objective.getParameters(KindOfParameter.LAST);
			objective.setParams(p);
			Optimizer.optimize( algo, neg, p, stop, eps*1E-1, start, SafeOutputStream.getSafeOutputStream( null ) );
			
			protocol.append("consensus: "+getConsensus( DNAAlphabetContainer.SINGLETON, (((MarkovModelDiffSM)((AbstractSingleMotifChIPper) score[0]).getFunction( 0 ) ).getPWM()))+"\n");
			
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
		
		LinkedList<ResultSetResult> results = new LinkedList<ResultSetResult>();
		
		//write final xmls and sequence logos
		for(int m=0,n=0;m<best.length;m++){
			if(use[m]){
				LinkedList<Result> result = new LinkedList<Result>();
				
				//result.add(new StorableResult( "Dimont "+(n+1), "The Dimont classifier", storables[index[m]] ) );//TODO
				result.add(new TextResult("Dimont "+(n+1), "The Dimont classifier", new FileRepresentation("", storables[index[m]].toXML().toString()), "xml", "Dimont", null, true));
				
				//result.add(getListResult(fgData, completeWeight[0],pairs[index[m]], ((ThresholdedStrandChIPper)((GenDisMixClassifier)storables[index[m]]).getDifferentiableSequenceScore( 0 )).getMotifLength( 0 ), n ));//TODO
				result.add(getTextResult(fgData, completeWeight[0],pairs[index[m]], ((ThresholdedStrandChIPper)((GenDisMixClassifier)storables[index[m]]).getDifferentiableSequenceScore( 0 )).getMotifLength( 0 ), n ));//TODO
				
				pwm = pairs[index[m]].getFirstElement()[0];//TODO
				
				/*StringBuffer temp = new StringBuffer();
				XMLParser.appendObjectWithTags( temp, pwm, "pwm" );
				
				result.add( new CategoricalResult( "pwm", "pwm", temp.toString() ) );
				*/
				
				if(!Double.isNaN( pwm[0][0] )){
					try{
						//int height = SeqLogoPlotter.getHeight( 750, pwm );
						
						result.add(new PlotGeneratorResult( "Sequence logo "+(n+1), "Sequence logo of motif "+(n+1), 
								new SeqLogoPlotGenerator(pwm, 1000), true));
						
						result.add(new PlotGeneratorResult( "Sequence logo "+(n+1)+" (rc)", "Sequence logo of the reverse complement of motif "+(n+1), 
								new SeqLogoPlotGenerator(PFMComparator.getReverseComplement( DNAAlphabet.SINGLETON, pwm ), 1000), true));
						
											}catch(Exception e){
						
					}catch(InternalError err){
						
					}
				}
				
				if(fgOrder == 0) {
					double[][] modelPwm = ((MarkovModelDiffSM)((ThresholdedStrandChIPper)((GenDisMixClassifier)storables[index[m]]).getDifferentiableSequenceScore(0)).getMotifModel()).getPWM();

					StringBuffer sb = new StringBuffer();
					sb.append(">Motif"+(n+1)+"\n");
					for(int i=0;i<modelPwm.length;i++) {
						for(int j=0;j<modelPwm[i].length;j++) {
							if(j>0) {
								sb.append("\t");
							}
							sb.append(modelPwm[i][j]);
						}
						sb.append("\n");
					}

					TextResult trPwm = new TextResult("Model PWM", "The model PWM in HOCOMOCO format", new FileParameter.FileRepresentation("", sb.toString()), "pwm", this.getToolName(), null, true);

					result.add(trPwm);
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

	@Override
	public String getToolName() {
		return "Dimont";
	}

	@Override
	public String getToolVersion() {
		return "1.2";
	}

	@Override
	public String getShortName() {
		return "dimont";
	}

	@Override
	public String getDescription() {
		return "a universal tool for de-novo motif discovery";
	}

	@Override
	public String getHelpText() {
			try {
				return FileManager.readInputStream( DimontTool.class.getClassLoader().getResourceAsStream( "projects/dimont/help.txt" ) ).toString();
			} catch ( IOException e ) {
				e.printStackTrace();
				return "";
			}
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[]{
				new ResultEntry(TextResult.class, "xml", "Dimont 1"),
				new ResultEntry(TextResult.class, "tsv", "Predictions for motif 1")
		};
	}

	
	public static TextResult getTextResult(DataSet data, double[] weights, Pair<double[][][], int[][]> pair, int motifLength, int motifIndex) throws OperationNotSupportedException {
		SplitSequenceAnnotationParser pars = new SplitSequenceAnnotationParser( ":", ";" );
		
		StringBuffer sb = new StringBuffer();
		
		sb.append("Sequence index\tPosition\tStrand\tp-value\tBinding site\tAdjusted binding site\tSequence annotation\n");
		
		int[][] pos = pair.getSecondElement();
		double[][] pvals = pair.getFirstElement()[1];
		for(int i=0;i<pos.length;i++){
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
				
				sb.append((i+1)+"\t"+(curr+1)+"\t"+(rc ? "-" : "+")+"\t"+pvals[i][j]+"\t"+sub.toString()+"\t"+sub2.toString()+"\t"+pars.parseAnnotationToComment( ' ', data.getElementAt( i ).getAnnotation() ).substring( 1 )+"\n");
				
			}
		}
		
		FileParameter.FileRepresentation file = new FileRepresentation("",sb.toString());
		
		TextResult tr = new TextResult("Predictions for motif "+(motifIndex+1), "", file, "tsv", "Dimont", null, true);
		
		return tr;
	}
	
	
	

	public static ListResult getListResult( DataSet data, double[] weights, Pair<double[][][], int[][]> pair, int motifLength, int motifIndex ) throws Exception {
		
		SplitSequenceAnnotationParser pars = new SplitSequenceAnnotationParser( ":", ";" );
		
		LinkedList<ResultSet> set = new LinkedList<ResultSet>();
		int[][] pos = pair.getSecondElement();
		double[][] pvals = pair.getFirstElement()[1];
		for(int i=0;i<pos.length;i++){
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
				ResultSet rs = new ResultSet( new Result[]{ 
				                                           new NumericalResult( "Sequence index", "The index of the sequence", i+1 ),
				                                           new NumericalResult( "Position", "The starting position of the motif within the sequence", curr+1 ),
				                                           new CategoricalResult( "Strand", "The strand of the predicted BS", rc ? "-" : "+" ),
				                                           new NumericalResult( "p-value", "The p-value of the predicted BS", pvals[i][j] ),
				                                           new CategoricalResult( "Binding site", "The binding site as in the sequence", sub.toString() ),
				                                           new CategoricalResult( "Adjusted binding site", "The binding site in predicted orientation", sub2.toString() ),
				                                           new CategoricalResult( "Sequence annotation", "The annotation of the original sequence", pars.parseAnnotationToComment( ' ', data.getElementAt( i ).getAnnotation() ).substring( 1 ) )
				} );
				set.add( rs );
			}
		}
		
		ListResult lr = new ListResult( "Predictions for motif "+(motifIndex+1), "", null, set.toArray( new ResultSet[0] ) );
		lr.setExport(true);
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

	private static DataSet annotate( Sequence[] annotated, double[][] weights, double[] mean, double sd, boolean[][] allowed ) throws WrongAlphabetException, WrongSequenceTypeException, EmptyDataSetException {
		AlphabetContainer ref = new AlphabetContainer( new ContinuousAlphabet() );
		float[][] histogram = new float[annotated.length][];
		for( int j = 0; j < weights[0].length; j++ ) {
			histogram[j] = new float[annotated[j].getLength()];
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
					histBg[i] = max - histBg[i];
					sum += histBg[i];
				}
			}
			for(int i=0;i<histBg.length;i++){
				histBg[i] /= sum;
			}
			
			for( int i = 0; i < histogram[j].length; i++ ) {
				histogram[j][i] = (float)(weights[0][j] * histogram[j][i] + weights[1][j] * histBg[i]);
			}
			annotated[j] = annotated[j].annotate( false, new ReferenceSequenceAnnotation( "reads", new ArbitraryFloatSequence( ref, histogram[j] ) ) );
		}
		return new DataSet( "", annotated );
	}

	private static boolean heuristic( MutableMotifDiscoverer best, DataSet completeData, double[][] completeWeight, LogGenDisMixFunction objective, double[] mean, double[] sds, Protocol protocol ) throws Exception {
		SignificantMotifOccurrencesFinder smof = new SignificantMotifOccurrencesFinder(best,completeData,completeWeight[1],ALPHA);//XXX Jan
		boolean modified = false;
		double log2 = Math.log( 2 );
		for(int im=0;im<best.getNumberOfMotifs();im++){

			Pair<double[][],double[]> pair = smof.getPWMAndPosDist( im, completeData, completeWeight[0], mean, 0, 0 );
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
			
			int left = 0, right = 0;
			double thresh = 0.2;
			while( left < kl.length &&  kl[left] < thresh ){
				left++;
			}			
			while( right > -kl.length && kl[kl.length-1+right] < thresh ){
				right--;
			}
			if(left == 0 && kl[0]>4*thresh){
				left--;
			}
			if(right == 0 && kl[kl.length-1]>4*thresh){
				right++;
			}
			protocol.append( "left: "+ left + ", right: " + right+"\n" );
			if( left == kl.length || right == -kl.length ) {
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
					protocol.append( "no modifications for the motif\n");
				}
			}
		}
		return modified;
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
			current = ((MarkovModelDiffSM) chipper.getFunction(0)).getPWM();
			
			//Filter for Poly-N
			double max = Double.NEGATIVE_INFINITY;
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
			protocol.append("added: "+getConsensus(chipper.getAlphabetContainer(), current)+"\n");
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
			current = ((MarkovModelDiffSM) chipper.getFunction(0)).getPWM();
			protocol.append("added: "+getConsensus(chipper.getAlphabetContainer(), current)+"\n");
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
			current = ((MarkovModelDiffSM) ((AbstractSingleMotifChIPper)disc[order[i]]).getFunction(0)).getPWM();
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
	public ToolResult[] getTestCases(String path) {
		try {
			return new ToolResult[]{
					new ToolResult(FileManager.readFile(path+File.separator+"xml/dimont.xml"))};
		} catch( Exception e ) {
			e.printStackTrace();
			return null;
		}
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
				"	Year = {2013}}"};
	}

	
	
}
