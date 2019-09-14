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

import java.io.IOException;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.CompositeLogPrior;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.FileManager;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ImageResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.UniformDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.MarkovModelDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.btMeasures.BTExplainingAwayResidual;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.btMeasures.BTMutualInformation;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.btMeasures.BTMutualInformation.DataSource;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.UniformHomogeneousDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.localMixture.LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder;
import de.jstacs.sequenceScores.statisticalModels.differentiable.localMixture.LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder.PriorType;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM.InitMethod;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter;
import de.jstacs.utils.graphics.GraphicsAdaptor;
import projects.dimont.Interpolation;
import projects.dimont.ThresholdedStrandChIPper;
import projects.slim.LearnDependencyModelWebParameterSet.ModelType;


public class LearnDependencyModelTool implements JstacsTool {

	
	public static class DepLogoPlotGenerator implements PlotGenerator{

		private DataSet data;
		private double[] weights;
		
		public DepLogoPlotGenerator(DataSet data, double[] weights){
			this.data = data;
			this.weights = weights;
		}
		
		public DepLogoPlotGenerator(StringBuffer xml) throws NonParsableException {
			xml = XMLParser.extractForTag(xml, "DepLogoPlotGenerator");
			try {
				data = new DataSet("",XMLParser.extractSequencesWithTags(xml, "data"));
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
			weights = (double[]) XMLParser.extractObjectForTags(xml, "weights");
		}
		
		@Override
		public StringBuffer toXML() {
			StringBuffer xml = new StringBuffer();
			XMLParser.appendSequencesWithTags(xml, "data", data.getAllElements());
			XMLParser.appendObjectWithTags(xml, weights, "weights");
			XMLParser.addTags(xml, "DepLogoPlotGenerator");
			return xml;
		}

		@Override
		public void generatePlot(GraphicsAdaptor ga) throws Exception {
			
			SeqLogoPlotter.plotDefaultDependencyLogoToGraphicsAdaptor(ga, data, weights, 600);
			
		}
		
	}
	
	
	
	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		
		CLI cli = new CLI(new boolean[]{true,false},new LearnDependencyModelTool(), new MotifScanningTool());
		
		cli.run(args);

	}



	@Override
	public ToolParameterSet getToolParameters() {
		try {
			return new LearnDependencyModelWebParameterSet();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}


	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		LearnDependencyModelWebParameterSet params = (LearnDependencyModelWebParameterSet) parameters;		
		
		Pair<DataSet,double[]> pair = params.getData();
		
		
		DataSet data = pair.getFirstElement();
		double[] signals = pair.getSecondElement();
		
		double ess = params.getESS();
		ModelType modelType = params.getModelType();
		int order = 0;
		if(modelType != ModelType.BT_EAR && modelType != ModelType.BT_MI && modelType != ModelType.SLIM){
			order = params.getOrder();
		}
		int bgO = params.getBgOrder();
		
		
		LinkedList<Result> result = new LinkedList<Result>();
		
		double[] weights = Interpolation.getWeight( data, signals, 0.5, Interpolation.PERCENTILE_LOGISTIC );
		
		DifferentiableStatisticalModel model = null;
		if(modelType == ModelType.IMM){
			model = new MarkovModelDiffSM( DNAAlphabetContainer.SINGLETON, data.getElementLength(), ess, true, new InhomogeneousMarkov( order ) );
		}else if(modelType == ModelType.BT_EAR){
			model = new BayesianNetworkDiffSM( DNAAlphabetContainer.SINGLETON, data.getElementLength(), ess, true, new BTExplainingAwayResidual( new double[]{ess,ess} ) );
		}else if(modelType == ModelType.BT_MI){
			model = new BayesianNetworkDiffSM( DNAAlphabetContainer.SINGLETON, data.getElementLength(), ess, true, new BTMutualInformation( DataSource.FG, new double[]{ess,ess} ) );
		}else if(modelType == ModelType.SLIM){
			model = new LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder(DNAAlphabetContainer.SINGLETON, data.getElementLength(), 1, data.getElementLength(), ess, 0.9, PriorType.BDeu);
		}else if(modelType == ModelType.LSLIM){
			model = new LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder(DNAAlphabetContainer.SINGLETON, data.getElementLength(), 1, -order, ess, 0.9, PriorType.BDeu);
		}else{
			throw new RuntimeException( "Model type unknown." );
		}
		
		model = new StrandDiffSM( model, 0.5, 1, true, InitMethod.INIT_FORWARD_STRAND );
		
		DifferentiableStatisticalModel bg = null;
		if(bgO >= 0){
			bg = new HomogeneousMMDiffSM( DNAAlphabetContainer.SINGLETON, bgO, ess, data.getElementLength() );
		}else if(bgO == -1){
			bg = new UniformDiffSM( DNAAlphabetContainer.SINGLETON, data.getElementLength(), ess );
		}else{
			throw new RuntimeException( "Illegal background order." );
		}
		
		
		GenDisMixClassifierParameterSet gdmp = new GenDisMixClassifierParameterSet( DNAAlphabetContainer.SINGLETON, data.getElementLength(), Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1E-6, 1E-4, false, KindOfParameter.PLUGIN, true, threads );
		
		GenDisMixClassifier cl = null;
		cl = new GenDisMixClassifier( gdmp, new CompositeLogPrior(), LearningPrinciple.MSP, model, bg );
		
		cl.train( new DataSet[]{data,data}, new double[][]{weights,Interpolation.getBgWeight( weights )} );
		
		
		//result.add( new StorableResult( "Classifier", "The classifier trained from the input data", cl ) );
		
		
		SplitSequenceAnnotationParser pars = new SplitSequenceAnnotationParser( ":", ";" );

		LinkedList<ResultSet> set = new LinkedList<ResultSet>();
		
		LinkedList<Sequence> bs = new LinkedList<Sequence>();
		DoubleList bsWeights = new DoubleList();
		
		model = (DifferentiableStatisticalModel)cl.getDifferentiableSequenceScore( 0 );
		
		for(int i=0;i<data.getNumberOfElements();i++){
			
			Sequence seq = data.getElementAt( i );
			
			boolean rc = ((StrandDiffSM)model).getStrand( seq, 0 ) == Strand.REVERSE;
			
			double score = model.getLogScoreFor( seq );
			
			Sequence sub2 = seq;
			if(rc){
				sub2 = seq.reverseComplement();
			}
			
			bs.add( sub2 );
			bsWeights.add( score );
			
			ResultSet rs = new ResultSet( new Result[]{ 
			                                           new NumericalResult( "Sequence index", "The index of the sequence", i+1 ),
			                                           new NumericalResult( "Position", "The starting position of the motif within the sequence", 0 ),
			                                           new CategoricalResult( "Strand", "The strand of the predicted BS", rc ? "-" : "+" ),
			                                           new NumericalResult( "Score", "The model score of the predicted BS", score ),
			                                           new CategoricalResult( "Binding site", "The binding site as in the sequence", seq.toString() ),
			                                           new CategoricalResult( "Adjusted binding site", "The binding site in predicted orientation", sub2.toString() ),
			                                           new NumericalResult( "Signal", "The signal of the sequence annotation", signals[i] ),
			                                           new CategoricalResult( "Sequence annotation", "The annotation of the original sequence", pars.parseAnnotationToComment( ' ', seq.getAnnotation() ).substring( 1 )  )
			} );
			
			set.add( rs );
		}
		
		
		ListResult lr = new ListResult( "Predicted sequence orientations and scores", "", null, set.toArray( new ResultSet[0] ) );
		
		result.add( lr );
		
		
		//result.add( new ImageResult( "Dependency logo", "Dependency logo of sequences", SeqLogoPlotter.plotDefaultDependencyLogoToBufferedImage( new DataSet("",bs), bsWeights.toArray(), 600 ) ) );

		result.add( new PlotGeneratorResult("Dependency logo", "Dependency logo of sequences", new DepLogoPlotGenerator(new DataSet("",bs), bsWeights.toArray()), true) );
		
		
		
		StrandDiffSM sd = (StrandDiffSM)cl.getDifferentiableSequenceScore( 0 );
		
		model = sd.getFunction( 0 );
		
		
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		byte algo = Optimizer.CONJUGATE_GRADIENTS_PRP;
		double eps = 1E-4;
		
		boolean free = false;
		
		GenDisMixClassifierParameterSet genDisMixParams = new GenDisMixClassifierParameterSet( con, 0, algo, eps, eps*1E-1, 1, free, KindOfParameter.PLUGIN, true, 1);
		
		bg = new UniformHomogeneousDiffSM( con, ess );
		
		
		
		DifferentiableStatisticalModel fg =	new ThresholdedStrandChIPper( 1, 0.5, model );
		//System.out.println(Arrays.toString( pars ));
		//System.out.println(Arrays.toString( fg.getCurrentParameterValues() ));
		
		DifferentiableStatisticalModel[] score = { fg, bg };
		
		cl = new GenDisMixClassifier( genDisMixParams, new CompositeLogPrior(), Double.NaN, LearningPrinciple.getBeta(LearningPrinciple.MSP), score );
		
		
		result.add( new StorableResult( "SlimDimont classifier", "The SlimDimont classifier built from the trained motif model.", cl ) );
		
		
		/*result.add(new ListResult( "Description of columns of binding site output (see history)", "You can download the predictions for the motifs discovered as tab-separated file from the history.", null, 
				new ResultSet[]{
				                new ResultSet( new Result[]{new CategoricalResult("Column","","Sequence index"),
				                                            new CategoricalResult( "Description", "", "The index of the sequence" )} ),
	                            new ResultSet( new Result[]{new CategoricalResult("Column","","Position"),
	                                                        new CategoricalResult( "Description", "", "The start position of predicted binding site (BS) within the sequence" )} ),
				                new ResultSet( new Result[]{new CategoricalResult("Column","","Strand"),
				                                            new CategoricalResult( "Description", "", "The strand of the predicted BS" )} ),
				                new ResultSet( new Result[]{new CategoricalResult("Column","","p-value"),
				            				                new CategoricalResult( "Description", "", "The p-value of the predicted BS" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","-log10(p-value)"),
				            	                            new CategoricalResult( "Description", "", "The negative logarithm of the p-value of the predicted BS" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Score"),
				            				            	new CategoricalResult( "Description", "", "The model score of the predicted BS" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Binding site"),
				            				                new CategoricalResult( "Description", "", "The binding site as in the sequence" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Adjusted binding site"),
				            				                new CategoricalResult( "Description", "", "The binding site in predicted orientation" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Signal"),
				            				            	new CategoricalResult( "Description", "", "The signal of the sequence annotation" )} ),
				            	new ResultSet( new Result[]{new CategoricalResult("Column","","Sequence annotation"),
				            				                new CategoricalResult( "Description", "", "The annotation of the original sequence" )} )
			}
		));*/
		
		
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(result.toArray(new Result[0])), parameters, getToolName(), new Date(System.currentTimeMillis()));
		
	}


	@Override
	public String getToolName() {
		return "LearnDependencyModel";
	}


	@Override
	public String getToolVersion() {
		return "1.0";
	}


	@Override
	public String getShortName() {
		return "learn";
	}


	@Override
	public String getDescription() {
		return "Learn a dependency model from aligned input sequences";
	}


	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( LearnDependencyModelTool.class.getClassLoader().getResourceAsStream( "projects/slim/helpLearn.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
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
