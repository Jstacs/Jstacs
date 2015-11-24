package projects.slim;

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
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ImageResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
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
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter;
import projects.dimont.Interpolation;
import projects.dimont.ThresholdedStrandChIPper;
import projects.slim.LearnDependencyModelWebParameterSet.ModelType;


public class LearnDependencyModelWeb {

	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		
		LearnDependencyModelWebParameterSet params = new LearnDependencyModelWebParameterSet();
		
		boolean[] lines = new boolean[params.getNumberOfParameters()];
		//lines[0] = lines[5] = lines[6] = lines[9] = true;
		
		
		GalaxyAdaptor ga = new GalaxyAdaptor( params,lines,"LearnDependencyModel", "- Learn a dependency model from aligned input sequences", "0.1", "java -Xms256M -Xmx2G -jar "+System.getProperty( "user.dir" )+System.getProperty( "file.separator" )+"LearnDependencyModelWeb.jar", "jobname" );
		ga.setHelp( FileManager.readInputStream( SlimDimontWeb.class.getClassLoader().getResourceAsStream( "projects/slim/helpLearn.txt" ) ).toString() );//TODO
		
		if(!ga.parse( args, false )){
			System.exit( 1 );
		}	
		
		
		Pair<DataSet,double[]> pair = params.getData();
		
		Result[] res = run( pair.getFirstElement(), pair.getSecondElement(), params.getModelType(), params.getOrder(), params.getBgOrder(), params.getESS(), 1 );
		
		//ga.addResult( res[0], true, false );
		ga.addResult( res[0], true, false );
		ga.addResult( res[1], false, true );
		ga.addResult( res[2], true, false );
		
		
		ga.addResult( new ListResult( "Description of columns of binding site output (see history)", "You can download the predictions for the motifs discovered as tab-separated file from the history.", null, 
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
		), false, true );
		
		ga.writeOutput();

	}

	
	private static Result[] run(DataSet data, double[] signals, ModelType modelType, Integer order, int bgO, double ess, int threads) throws Exception {
		
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
			model = new LimitedSparseLocalInhomogeneousMixtureDiffSM_higherOrder(DNAAlphabetContainer.SINGLETON, data.getElementLength(), 1, order, ess, 0.9, PriorType.BDeu);
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
		
		
		GenDisMixClassifierParameterSet params = new GenDisMixClassifierParameterSet( DNAAlphabetContainer.SINGLETON, data.getElementLength(), Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1E-6, 1E-4, false, KindOfParameter.PLUGIN, true, threads );
		
		GenDisMixClassifier cl = null;
		cl = new GenDisMixClassifier( params, new CompositeLogPrior(), LearningPrinciple.MSP, model, bg );
		
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
		
		
		result.add( new ImageResult( "Dependency logo", "Dependency logo of sequences", SeqLogoPlotter.plotDefaultDependencyLogoToBufferedImage( new DataSet("",bs), bsWeights.toArray(), 600 ) ) );

		
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
		
		
		return result.toArray( new Result[0] );
		
	}
	
}
