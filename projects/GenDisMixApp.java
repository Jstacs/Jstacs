package projects;

import java.io.File;
import java.util.Arrays;

import de.jstacs.DataType;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.CompositeLogPrior;
import de.jstacs.data.DNADataSet;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetTagger;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;



public class GenDisMixApp {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main( String[] args ) throws Exception {
	
		ParameterSetTagger params = new ParameterSetTagger( strs, new GDMParameters() );
		try{
			params.fillParameters( "=", args );
		}catch(Exception e){
			System.out.println( "Some of the required parameters are not specified or unknown parameter specified." );
			System.out.println("You provided: "+Arrays.toString( args ));
			System.out.println("Allowed parameters are: "+Arrays.toString( strs ));
			System.exit( 1 );
		}
		System.out.println( "parameters:" );
		System.out.println( params );
		System.out.println("_________________________________");
		if( !params.hasDefaultOrIsSet() ) {
			System.out.println( "Some of the required parameters are not specified." );
			System.exit( 1 );
		}
		
		// Reading the training data from FastA-files
		DNADataSet fg = new DNADataSet( params.getValueFromTag( strs[0] )+System.getProperty( "file.separator" )+params.getValueFromTag( strs[1] ) );
		DNADataSet bg = new DNADataSet( params.getValueFromTag( strs[0] )+System.getProperty( "file.separator" )+params.getValueFromTag( strs[2] ) );

		// setting the weights for the unified learning principle
		double[] beta = new double[3];
		beta[LearningPrinciple.LIKELIHOOD_INDEX] = (Double)params.getValueFromTag( strs[3] );
		beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] = (Double)params.getValueFromTag( strs[4] );
		if(beta[LearningPrinciple.LIKELIHOOD_INDEX] + beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] > 1){
			System.out.println("The weights for the generative and discriminative components are greater than 1!");
			System.exit( 1 );
		}
		beta[LearningPrinciple.PRIOR_INDEX] = 1.0 - ( beta[LearningPrinciple.LIKELIHOOD_INDEX] + beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] );
		
		// setting the stopping criterion for numerical optimization
		double eps = (Double) params.getValueFromTag( strs[5] );
		
		GenDisMixClassifierParameterSet pars = new GenDisMixClassifierParameterSet(
				fg.getAlphabetContainer(), //alphabet from data
				fg.getElementLength(), //sequence length from data
				Optimizer.QUASI_NEWTON_BFGS, //numerical optimization algorithm 
				eps, //epsilon
				eps, // lin eps
				eps*1E2, //start distance
				false, // use only free parameters
				KindOfParameter.PLUGIN, // how to start for class parameters
				false, // normalize objective function
				(Integer)params.getValueFromTag( strs[6] ) //threads 
			);
		
		
		// setting the equivalent sample sizes
		double essFg = (Double) params.getValueFromTag( strs[7] );
		double essBg = (Double) params.getValueFromTag( strs[8] );
		
		BayesianNetworkDiffSMParameterSet modelParsFg = new BayesianNetworkDiffSMParameterSet(fg.getAlphabetContainer(),fg.getElementLength(),essFg,true,new InhomogeneousMarkov( 0 ));
		BayesianNetworkDiffSMParameterSet modelParsBg = new BayesianNetworkDiffSMParameterSet(fg.getAlphabetContainer(),fg.getElementLength(),essBg,true,new InhomogeneousMarkov( 0 ));
		
		// instantiating the classifier learned by the unified learning principle
		GenDisMixClassifier cl = new GenDisMixClassifier( pars, new CompositeLogPrior(), beta, 
				new BayesianNetworkDiffSM( modelParsFg ),
				new BayesianNetworkDiffSM( modelParsBg )
		);
		
		// training the classifier on the training data
		cl.train( fg,bg );
		
		// storing the trained classifier to disk
		FileManager.writeFile( new File( params.getValueFromTag( strs[0] )+System.getProperty( "file.separator" )+params.getValueFromTag( strs[9] ) ), cl.toXML() );
		
		String clname = (String) params.getValueFromTag( strs[10] );
		
		// classification of unclassified data
		if(clname != null) {
			
			DNADataSet uk = new DNADataSet( params.getValueFromTag( strs[0] )+System.getProperty( "file.separator" )+clname );
			
			double[] scores = cl.getScores( uk );
			byte[] cls = cl.classify( uk );
			
			System.out.println();
			System.out.println("_________________________________");
			System.out.println("index\tsequence\tscore\tclass");
			for(int i=0;i<scores.length;i++){
				System.out.println(i+"\t"+uk.getElementAt( i )+"\t"+scores[i]+"\t"+cls[i]);
			}
			System.out.println();
		}		
	}
	
	private static String[] strs = new String[]{
	                        "home",
	                        "fg",
	                        "bg",
	                        "gen",
	                        "dis",
	                        "eps",
	                        "threads",
	                        "essFG",
	                        "essBG",
	                        "outfile",
	                        "uk"
	};
	
	private static class GDMParameters extends ParameterSet {

		public GDMParameters() throws Exception {
			super();
			initParameterList(10);
			parameters.add( new SimpleParameter( DataType.STRING, "home directory", "the path to the data directory", true, "./" ) );
			parameters.add( new SimpleParameter( DataType.STRING, "foreground file", "the file name of the foreground data file in FastA-format", true ) );
			parameters.add( new SimpleParameter( DataType.STRING, "background file", "the file name of the background data file in FastA-format", true ) );
			parameters.add( new SimpleParameter( DataType.DOUBLE, "generative weight", "the weight of the generative component", true ,new NumberValidator<Comparable<? extends Number>>( 0.0, 1.0 )) );
			parameters.add( new SimpleParameter( DataType.DOUBLE, "discriminative weight", "the weight of the discriminative component", true ,new NumberValidator<Comparable<? extends Number>>( 0.0, 1.0 )) );
			parameters.add( new SimpleParameter( DataType.DOUBLE, "epsilon", "numerical optimization is stopped if gain less than epsilon", true, 1E-6 ) );
			parameters.add( new SimpleParameter( DataType.INT, "threads", "the number of threads used for the computation", true, 1 ) );
			parameters.add( new SimpleParameter( DataType.DOUBLE, "essFG", "the equivalent sample size used for the foreground class", true, 4.0 ) );
			parameters.add( new SimpleParameter( DataType.DOUBLE, "essFG", "the equivalent sample size used for the background class", true, 4.0 ) );
			parameters.add( new SimpleParameter( DataType.STRING, "outfile", "the name of the file where to store the classifier in XML-format", true , "gendismix.xml" ) );
			parameters.add( new SimpleParameter( DataType.STRING, "unkown file", "the file name of the data in FastA-format that shall be classified", false ) );
		}		
	}
}
