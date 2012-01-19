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

package supplementary.codeExamples;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;

import javax.swing.JFrame;

import org.biojavax.bio.db.ncbi.GenbankRichSequenceDB;
import org.rosuda.REngine.REXP;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.classifiers.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.classifiers.assessment.KFoldCrossValidation;
import de.jstacs.classifiers.assessment.KFoldCrossValidationAssessParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.AbstractMultiThreadedOptimizableFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.CompositeLogPrior;
import de.jstacs.classifiers.differentiableSequenceScoreBased.msp.MSPClassifier;
import de.jstacs.classifiers.performanceMeasures.AbstractPerformanceMeasure;
import de.jstacs.classifiers.performanceMeasures.NumericalPerformanceMeasureParameterSet;
import de.jstacs.classifiers.performanceMeasures.PRCurve;
import de.jstacs.classifiers.performanceMeasures.PerformanceMeasureParameterSet;
import de.jstacs.classifiers.performanceMeasures.ROCCurve;
import de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.DataSet.PartitionMethod;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.GenericComplementableDiscreteAlphabet;
import de.jstacs.data.bioJava.BioJavaAdapter;
import de.jstacs.data.bioJava.SimpleSequenceIterator;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.FileManager;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.StringExtractor;
import de.jstacs.results.ImageResult;
import de.jstacs.results.ResultSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.MixtureDiffSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.BayesianNetworkTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.BayesianNetworkTrainSMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.MixtureTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization;
import de.jstacs.utils.REnvironment;


/**
 * This class contains some code examples for the Jstacs homepage.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class CodeExampleTest {

	private static String server="localhost";
	private static String login="";
	private static String password="";
	private static String home = "./supplementary/codeExamples/";
	
	/**
	 * @param args 0: fg file, 1: bg file, 2: test file
	 */
	public static void main( String[] args ) throws Exception {
		/*
		loadData();
		System.out.println("loadData");
		trainClassifier( args );
		System.out.println("trainClassifier");
		saveModel( args );
		System.out.println("saveModel");
		crossValidation( args );
		System.out.println("crossValidation");
		testREnvironment();
		System.out.println("testREnvironment");
		plotCurves( args );
		System.out.println("plotCurves");
		parameterLearning_MCL_and_MSP( args );
		System.out.println("discriminativeLearning_old");
		parameterLearning_genDisMix( args );
		System.out.println( "parameterLearning_genDis" );
		*/
		createAlphabets();
		System.out.println( "alphabets created" );
	}
	

	
	public static void createAlphabets() throws Exception {
		String[] symbols = {"A", "C", "G", "T", "-"};
		//new alphabet
		DiscreteAlphabet abc = new DiscreteAlphabet( true, symbols );
		
		//new alphabet that allows to build the reverse complement of a sequence
		int[] revComp = new int[symbols.length];
		revComp[0] = 3; //symbols[0]^rc = symbols[3]
		revComp[1] = 2; //symbols[1]^rc = symbols[2]
		revComp[2] = 1; //symbols[2]^rc = symbols[1]
		revComp[3] = 0; //symbols[3]^rc = symbols[0]
		revComp[4] = 4; //symbols[4]^rc = symbols[4]
		GenericComplementableDiscreteAlphabet abc2 = new GenericComplementableDiscreteAlphabet( true, symbols, revComp );
		
		Sequence seq = Sequence.create( new AlphabetContainer( abc2 ), "ACGT-" );
		Sequence rc = seq.reverseComplement();
		System.out.println( seq );
		System.out.println( rc );		
	}
	
	public static void parameterLearning_MCL_and_MSP( String[] args ) throws Exception {
		//read FastA-files
		DataSet[] data = {
		         new DNADataSet( args[0] ),
		         new DNADataSet( args[1] )
		};
		AlphabetContainer container = data[0].getAlphabetContainer();
		
		//equivalent sample size =^= ESS
		double essFg = 4, essBg = 4;
		//create DifferentiableSequenceScore, here PWM
		DifferentiableStatisticalModel pwmFg = new BayesianNetworkDiffSM( container, data[0].getElementLength(), essFg, true, new InhomogeneousMarkov(0) );
		DifferentiableStatisticalModel pwmBg = new BayesianNetworkDiffSM( container, data[1].getElementLength(), essBg, true, new InhomogeneousMarkov(0) );
		
		//create parameters of the classifier
		int threads = AbstractMultiThreadedOptimizableFunction.getNumberOfAvailableProcessors();
		GenDisMixClassifierParameterSet cps = new GenDisMixClassifierParameterSet( container, data[0].getElementLength(),
				//optimization parameter
				Optimizer.QUASI_NEWTON_BFGS, 1E-9, 1E-11, 1, false, KindOfParameter.ZEROS, true, threads
		);
		//create classifiers
		MSPClassifier[] cl = {
		         //MCL
		         new MSPClassifier( cps, pwmFg, pwmBg ),
		         //MSP with composite prior (here this equivalent to a transformed product-Dirichlet)
		         new MSPClassifier( cps, new CompositeLogPrior(), pwmFg, pwmBg ),
		};
		
		//do what ever you like
		
		//e.g., train
		for( int i = 0; i < cl.length; i++ ){
			cl[i].train( data );
		}
		
		//e.g., evaluate (normally done on a test data set)
		PerformanceMeasureParameterSet mp = PerformanceMeasureParameterSet.createFilledParameters();
		for( int i = 0; i < cl.length; i++ ){
			System.out.println( cl[i].evaluate( mp, true, data ) );
		}			
	}
	
	public static void parameterLearning_genDisMix( String[] args ) throws Exception {
		//read FastA-files
		DataSet[] data = {
		         new DNADataSet( args[0] ),
		         new DNADataSet( args[1] )
		};
		AlphabetContainer container = data[0].getAlphabetContainer();
		int length = data[0].getElementLength();
		
		//equivalent sample size =^= ESS
		double essFg = 4, essBg = 4;
		//create DifferentiableSequenceScore, here PWM
		DifferentiableStatisticalModel pwmFg = new BayesianNetworkDiffSM( container, length, essFg, true, new InhomogeneousMarkov(0) );
		DifferentiableStatisticalModel pwmBg = new BayesianNetworkDiffSM( container, length, essBg, true, new InhomogeneousMarkov(0) );
		
		//create parameters of the classifier
		GenDisMixClassifierParameterSet cps = new GenDisMixClassifierParameterSet(
				container,//the used alphabets
				length,//sequence length that can be modeled/classified
				Optimizer.QUASI_NEWTON_BFGS, 1E-9, 1E-11, 1,//optimization parameter
				false,//use free parameters or all
				KindOfParameter.PLUGIN,//how to start the numerical optimization
				true,//use a normalized objective function
				AbstractMultiThreadedOptimizableFunction.getNumberOfAvailableProcessors()//number of compute threads		
		);
		
		//create classifiers
		LearningPrinciple[] lp = LearningPrinciple.values();
		GenDisMixClassifier[] cl = new GenDisMixClassifier[lp.length+1];
		//elementary learning principles
		int i = 0;
		for( ; i < cl.length-1; i++ ){
			System.out.println( "classifier " + i + " uses " + lp[i] );
			cl[i] = new GenDisMixClassifier( cps, new CompositeLogPrior(), lp[i], pwmFg, pwmBg );
		}
		
		//use some weighted version of log conditional likelihood, log likelihood, and log prior
		double[] beta = {0.3,0.3,0.4};
		System.out.println( "classifier " + i + " uses the weights " + Arrays.toString( beta ) );
		cl[i] = new GenDisMixClassifier( cps, new CompositeLogPrior(), beta, pwmFg, pwmBg );
		
		//do what ever you like
		
		//e.g., train
		for( i = 0; i < cl.length; i++ ){
			cl[i].train( data );
		}
		
		//e.g., evaluate (normally done on a test data set)
		PerformanceMeasureParameterSet mp = PerformanceMeasureParameterSet.createFilledParameters();
		for( i = 0; i < cl.length; i++ ){
			System.out.println( cl[i].evaluate( mp, true, data ) );
		}			
	}
	
	
	public static void plotCurves(String[] args) throws Exception {

		//read XML-representation from disk
		StringBuffer buf2 = FileManager.readFile( new File(home+"myClassifier.xml") );
		 
		//create new classifier from read StringBuffer containing XML-code
		TrainSMBasedClassifier trainedClassifier = new TrainSMBasedClassifier(buf2);	

		//create a DataSet for each class from the input data, using the DNA alphabet
		DataSet[] test = new DataSet[2];
		test[0] = new DNADataSet( args[0] );
		
		//the length of our input sequences
		int length = test[0].getElementLength();

		test[1] = new DataSet( new DNADataSet( args[1] ), length );
		
		 
		AbstractPerformanceMeasure[] m = { new PRCurve(), new ROCCurve() };
		PerformanceMeasureParameterSet mp = new PerformanceMeasureParameterSet( 2, m );
		ResultSet rs = trainedClassifier.evaluate( mp, true, test );
		 
		REnvironment r = null;
		try {
			r = new REnvironment( server, login, password );
			for( AbstractPerformanceMeasure s : m )  {
				DoubleTableResult dtr = (DoubleTableResult) rs.getResultAt( 1 );
				ImageResult ir = DoubleTableResult.plot( r, dtr );
				REnvironment.showImage( dtr.getName(), ir.getValue() );
			}
		} catch( Exception e ) {
			e.printStackTrace();
		} finally {
			if( r != null ) {
				r.close();
			}
		}
	}
	
	public static void testREnvironment(){
		REnvironment e = null;
		try {
			//create a connection to R with YOUR server name, login and password
			e = new REnvironment( server, login, password );
		 
			System.out.println( "java: " + System.getProperty( "java.version" ) );
			System.out.println();
			System.out.println( e.getVersionInformation() );
		 
			// compute something in R
			REXP erg = e.eval( "sin(10)" );
			System.out.println( erg.asDouble() );
		 
			//create a histrgram in R in 3 steps
			//1) create the data
			e.voidEval( "a = 100;" );
			e.voidEval( "n = rnorm(a)" );
			//2) create the plot command
			String plotCmd = "hist(n,breaks=a/5)";
			//3a) plot as pdf
			e.plotToPDF( plotCmd, home+"test.pdf", true );
			//or
			//3b) create an image and show it
			BufferedImage i = e.plot( plotCmd, 640, 480 );
			REnvironment.showImage( "histogramm", i, JFrame.EXIT_ON_CLOSE );
		 
		} catch ( Exception ex ) {
			ex.printStackTrace();
		} finally {
			if( e != null ) {
				try {
					//close REnvironment correctly
					e.close();
				} catch ( Exception e1 ) {
					System.err.println( "could not close REnvironment." );
					e1.printStackTrace();
				}
			}
		}
	}
	
	public static void crossValidation(String[] args) throws Exception {
		//create a DataSet for each class from the input data, using the DNA alphabet
		DataSet[] data = new DataSet[2];
		data[0] = new DNADataSet( args[0] );
		
		//the length of our input sequences
		int length = data[0].getElementLength();

		data[1] = new DataSet( new DNADataSet( args[1] ), length );
		 
		AlphabetContainer container = data[0].getAlphabetContainer();
		
		//create a new PWM
		BayesianNetworkTrainSM pwm = new BayesianNetworkTrainSM( new BayesianNetworkTrainSMParameterSet(
				//the alphabet and the length of the model:
				container, length, 
				//the equivalent sample size to compute hyper-parameters
				4, 
				//some identifier for the model
				"my PWM", 
				//we want a PWM, which is an inhomogeneous Markov model (IMM) of order 0
				ModelType.IMM, (byte) 0, 
				//we want to estimate the MAP-parameters
				LearningType.ML_OR_MAP ) );
		 
		//create a new mixture model using 2 PWMs
		MixtureTrainSM mixPwms = new MixtureTrainSM(
				//the length of the mixture model
				length, 
				//the two components, which are PWMs
				new TrainableStatisticalModel[]{pwm,pwm},
				//the number of starts of the EM
				10,
				//the equivalent sample sizes
				new double[]{pwm.getESS(),pwm.getESS()},
				//the hyper-parameters to draw the initial sequence-specific component weights (hidden variables)
				1,
				//stopping criterion
				new SmallDifferenceOfFunctionEvaluationsCondition(1E-6),
				//parameterization of the model, LAMBDA complies with the
				//parameterization by log-probabilities
				Parameterization.LAMBDA);
		 
		//create a new inhomogeneous Markov model of order 3
		BayesianNetworkTrainSM mm = new BayesianNetworkTrainSM( 
				new BayesianNetworkTrainSMParameterSet( container, length, 256, "my iMM(3)", ModelType.IMM, (byte) 3, LearningType.ML_OR_MAP ) );
		 
		//create a new PWM scoring function
		BayesianNetworkDiffSM dPwm = new BayesianNetworkDiffSM(
				//the alphabet and the length of the scoring function
				container, length, 
				//the equivalent sample size for the plug-in parameters
				4, 
				//we use plug-in parameters
				true, 
				//a PWM is an inhomogeneous Markov model of order 0
				new InhomogeneousMarkov(0));
		 
		//create a new mixture scoring function
		MixtureDiffSM dMixPwms = new MixtureDiffSM(
				//the number of starts
				2,
				//we use plug-in parameters
				true,
				//the two components, which are PWMs
				dPwm,dPwm);
		 
		//create a new scoring function that is an inhomogeneous Markov model of order 3
		BayesianNetworkDiffSM dMm = new BayesianNetworkDiffSM(container, length, 4, true, new InhomogeneousMarkov(3));
		 
		//create the classifiers
		int threads = AbstractMultiThreadedOptimizableFunction.getNumberOfAvailableProcessors();
		AbstractScoreBasedClassifier[] classifiers = new AbstractScoreBasedClassifier[]{
									   //model based with mixture model and Markov model
									   new TrainSMBasedClassifier( mixPwms, mm ),
									   //conditional likelihood based classifier
									   new MSPClassifier( new GenDisMixClassifierParameterSet(container, length, 
											   //method for optimizing the conditional likelihood and 
											   //other parameters of the numerical optimization
											   Optimizer.QUASI_NEWTON_BFGS, 1E-2, 1E-2, 1, true, KindOfParameter.PLUGIN, false, threads),
											   //mixture scoring function and Markov model scoring function
											   dMixPwms,dMm )
		};
		 
		//create an new k-fold cross validation using above classifiers
		KFoldCrossValidation cv = new KFoldCrossValidation( classifiers );
		 
		//we use a specificity of 0.999 to compute the sensitivity and a sensitivity of 0.95 to compute FPR and PPV
		NumericalPerformanceMeasureParameterSet mp = (NumericalPerformanceMeasureParameterSet) PerformanceMeasureParameterSet.createFilledParameters();
		//we do a 10-fold cross validation and partition the data by means of the number of symbols
		KFoldCrossValidationAssessParameterSet cvpars = new KFoldCrossValidationAssessParameterSet(PartitionMethod.PARTITION_BY_NUMBER_OF_SYMBOLS, length, true, 2);
		 
		//compute the result of the cross validation and print them to System.out
		System.out.println( cv.assess( mp, cvpars, data ) );
	}
	
	public static void saveModel(String[] args) throws Exception {
		//create a DataSet for each class from the input data, using the DNA alphabet
		DataSet[] data = new DataSet[2];
		data[0] = new DNADataSet( args[0] );
		
		//the length of our input sequences
		int length = data[0].getElementLength();

		data[1] = new DataSet( new DNADataSet( args[1] ), length );
		 
		//create a new PWM
		BayesianNetworkTrainSM pwm = new BayesianNetworkTrainSM( new BayesianNetworkTrainSMParameterSet(
				//the alphabet and the length of the model:
				data[0].getAlphabetContainer(), length, 
				//the equivalent sample size to compute hyper-parameters
				4, 
				//some identifier for the model
				"my PWM", 
				//we want a PWM, which is an inhomogeneous Markov model (IMM) of order 0
				ModelType.IMM, (byte) 0, 
				//we want to estimate the MAP-parameters
				LearningType.ML_OR_MAP ) );
		 
		//create a new classifier
		TrainSMBasedClassifier classifier = new TrainSMBasedClassifier( pwm, pwm );
		 
		//train the classifier
		classifier.train( data );
		 
		//create the XML-representation of the classifier
		StringBuffer buf = classifier.toXML();
		 
		//write it to disk
		FileManager.writeFile( new File(home+"myClassifier.xml"), buf );
		 
		//read XML-representation from disk
		StringBuffer buf2 = FileManager.readFile( new File(home+"myClassifier.xml") );
		 
		//create new classifier from read StringBuffer containing XML-code
		TrainSMBasedClassifier loaded = new TrainSMBasedClassifier(buf2);
	}
	
	public static void trainClassifier(String[] args) throws Exception{
				 
		//create a DataSet for each class from the input data, using the DNA alphabet
		DataSet[] data = new DataSet[2];
		data[0] = new DNADataSet( args[0] );
		
		//the length of our input sequences
		int length = data[0].getElementLength();

		data[1] = new DataSet( new DNADataSet( args[1] ), length );

		
		//sequences that will be classified
		DataSet toClassify = new DNADataSet( args[2] );
		 
		//create a new PWM
		BayesianNetworkTrainSM pwm = new BayesianNetworkTrainSM( new BayesianNetworkTrainSMParameterSet(
				//the alphabet and the length of the model:
				data[0].getAlphabetContainer(), length, 
				//the equivalent sample size to compute hyper-parameters
				4, 
				//some identifier for the model
				"my PWM", 
				//we want a PWM, which is an inhomogeneous Markov model (IMM) of order 0
				ModelType.IMM, (byte) 0, 
				//we want to estimate the MAP-parameters
				LearningType.ML_OR_MAP ) );
		 
		//create a classifier with a PWM in the foreground and a PWM in the background
		TrainSMBasedClassifier classifier = new TrainSMBasedClassifier( pwm, pwm );
		 
		//train the classifier
		classifier.train( data );
		 
		//use the trained classifier to classify new sequences
		byte[] result = classifier.classify( toClassify );
		 
		System.out.println( Arrays.toString( result ) );
	}
	
	public static void loadData() throws Exception {
		//load DNA sequences in FastA-format
		DataSet data = new DNADataSet( home+"myfile.fa" ); 
		
		//create a DNA-alphabet
		AlphabetContainer container = DNAAlphabetContainer.SINGLETON;
		
		//create a DataSet using the alphabet from above in FastA-format
		data = new DataSet( container, new SparseStringExtractor( home+"myfile.fa", StringExtractor.FASTA ));
		
		//create a DataSet using the alphabet from above
		data = new DataSet( container, new SparseStringExtractor( home+"myfile.txt" ));
		
		//defining the ids, we want to obtain from NCBI Genbank:
		GenbankRichSequenceDB db = new GenbankRichSequenceDB();
		
		SimpleSequenceIterator it = new SimpleSequenceIterator(
				db.getRichSequence( "NC_001284.2" ),
				db.getRichSequence( "NC_000932.1" )
				);
		 
		//conversion to Jstacs DataSet
		data = BioJavaAdapter.sequenceIteratorToDataSet( it, null );
	}
}
