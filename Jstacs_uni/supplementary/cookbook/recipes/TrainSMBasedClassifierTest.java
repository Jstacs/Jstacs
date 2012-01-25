package supplementary.cookbook.recipes;

import java.io.File;
import java.util.Arrays;

import de.jstacs.classifiers.AbstractClassifier;
import de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.BayesianNetworkTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.BayesianNetworkTrainSMParameterSet;

/**
 * This class exemplarily shows how to train, use, save, and load a {@link TrainSMBasedClassifier}.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see GenDisMixClassifierTest
 */
public class TrainSMBasedClassifierTest {

	/**
	 * @param args
	 * <ul>
	 * <li>args[0] contains the path to the classifier</li>
	 * <li>args[1] contains the path to the foreground data set</li>
	 * <li>args[2] contains the path to the background data set</li>
	 * <li>args[3] contains the path to the data set to be classified</li>
	 * </ul>
	 */
	public static void main(String[] args) throws Exception {
		String home = args[0];
		
		//create a DataSet for each class from the input data, using the DNA alphabet
		DataSet[] data = new DataSet[2];
		data[0] = new DNADataSet( args[1] );
		
		//the length of our input sequences
		int length = data[0].getElementLength();

		data[1] = new DataSet( new DNADataSet( args[2] ), length );
		 
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
		
		//sequences that will be classified
		DataSet toClassify = new DNADataSet( args[2] );
		 
		//use the trained classifier to classify new sequences
		byte[] result = classifier.classify( toClassify ); 
		System.out.println( Arrays.toString( result ) );
		 
		//create the XML-representation of the classifier
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags( buf, classifier, "classifier" );
		 
		//write it to disk
		FileManager.writeFile( new File(home+"myClassifier.xml"), buf );
		 
		//read XML-representation from disk
		StringBuffer buf2 = FileManager.readFile( new File(home+"myClassifier.xml") );
		 
		//create new classifier from read StringBuffer containing XML-code
		AbstractClassifier trainedClassifier = (AbstractClassifier) XMLParser.extractObjectForTags(buf2, "classifier");	
	}
}