package supplementary.cookbook.recipes;

import de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModelFactory;


public class TrainClassifier {

	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		//read data from FastA files
		DataSet[] data = new DataSet[2];
		data[0] = new DNADataSet( args[0] );
		data[1] = new DNADataSet( args[1] );

		//create position weight matrix model
		TrainableStatisticalModel pwm = TrainableStatisticalModelFactory.createPWM( data[0].getAlphabetContainer(), data[0].getElementLength(), 4 );
		//create homogeneous Markov model of order 1
		TrainableStatisticalModel hmm = TrainableStatisticalModelFactory.createHomogeneousMarkovModel( data[1].getAlphabetContainer(), 4, (byte)1 );
		//build a classifier using these models
		TrainSMBasedClassifier cl = new TrainSMBasedClassifier( pwm, hmm );
		//train it on the training data
		cl.train( data );
		
		//print the trained classifier
		System.out.println(cl);
		//classify one of the sequences
		Sequence seq = data[0].getElementAt( 0 );
		byte res = cl.classify( seq );
		//print sequence and classification result
		System.out.println(seq+" -> "+res);
	}

}
