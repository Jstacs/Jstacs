package supplementary.cookbook.recipes;

import de.jstacs.data.DNADataSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModelFactory;


public class TrainHomogeneousMM {

	/**
	 * @param args 
	 * <ul>
	 * <li>args[0] contains the path to the training data set</li>
	 * </ul>
	 */
	public static void main( String[] args ) throws Exception {
		
		//read data from FastA file
		DNADataSet ds = new DNADataSet( args[0] );
		//create homogeneous Markov model of order 1
		TrainableStatisticalModel hmm = TrainableStatisticalModelFactory.createHomogeneousMarkovModel( ds.getAlphabetContainer(), 4, (byte)1 );
		//train it on the input data
		hmm.train( ds );
		//print the trained model
		System.out.println(hmm);
	}

}
