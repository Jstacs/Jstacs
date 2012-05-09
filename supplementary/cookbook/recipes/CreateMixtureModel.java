package supplementary.cookbook.recipes;

import de.jstacs.data.DNADataSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModelFactory;


public class CreateMixtureModel {

	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		
		//read data from FastA file
		DNADataSet ds = new DNADataSet( args[0] );
		//create position weight matrix model
		TrainableStatisticalModel pwm = TrainableStatisticalModelFactory.createPWM( ds.getAlphabetContainer(), ds.getElementLength(), 4 );
		//create mixture model of two position weight matrices
		TrainableStatisticalModel mixture = TrainableStatisticalModelFactory.createMixtureModel( new double[]{4,4}, new TrainableStatisticalModel[]{pwm,pwm} );
		//train it on the input data using EM
		mixture.train( ds );
		//print the trained model
		System.out.println(mixture);
		
	}

}
