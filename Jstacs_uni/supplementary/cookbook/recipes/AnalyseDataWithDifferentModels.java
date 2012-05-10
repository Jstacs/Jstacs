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
package supplementary.cookbook.recipes;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModelFactory;


public class AnalyseDataWithDifferentModels {

	/**
	 * @param args 
	 * <ul>
	 * <li>args[0] contains the path to the training data set</li>
	 * </ul>
	 */
	public static void main( String[] args ) throws Exception {
		//read data from FastA file
		DNADataSet ds = new DNADataSet( args[0] );
		
		//get alphabet, length from data
		AlphabetContainer alphabet = ds.getAlphabetContainer();
		int length = ds.getElementLength();
		//set ESS used for all models
		double ess = 4;
		
		TrainableStatisticalModel[] models = new TrainableStatisticalModel[4];
		//create position weight matrix
		models[0] = TrainableStatisticalModelFactory.createPWM( alphabet, length, 4 );
		//create inhomogeneous Markov model of order 1 (WAM)
		models[1] = TrainableStatisticalModelFactory.createInhomogeneousMarkovModel( alphabet, length, ess, (byte)1 );
		//create Bayesian tree
		models[2] = TrainableStatisticalModelFactory.createBayesianNetworkModel( alphabet, length, ess, (byte)1 );
		//create homogeneous Markov model of order 2
		models[3] = TrainableStatisticalModelFactory.createHomogeneousMarkovModel( alphabet, ess, (byte)2 );
		
		//train and print all models
		for(int i=0;i<models.length;i++){
			models[i].train( ds );
			System.out.println(models[i]);
		}
	}

}
