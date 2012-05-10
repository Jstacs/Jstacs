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

import de.jstacs.classifiers.performanceMeasures.NumericalPerformanceMeasureParameterSet;
import de.jstacs.classifiers.performanceMeasures.PerformanceMeasureParameterSet;
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
		
		//evaluate
		NumericalPerformanceMeasureParameterSet params = PerformanceMeasureParameterSet.createFilledParameters();
		System.out.println( cl.evaluate(params, true, data) );
	}

}
