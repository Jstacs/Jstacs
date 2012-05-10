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


public class TrainPWM {

	/**
	 * @param args 
	 * <ul>
	 * <li>args[0] contains the path to the training data set</li>
	 * </ul>
	 */
	public static void main( String[] args ) throws Exception {
		
		//read data from FastA file
		DNADataSet ds = new DNADataSet( args[0] );
		AlphabetContainer con = ds.getAlphabetContainer();
		//create position weight matrix model
		TrainableStatisticalModel pwm = TrainableStatisticalModelFactory.createPWM( con, ds.getElementLength(), 4 );
		//train it on the input data
		pwm.train( ds );
		//print the trained model
		System.out.println(pwm);
		
		
	}

}
