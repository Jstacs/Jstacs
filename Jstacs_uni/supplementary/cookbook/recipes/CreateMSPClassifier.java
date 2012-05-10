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

import de.jstacs.classifiers.AbstractClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.msp.MSPClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModelFactory;


public class CreateMSPClassifier {

	/**
	 * @param args 
	 * <ul>
	 * <li>args[0] contains the path to the foreground data set</li>
	 * <li>args[1] contains the path to the background data set</li>
	 * </ul>
	 */
	public static void main( String[] args ) throws Exception {
		//read data from FastA files
		DataSet[] data = new DataSet[2];
		data[0] = new DNADataSet( args[0] );
		data[1] = new DNADataSet( args[1] );
		AlphabetContainer con = data[0].getAlphabetContainer();
		
		//define differentiable PWM model
		DifferentiableStatisticalModel pwm = DifferentiableStatisticalModelFactory.createPWM(con, 10, 4);
		//parameters for numerical optimization
		GenDisMixClassifierParameterSet pars = new GenDisMixClassifierParameterSet(con,10,(byte)10,1E-9,1E-10,1, false,KindOfParameter.PLUGIN,true,1);
		//define and train classifier
		AbstractClassifier cl = new MSPClassifier( pars, pwm, pwm );
		cl.train( data );
		
		System.out.println(cl);
	}
}