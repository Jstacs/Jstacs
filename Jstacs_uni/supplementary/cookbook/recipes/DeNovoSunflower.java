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

import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.HMMFactory;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.BaumWelchParameterSet;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;


public class DeNovoSunflower {

	/**
	 * @param args 
	 * <ul>
	 * <li>args[0] contains the path to the training data set, which typicall contains longer sequences with hidden motifs</li>
	 * </ul>
	 */
	public static void main(String[] args) throws Exception {
		//load data
		DataSet data = new DNADataSet(args[0]);
		//define parameters of Baum-Welch training
		BaumWelchParameterSet pars = new BaumWelchParameterSet(10, new SmallDifferenceOfFunctionEvaluationsCondition(1E-6), 2);
		//create sunflower HMM with motifs of length 8 and 12
		AbstractHMM hmm = HMMFactory.createSunflowerHMM(pars, data.getAlphabetContainer(), 0, data.getElementLength(), true, 8,12);
		//train the HMM using Baum-Welch
		hmm.train(data);
		//print the trained HMM
		System.out.println(hmm);
		//print Viterbi paths of all sequences
		for(int i=0;i<data.getNumberOfElements();i++){
			Pair<IntList,Double> p = hmm.getViterbiPathFor(data.getElementAt(i));
			System.out.println(p.getSecondElement()+"\t"+p.getFirstElement());
		}
		
	}

}
