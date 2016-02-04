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

package projects.dimont;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.OneDataSetLogGenDisMixFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;

public class HeuristicOneDataSetLogGenDisMixFunction extends OneDataSetLogGenDisMixFunction {
	
	public HeuristicOneDataSetLogGenDisMixFunction(int threads, DifferentiableSequenceScore[] score, DataSet data, double[][] weights, LogPrior prior, double[] beta, boolean norm, boolean freeParams) throws IllegalArgumentException {
		super(threads, score, data, weights, prior, beta, norm, freeParams);
	}

	public void resetHeuristics() {
		for( int i = 0; i < score.length; i++ ) {
			for( int j = 0; j < score[i].length; j++ ) {
				if( score[i][j] instanceof AbstractSingleMotifChIPper ) {
					//System.out.println("reset " + i + " " + j );
					((AbstractSingleMotifChIPper) score[i][j]).reset();
				}
			}
		}
	}	
}
