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

package de.jstacs.sequenceScores.statisticalModels.differentiable;


/**
 * Interface for {@link DifferentiableStatisticalModel}s that can be used for
 * Metropolis-Hastings sampling in a {@link de.jstacs.classifiers.differentiableSequenceScoreBased.sampling.SamplingScoreBasedClassifier}.
 * 
 * @author Jan Grau
 *
 */
public interface SamplingDifferentiableStatisticalModel extends DifferentiableStatisticalModel {

	/**
	 * Returns groups of indexes of parameters that shall be drawn
	 * together in a sampling procedure
	 * @param parameterOffset a global offset on the parameter indexes
	 * @return the groups of indexes. The first dimension represents the different groups while the second
	 * 			dimension contains the parameters that shall be sampled together. Internal parameter indexes
	 * 			need to be increased by <code>parameterOffset</code>.
	 */
	public int[][] getSamplingGroups(int parameterOffset);
	
}
