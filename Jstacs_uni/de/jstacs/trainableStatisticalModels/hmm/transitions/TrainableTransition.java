/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.trainableStatisticalModels.hmm.transitions;

/**
 * This class declares methods that allow for estimating the parameters from a sufficient statistic,
 * as for instance done in the (modified) Baum-Welch algorithm, viterbi training, or Gibbs sampling.
 * 
 * @author Jens Keilwagen
 */
public interface TrainableTransition extends TransitionWithSufficientStatistic {

	/**
	 * This method estimates the parameter of the transition using the internal sufficient statistic.
	 */
	public abstract void estimateFromStatistic();
}
