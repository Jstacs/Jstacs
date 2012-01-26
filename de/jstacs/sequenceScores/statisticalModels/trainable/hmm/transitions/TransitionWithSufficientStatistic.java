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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions;

import de.jstacs.data.sequences.Sequence;

/**
 * This interface defines method for reseting and filling an internal sufficient statistic.
 * Using this statistic the interfaces {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.TrainableTransition} and {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.SamplingTransition}
 * can be used to estimate the parameters (cf. {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.TrainableTransition#estimateFromStatistic()})
 * or to draw new parameters (cf. {@link de.jstacs.sampling.SamplingFromStatistic#drawParametersFromStatistic()}).
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.TrainableTransition
 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.SamplingTransition
 */
public interface TransitionWithSufficientStatistic extends Transition {

	/**
	 * This method reset the internal sufficient statistics that can be used for estimating the parameters.
	 */
	public abstract void resetStatistic();

	/**
	 * This method allows to add a certain <code>weight</code> to the sufficient statistic of a specific transition. 
	 * 
	 * @param layer the layer of the matrix
	 * @param index the index encoding the context
	 * @param childIdx the index of the child
	 * @param weight the weight added to the sufficient statistic
	 * @param sequence the sequence
	 * @param sequencePosition the position within the sequence
	 */
	public abstract void addToStatistic( int layer, int index, int childIdx, double weight, Sequence sequence, int sequencePosition );
	
	/**
	 * This method joins the statistics of different instances and sets this joined statistic as statistic of each instance.
	 * 
	 * This method might be used for instance in a multi-threaded optimization to join partial statistics.
	 * 
	 * @param transitions the transitions to be joined
	 */
	public abstract void joinStatistics(Transition... transitions);
	
	/**
         * This method calculates a score for the current statistics, which is independent from the current parameters
         *
         * In general the gamma-score is a product of gamma-functions parameterized with the current statistics
         *
         * @return the logarithm of the gamma-score for the current statistics
         */
	public double getLogGammaScoreFromStatistic();
}