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