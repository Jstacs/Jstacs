package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions;

import java.text.NumberFormat;

import javax.naming.OperationNotSupportedException;

import de.jstacs.Storable;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.Transition;

/**
 * This interface declares all method for an emission of a state.
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.State
 */
public interface Emission extends Storable, Cloneable {

	/**
	 * This method returns the {@link AlphabetContainer} of this emission.
	 * 
	 * @return the {@link AlphabetContainer} of this emission
	 */
	public AlphabetContainer getAlphabetContainer();

	/**
	 * This method initializes the emission randomly.
	 */
	public void initializeFunctionRandomly();

	/**
	 * This method computes the logarithm of the likelihood.
	 * 
	 * @param forward whether to use the forward or the reverse strand
	 * @param startPos the start position
	 * @param endPos the end position
	 * @param seq the sequence
	 * 
	 * @return the logarithm of the probability
	 * 
	 * @throws OperationNotSupportedException if <code>forward=false</code> and the reverse complement of the sequence <code>seq</code> is not defined  
	 */
	public double getLogProbFor( boolean forward, int startPos, int endPos, Sequence seq ) throws OperationNotSupportedException;

	/**
	 * Returns a value that is proportional to the log of the prior. For maximum likelihood (ML) 0
	 * should be returned.
	 * 
	 * @return a value that is proportional to the log of the prior
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorTerm();

	/**
	 * This method resets the internal sufficient statistic.
	 */
	public void resetStatistic();

	/**
	 * This method adds the <code>weight</code> to the internal sufficient statistic.
	 * 
	 * @param forward whether to use the forward or the reverse strand
	 * @param startPos the start position
	 * @param endPos the end position
	 * @param weight the weight of the sequence
	 * @param seq the sequence
	 *  
	 * @throws OperationNotSupportedException if <code>forward=false</code> and the reverse complement of the sequence <code>seq</code> is not defined 
	 */
	public void addToStatistic( boolean forward, int startPos, int endPos, double weight, Sequence seq ) throws OperationNotSupportedException;

	/**
	 * This method joins the statistics of different instances and sets this joined statistic as statistic of each instance.
	 * 
	 * This method might be used for instance in a multi-threaded optimization to join partial statistics.
	 * 
	 * @param emissions the emissions to be joined
	 */
	public void joinStatistics(Emission... emissions);
	
	/**
	 * This method estimates the parameters from the internal sufficient statistic.
	 */
	public void estimateFromStatistic();
	
	/**
	 * Returns the graphviz string for the shape of the node.
	 * @param forward if this emission is used on the forward strand
	 * @return the shape
	 */
	public String getNodeShape(boolean forward);
	
	
	/**
	 * Returns the graphviz label of the node containing this emission.
	 * @param weight the weight of the node which is represented by 
	 * 			the color of the node, or -1 for no representation, i.e.,
	 * 			white background
	 * @param name the name of the state using this emission
	 * @param nf the {@link NumberFormat} for formatting the textual representation of this emission
	 * @return the label
	 */
	public String getNodeLabel(double weight, String name, NumberFormat nf);

	/**
	 * Set values of parameters of the instance to the value of the parameters of the given instance.
	 * It can be assumed that the given instance and the current instance are from the same class.
	 * 
	 * This method might be used for instance in a multi-threaded optimization to broadcast the parameters. 
	 * 
	 * @param t the emission with the parameters to be set 
	 * 
	 * @throws IllegalArgumentException if the assumption about the same class for given and current instance is wrong
	 */
	public void setParameters( Emission t ) throws IllegalArgumentException;
}