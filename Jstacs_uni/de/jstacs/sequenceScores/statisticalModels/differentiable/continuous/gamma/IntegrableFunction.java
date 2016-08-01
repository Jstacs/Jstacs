package de.jstacs.sequenceScores.statisticalModels.differentiable.continuous.gamma;

import de.jstacs.algorithms.optimization.Function;

/**
 * Interface for a function that can be integrated using the class {@link NumericalIntegration}.
 * @author Jan Grau
 *
 */
public interface IntegrableFunction extends Function {

	/**
	 * Returns the limits of the domain of this function. The value in <code>limits[i][0]</code> is 
	 * the lower limit and <code>limits[i][1] is the upper limit of dimension <code>i</code>.
	 * 
	 * @return the limits
	 */
	public abstract double[][] getLimits();
	
	/**
	 * Returns <code>true</code> if this function is unimodal.
	 * @return if this function is unimodal
	 */
	public abstract boolean isUnimodal();
	
}
