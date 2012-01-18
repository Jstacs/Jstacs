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

package de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix;


/**
 * This enum can be used to obtain the weights for well-known optimization tasks.
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LogGenDisMixFunction
 */
public enum LearningPrinciple
{
	/**
	 * Maximum Likelihood.
	 */
	ML,
	/**
	 * Maximum a posteriori.
	 */
	MAP,
	/**
	 * Maximum conditional likelihood.
	 */
	MCL,
	/**
	 * Maximum supervised posterior.
	 */
	MSP;
	
	/**
	 * This constant is the array index of the weighting factor for the conditional likelihood.
	 * @see LearningPrinciple#getBeta(LearningPrinciple)
	 */
	public static final int CONDITIONAL_LIKELIHOOD_INDEX = 0;

	/**
	 * This constant is the array index of the weighting factor for the likelihood.
	 * @see LearningPrinciple#getBeta(LearningPrinciple)
	 */
	public static final int LIKELIHOOD_INDEX = 1;
	
	/**
	 * This constant is the array index of the weighting factor for the prior.
	 * @see LearningPrinciple#getBeta(LearningPrinciple)
	 */
	public static final int PRIOR_INDEX = 2;
	
	/**
	 * This method returns the standard weights for a predefined key.
	 *  
	 * @param key the key
	 * 
	 * @return the weights array
	 * 
	 * @see LearningPrinciple#ML
	 * @see LearningPrinciple#MAP
	 * @see LearningPrinciple#MCL
	 * @see LearningPrinciple#MSP
	 */
	public static double[] getBeta( LearningPrinciple key )
	{
		double[] beta = new double[3];
		switch( key )
		{
			case ML:
				beta[LIKELIHOOD_INDEX] = 1;
				break;
			case MAP:
				beta[LIKELIHOOD_INDEX] = beta[PRIOR_INDEX] = 0.5;
				break;
			case MCL:
				beta[CONDITIONAL_LIKELIHOOD_INDEX] = 1;
				break;
			case MSP:
				beta[CONDITIONAL_LIKELIHOOD_INDEX] = beta[PRIOR_INDEX] = 0.5;
				break;
			default:
				throw new IllegalArgumentException( "Unknown key" );
		}
		return beta;
	}
	/**
	 * This method checks the values of the <code>weights</code> array. If everything is okay it returns a deep copy
	 * of the array, otherwise it throws an exception
	 * 
	 * @param weights and array of length 3 with non-negative entries that sum to 1
	 * 
	 * @return a deep copy of the array
	 * 
	 * @throws IllegalArgumentException if the weights array is not correct
	 */
	public static double[] checkWeights( double[] weights ) throws IllegalArgumentException
	{
		if( weights.length != 3 )
		{
			throw new IllegalArgumentException( "Wrong dimension for the weights." );
		}
		double sum = weights[0] + weights[1] + weights[2];
		if( Math.abs( 1d - sum ) > 1E-9 )
		{
			throw new IllegalArgumentException( "The weights have to be normalized to 1." );
		}
		return new double[]{ weights[0], weights[1], weights[2] };
	}
}