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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.differentiableStatisticalModels;

import de.jstacs.StatisticalModel;
import de.jstacs.differentiableSequenceScores.DifferentiableSequenceScore;

/**
 * The interface for normalizable {@link DifferentiableSequenceScore}s.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public interface DifferentiableStatisticalModel extends DifferentiableSequenceScore, StatisticalModel {
	/**
	 * Returns the size of the event space of the random variables that are
	 * affected by parameter no. <code>index</code>, i.e. the product of the
	 * sizes of the alphabets at the position of each random variable affected
	 * by parameter <code>index</code>. For DNA alphabets this corresponds to 4
	 * for a PWM, 16 for a WAM except position 0, ...
	 * 
	 * @param index
	 *            the index of the parameter
	 * 
	 * @return the size of the event space
	 */
	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index );

	/**
	 * Returns the logarithm of the sum of the scores over all sequences of the event space.
	 * 
	 * @return the logarithm of the normalization constant Z
	 */
	public double getLogNormalizationConstant();

	/**
	 * Returns the logarithm of the partial normalization constant for the parameter with index
	 * <code>parameterIndex</code>. This is the logarithm of the partial derivation of the
	 * normalization constant for the parameter with index
	 * <code>parameterIndex</code>, 
	 * {@latex.ilb \\[\\log \\frac{\\partial Z(\\underline{\\lambda})}{\\partial \\lambda_{parameterindex}}\\]}.
	 * 
	 * @param parameterIndex
	 *            the index of the parameter
	 * 
	 * @return the logarithm of the partial normalization constant
	 * 
	 * @throws Exception
	 *             if something went wrong with the normalization
	 *             
	 * @see DifferentiableStatisticalModel#getLogNormalizationConstant()
	 */
	public double getLogPartialNormalizationConstant(int parameterIndex)
			throws Exception;

	/**
	 * This method computes a value that is proportional to
	 * 
	 * <p><code>
	 * {@link #getESS()} * {@link #getLogNormalizationConstant()} + Math.log( prior )
	 * </code></p>
	 * 
	 * where <code>prior</code> is the prior for the parameters of this model.
	 * 
	 * @return a value that is proportional to <code>{@link #getESS()} * {@link #getLogNormalizationConstant()} + Math.log( prior ).</code>
	 * 
	 * @see DifferentiableStatisticalModel#getESS()
	 * @see DifferentiableStatisticalModel#getLogNormalizationConstant()
	 */
	public double getLogPriorTerm();

	/**
	 * This method computes the gradient of {@link #getLogPriorTerm()} for each
	 * parameter of this model. The results are added to the array
	 * <code>grad</code> beginning at index <code>start</code>.
	 * 
	 * @param grad
	 *            the array of gradients
	 * @param start
	 *            the start index in the <code>grad</code> array, where the
	 *            partial derivations for the parameters of this models shall be
	 *            entered
	 * 
	 * @throws Exception
	 *             if something went wrong with the computing of the gradients
	 * 
	 * @see DifferentiableStatisticalModel#getLogPriorTerm()
	 */
	public void addGradientOfLogPriorTerm(double[] grad, int start)
			throws Exception;

	/**
	 * This method indicates whether the implemented score is already normalized
	 * to 1 or not. The standard implementation returns <code>false</code>.
	 * 
	 * @return <code>true</code> if the implemented score is already normalized
	 *         to 1, <code>false</code> otherwise
	 */
	public boolean isNormalized();
	
	/**
	 * Returns the equivalent sample size (ess) of this model, i.e. the
	 * equivalent sample size for the class or component that is represented by
	 * this model.
	 * 
	 * @return the equivalent sample size.
	 */
	public double getESS();
}
