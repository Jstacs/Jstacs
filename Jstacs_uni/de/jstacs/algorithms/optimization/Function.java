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

package de.jstacs.algorithms.optimization;

/**
 * This interface is the framework for any mathematical function
 * <code>f: R^n -> R</code>.
 * 
 * @author Jens Keilwagen
 */
public interface Function {

	/**
	 * Evaluates the function at a certain vector (in mathematical sense)
	 * <code>x</code>.
	 * 
	 * @param x
	 *            the current vector
	 * 
	 * @return the evaluation of the function
	 * 
	 * @throws DimensionException
	 *             if <code>dim(x) != n</code>, with <code>f: R^n -> R</code>
	 * @throws EvaluationException
	 *             if there was something wrong during the evaluation of the
	 *             function
	 */
	public double evaluateFunction( double[] x ) throws DimensionException, EvaluationException;

	/**
	 * Returns the dimension of the scope of the {@link Function}.
	 * 
	 * @return the dimension of the scope: <code>n</code> with
	 *         <code>f: R^n -> R</code>
	 */
	public int getDimensionOfScope();
}
