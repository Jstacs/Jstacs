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

package de.jstacs.classifier.scoringFunctionBased;

import de.jstacs.algorithms.optimization.DifferentiableFunction;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.data.DataSet;

/**
 * This is the main function for the {@link ScoreClassifier}.
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.algorithms.optimization.Optimizer
 */
public abstract class OptimizableFunction extends DifferentiableFunction {

	/**
	 * This <code>enum</code> defines the kinds of parameters that can be
	 * returned by the method
	 * {@link OptimizableFunction#getParameters(KindOfParameter)}.
	 * 
	 * @author Jens Keilwagen
	 * 
	 * @see OptimizableFunction#getParameters(KindOfParameter)
	 */
	public enum KindOfParameter {
		/**
		 * Indicates that all parameters will have value 0.
		 */
		ZEROS,
		/**
		 * Indicates that the last parameters will be (re-)used.
		 */
		LAST,
		/**
		 * Indicates that some plug-in parameters will be used (e.g.
		 * MAP-parameters).
		 */
		PLUGIN
	}

	/**
	 * Returns some parameters that can be used for instance as start
	 * parameters.
	 * 
	 * @param kind
	 *            the kind of the class parameters that will be returned
	 * 
	 * @return some start parameters
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public abstract double[] getParameters( KindOfParameter kind ) throws Exception;

	/**
	 * Sets the current values as parameters.
	 * 
	 * @param current
	 *            the current values
	 * 
	 * @throws DimensionException
	 *             if the dimension of the current values does not match with
	 *             the dimension of the internal parameters
	 */
	public abstract void setParams( double[] current ) throws DimensionException;

	/**
	 * Resets the all objects and pre-computed values.
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public abstract void reset() throws Exception;

	/**
	 * Returns the data for each class used in this {@link OptimizableFunction}.
	 * 
	 * @return the data for each class
	 * 
	 * @see OptimizableFunction#getSequenceWeights()
	 */
	public abstract DataSet[] getData();

	/**
	 * Returns the weights for each {@link de.jstacs.data.Sequence} for each
	 * class used in this {@link OptimizableFunction}.
	 * 
	 * @return the weights for each {@link de.jstacs.data.Sequence} and each
	 *         class
	 * 
	 * @see OptimizableFunction#getData()
	 */
	public abstract double[][] getSequenceWeights();

	/**
	 * This method sets the data set and the sequence weights to be used. It also allows to do further preparation for the computation on this data.
	 * 
	 * @param data the data sets
	 * @param weights the sequence weights for each sequence in each data set
	 * 
	 * @throws IllegalArgumentException if the data or the weights can not be used
	 */
	public abstract void setDataAndWeights( DataSet[] data, double[][] weights ) throws IllegalArgumentException;
}
