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
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.trainable;

import de.jstacs.data.DataSet;
import de.jstacs.sequenceScores.StatisticalModel;

/**
 * This interface defines all methods for a probabilistic model.
 * 
 * @author Andre Gohr, Jan Grau, Jens Keilwagen
 */
public interface TrainableStatisticalModel extends StatisticalModel {

	/**
	 * Creates a clone (deep copy) of the current {@link TrainableStatisticalModel} instance.
	 * 
	 * @return the cloned instance
	 * 
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 */
	public TrainableStatisticalModel clone() throws CloneNotSupportedException;

	/**
	 * Trains the {@link TrainableStatisticalModel} object given the data as {@link DataSet}. <br>
	 * This method should work non-incrementally. That means the result of the
	 * following series: <code>train(data1)</code>; <code>train(data2)</code>
	 * should be a fully trained model over <code>data2</code> and not over
	 * <code>data1+data2</code>. All parameters of the model were given by the
	 * call of the constructor.
	 * 
	 * @param data
	 *            the given sequences as {@link DataSet}
	 * @throws Exception
	 *             if the training did not succeed
	 * 
	 * @see DataSet#getElementAt(int)
	 * @see de.jstacs.data.DataSet.ElementEnumerator
	 */
	public void train(DataSet data) throws Exception;

	/**
	 * Trains the {@link TrainableStatisticalModel} object given the data as {@link DataSet} using
	 * the specified weights. The weight at position i belongs to the element at
	 * position i. So the array <code>weight</code> should have the number of
	 * sequences in the sample as dimension. (Optionally it is possible to use
	 * <code>weight == null</code> if all weights have the value one.)<br>
	 * This method should work non-incrementally. That means the result of the
	 * following series: <code>train(data1)</code>; <code>train(data2)</code>
	 * should be a fully trained model over <code>data2</code> and not over
	 * <code>data1+data2</code>. All parameters of the model were given by the
	 * call of the constructor.
	 * 
	 * @param data
	 *            the given sequences as {@link DataSet}
	 * @param weights
	 *            the weights of the elements, each weight should be
	 *            non-negative
	 * @throws Exception
	 *             if the training did not succeed (e.g. the dimension of
	 *             <code>weights</code> and the number of sequences in the
	 *             sample do not match)
	 * 
	 * @see DataSet#getElementAt(int)
	 * @see de.jstacs.data.DataSet.ElementEnumerator
	 */
	public void train(DataSet data, double[] weights) throws Exception;

	/**
	 * Should give a simple representation (text) of the model as {@link String}.
	 * 
	 * @return the representation as {@link String}
	 */
	public String toString();
}
