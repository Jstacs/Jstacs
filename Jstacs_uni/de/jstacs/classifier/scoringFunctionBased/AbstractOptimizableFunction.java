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

import de.jstacs.data.Sample;

/**
 * This class extends {@link OptimizableFunction} and implements some common
 * methods.
 * 
 * @author Jens Keilwagen
 */
public abstract class AbstractOptimizableFunction extends OptimizableFunction {

	/**
	 * The data that is used to evaluate this function.
	 */
	protected Sample[] data;

	/**
	 * The weights for the data.
	 * 
	 * @see AbstractOptimizableFunction#data
	 */
	protected double[][] weights;

	/**
	 * The class parameters.
	 */
	protected double[] clazz;

	/**
	 * The logarithm of the class parameters.
	 * 
	 * @see AbstractOptimizableFunction#clazz
	 */
	protected double[] logClazz;

	/**
	 * The sums of the weighted data per class and additional the total weight
	 * sum.
	 * 
	 * @see AbstractOptimizableFunction#data
	 * @see AbstractOptimizableFunction#weights
	 */
	protected double[] sum;

	/**
	 * The number of different classes.
	 */
	protected int cl;

	/**
	 * Indicates whether a normalization should be done or not.
	 */
	protected boolean norm;

	/**
	 * Indicates whether only the free parameters or all should be used.
	 */
	protected boolean freeParams;

	/**
	 * The constructor creates an instance using the given weighted data.
	 * 
	 * @param data
	 *            the data
	 * @param weights
	 *            the weights
	 * @param norm
	 *            the switch for using the normalization (division by the number
	 *            of sequences)
	 * @param freeParams
	 *            the switch for using only the free parameters
	 * 
	 * @throws IllegalArgumentException
	 *             if the number of classes or the dimension of the weights is not correct
	 */
	protected AbstractOptimizableFunction( Sample[] data, double[][] weights, boolean norm, boolean freeParams ) throws IllegalArgumentException {
		this.norm = norm;
		this.freeParams = freeParams;
		cl = weights.length;
		logClazz = new double[cl];
		clazz = new double[cl];
		sum = new double[cl + 1];
		setDataAndWeights(data, weights);
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.OptimizableFunction#setDataAndWeights(de.jstacs.data.Sample[], double[][])
	 */
	public void setDataAndWeights( Sample[] data, double[][] weights ) throws IllegalArgumentException {
		if( data.length != cl || weights == null || weights.length != cl ) {
			throw new IllegalArgumentException( "The dimension of the sample or weights (array) is not correct."  );
		}
		this.data = data;
		this.weights = weights;
		sum[cl] = 0;
		int i = 0, j;
		for( ; i < cl; i++ ) {
			sum[i] = 0;
			if( data[i].getNumberOfElements() != weights[i].length ) {
				throw new IllegalArgumentException( "The dimension of the " + i + "-th weights (array) is not correct."  );
			}
			for( j = 0; j < weights[i].length; j++ ) {
				sum[i] += weights[i][j];
			}
			sum[cl] += sum[i];
		}
	}

	/**
	 * This method enables the user to get the parameters without creating a new
	 * array.
	 * 
	 * @param kind
	 *            the kind of the class parameters to be returned in
	 *            <code>erg</code>
	 * @param erg
	 *            the array for the start parameters
	 * 
	 * @throws Exception
	 *             if the array is <code>null</code> or does not have the
	 *             correct length
	 * 
	 * @see OptimizableFunction#getParameters(KindOfParameter)
	 */
	public abstract void getParameters( KindOfParameter kind, double[] erg ) throws Exception;

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.OptimizableFunction#getParameters(KindOfParameter)
	 */
	@Override
	public final double[] getParameters( KindOfParameter kind ) throws Exception {
		double[] temp = new double[getDimensionOfScope()];
		getParameters( kind, temp );
		return temp;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.OptimizableFunction#getData()
	 */
	@Override
	public Sample[] getData() {
		return data.clone();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.scoringFunctionBased.OptimizableFunction#getSequenceWeights()
	 */
	@Override
	public double[][] getSequenceWeights() {
		double[][] res = new double[cl][];
		for( int i = 0; i < cl; i++ ) {
			res[i] = weights[i].clone();
		}
		return res;
	}
}
