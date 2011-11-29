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

package de.jstacs.algorithms.graphs.tensor;

import java.util.Arrays;

/**
 * This class can be used for {@link Tensor}s which are not symmetric, as
 * opposed to the symmetry defined in {@link SymmetricTensor}.
 * 
 * @author Jens Keilwagen
 * 
 * @see Tensor
 * @see SymmetricTensor
 */
public class AsymmetricTensor extends Tensor {

	/**
	 * The internal tensor.
	 */
	protected double[][][] tensor;

	/**
	 * This constructor creates an empty asymmetric tensor with given dimension.
	 * 
	 * @param n
	 *            the number of nodes
	 * @param k
	 *            the order
	 * 
	 * @see Tensor#Tensor(int, byte)
	 */
	public AsymmetricTensor( int n, byte k ) {
		super( n, k );
		tensor = new double[k + 1][n][];
		int i = 0, j, l;
		do {
			l = (int)Math.pow( n, i );
			for( j = 0; j < n; j++ ) {
				tensor[i][j] = new double[l];
				Arrays.fill( tensor[i][j], Double.NEGATIVE_INFINITY );
			}
		} while( ++i <= k );
	}

	/**
	 * This constructor creates and checks a filled asymmetric tensor with given
	 * dimension.
	 * 
	 * @param asym_tensor
	 *            the tensor weights
	 * @param n
	 *            the number of nodes
	 * @param k
	 *            the order
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>n</code> &lt; 0 or <code>k</code> &lt; 1 or the
	 *             given tensor has a wrong dimension
	 * 
	 * @see Tensor#Tensor(int, byte)
	 */
	public AsymmetricTensor( double[][][] asym_tensor, int n, byte k ) throws IllegalArgumentException {
		super( n, k );
		tensor = new double[k + 1][n][];
		int i, l;
		for( i = 0; i <= k; i++ ) {
			for( l = 0; l < n; l++ ) {
				if( asym_tensor[i][l].length == powers[i] ) {
					tensor[i][l] = (double[])asym_tensor[i][l].clone();
				} else {
					throw new IllegalArgumentException( "The given tensor has not the correct dimension." );
				}
			}
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.graphs.tensor.Tensor#getMaximalEdgeFor(byte, int, int[])
	 */
	@Override
	public int[] getMaximalEdgeFor( byte k, int child, int... parents ) {
		int[] erg = new int[k + 1];
		System.arraycopy( parents, 0, erg, 0, k );
		erg[k] = child;
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.graphs.tensor.Tensor#getRootValue(int)
	 */
	@Override
	public double getRootValue( int child ) {
		return tensor[0][child][0];
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.graphs.tensor.Tensor#getValue(byte, int, int[])
	 */
	@Override
	public double getValue( byte k, int child, int... parents ) {
		return tensor[k][child][getAsymIndex( child, parents, k )];
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.graphs.tensor.Tensor#setRootValue(int, double)
	 */
	@Override
	public void setRootValue( int child, double val ) {
		tensor[0][child][0] = val;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.graphs.tensor.Tensor#setValue(byte, double, int, int[])
	 */
	@Override
	public void setValue( byte k, double val, int child, int... parents ) {
		tensor[k][child][getAsymIndex( child, parents, k )] = val;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.graphs.tensor.Tensor#resetValue(byte, int, int[])
	 */
	@Override
	public void resetValue( byte k, int child, int... parents ) {
		tensor[k][child][getAsymIndex( child, parents, k )] = Double.NEGATIVE_INFINITY;
	}
}
