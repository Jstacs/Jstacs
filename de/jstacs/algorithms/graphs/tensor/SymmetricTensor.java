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
 * This class can be used for {@link Tensor}s with a special symmetry property.
 * This property states that for a given target node and for all permutations of
 * parent nodes the edge weight is constant. This is very helpful to save
 * storage and time in some algorithms. Since it is very useful it also enables
 * you to use this class if your tensor does not fulfill this property. For that
 * case the tensor saves only the maximum edge weight and the permutation of the
 * parent nodes.
 * 
 * @author Jens Keilwagen
 * 
 * @see Tensor
 * @see AsymmetricTensor
 */
public class SymmetricTensor extends Tensor {

	private int[][] pascal;

	private double[][][] tensor;

	private int[][][] trueEdgeName;

	private int[] indices;

	private int idxOfLastBest;

	/**
	 * This constructor creates an empty symmetric tensor with given dimension.
	 * 
	 * @param n
	 *            the number of nodes
	 * @param k
	 *            the order
	 * 
	 * @see Tensor#Tensor(int, byte)
	 */
	public SymmetricTensor( int n, byte k ) {
		super( n, k );
		int j = 1, i = 0, l = k + 1;
		indices = new int[k];

		pascal = new int[L][l];
		for( i = 0; i < L; i++ ) {
			pascal[i][0] = 1;
		}
		for( i = L - 2; i >= 0; i-- ) {
			for( j = 1; j < l; j++ ) {
				pascal[i][j] = pascal[i + 1][j - 1] + pascal[i + 1][j];
			}
		}
		/*
		for( i = 0; i < L; i++ )
		{
			System.out.println( Arrays.toString(pascal[i]) );
		}
		*/
		trueEdgeName = new int[l - 1][][];
		tensor = new double[l][][];
		tensor[0] = new double[n][1];
		for( j = 1; j < l; j++ ) {
			// tensor[j] = new double[n][pascal[L][j]];
			tensor[j] = new double[n][pascal[0][j - 1] + pascal[0][j]];
			trueEdgeName[j - 1] = new int[n][tensor[j][0].length];
			//System.out.println( j + ", " + n + ", " + tensor[j][0].length );
			for( i = 0; i < n; i++ ) {
				Arrays.fill( tensor[j][i], Double.NEGATIVE_INFINITY );
			}
		}
	}

	/**
	 * The constructor can be used creating a new {@link SymmetricTensor} as
	 * weighted sum of {@link SymmetricTensor}s.
	 * 
	 * @param parts
	 *            the {@link SymmetricTensor}s
	 * @param weights
	 *            the weight for each {@link SymmetricTensor}
	 * 
	 * @throws IllegalArgumentException
	 *             if something went wrong
	 * 
	 * @see SymmetricTensor#SymmetricTensor(int, byte)
	 */
	public SymmetricTensor( SymmetricTensor[] parts, double[] weights ) throws IllegalArgumentException {
		this( parts[0].getNumberOfNodes(), parts[0].getOrder() );
		int counter1 = 1, counter2, counter3, counter4, n = parts.length, l = parts[0].getNumberOfNodes(), k = parts[0].getOrder();
		if( n != weights.length ) {
			throw new IllegalArgumentException( "The parts and the weights have to have the same dimension." );
		}

		for( ; counter1 < n; counter1++ ) {
			if( parts[counter1].getNumberOfNodes() != l || parts[counter1].getOrder() != k ) {
				throw new IllegalArgumentException( "All parts have to have the same order and number of nodes." );
			}
		}

		for( counter1 = 0; counter1 <= k; counter1++ ) {
			for( counter2 = 0; counter2 < l; counter2++ ) {
				for( counter3 = 0; counter3 < tensor[counter1][counter2].length; counter3++ ) {
					tensor[counter1][counter2][counter3] = parts[0].tensor[counter1][counter2][counter3] * weights[0];
					if( counter1 > 0 ) {
						trueEdgeName[counter1 - 1][counter2][counter3] = parts[0].trueEdgeName[counter1 - 1][counter2][counter3];
					}
					for( counter4 = 1; counter4 < n; counter4++ ) {
						if( counter1 > 0 && trueEdgeName[counter1 - 1][counter2][counter3] != parts[counter4].trueEdgeName[counter1 - 1][counter2][counter3] ) {
							throw new IllegalArgumentException( "The tensors does not encode the same graph, since at least one edge has a differnet parent permutation." );
						}
						tensor[counter1][counter2][counter3] += parts[counter4].tensor[counter1][counter2][counter3] * weights[counter4];
					}
				}
			}
		}

	}

	/**
	 * This constructor creates and checks a filled asymmetric tensor from an
	 * {@link AsymmetricTensor} instance.
	 * 
	 * @param asym_tensor
	 *            the asymmetric instance
	 * 
	 * @see AsymmetricTensor
	 * @see SymmetricTensor#SymmetricTensor(double[][][], int, byte)
	 */
	public SymmetricTensor( AsymmetricTensor asym_tensor ) {
		this( asym_tensor.tensor, asym_tensor.getNumberOfNodes(), asym_tensor.getOrder() );
	}

	/**
	 * This constructor creates and checks a filled asymmetric tensor with given
	 * dimension.
	 * 
	 * @param asym_tensor
	 *            the tensor weights
	 * @param N
	 *            the number of nodes
	 * @param k
	 *            the order
	 * 
	 * @see SymmetricTensor#SymmetricTensor(int, byte)
	 */
	public SymmetricTensor( double[][][] asym_tensor, int N, byte k ) {
		this( N, k );

		int j = 0, idx, l;
		do {
			setRootValue( j, asym_tensor[0][j][0] );
		} while( ++j < N );
		int[] counter, array;
		int M = N - 1;
		boolean erg;
		for( byte i = 1; i <= order; i++ ) {
			counter = new int[i + 1];
			array = new int[i];
			for( l = i - 1, j = 0; l >= 0; l--, j++ ) {
				counter[l] = j;
			}
			do {
				System.arraycopy( counter, 0, array, 0, i );
				Arrays.sort( array );
				idx = 0;
				for( j = 0; j < i; j++ ) {
					idx *= N;
					idx += counter[j];
				}
				// System.out.println( K_DAG.toString( counter, i ) + " -> " +
				// idx );
				for( j = 0; j < N; j++ ) {
					erg = array[0] != j;
					l = 1;
					while( l < i && ( erg &= ( array[l] != j ) && array[l] != array[l - 1] ) ) {
						l++;
					}
					if( l == i && erg ) {
						setValue( i, asym_tensor[i][j][idx], j, counter );
					}
				}
				j = 0;
				while( counter[j] == M ) {
					counter[j++] = 0;
				}
				counter[j]++;
			} while( counter[i] == 0 );
		}
	}

	/**
	 * Returns the maximal weight which can be reached for an edge from
	 * <code>k</code> nodes from the (encoded) set <code>par</code> to the node
	 * <code>child</code>.
	 * 
	 * @param child
	 *            the target node
	 * @param par
	 *            an encoded set, elements are sorted and for all elements
	 *            <code>e &gt;= child</code> the real value is <code>(e+1)</code>,
	 *            e.g. <code>child = 3</code> and possible parents
	 *            <code> = {5,2,0,6,1}</code> then
	 *            <code>par = {0,1,2,4,5}</code>
	 * @param k
	 *            the number of parents used from par
	 * 
	 * @return the best (=maximal) weight
	 */
	public double getBest( int child, int[] par, byte k ) {
		int i = k, j = 0, b = 0;
		byte l = (byte)par.length;
		int[] current = new int[k], parents = new int[k];
		for( ; b < k; b++ ) {
			current[b] = b;
			parents[b] = par[current[b]];
		}
		double max = Double.NEGATIVE_INFINITY;
		do {
			// j = getIndex( parents, k );
			j = getIndices( parents, j, k );
			if( max < tensor[k][child][j] ) {
				max = tensor[k][child][j];
				idxOfLastBest = j;
			}
			// System.out.println( K_DAG.toString( parents, k ) + " " + j );

			// what's the maximum index at position i
			// - we have only l indices
			// - the positions bigger than i have to be filled
			// -> so at position k-1 => l-1, k-2 => l-2
			b = l;
			while( --i >= 0 && current[i] == --b );
			if( i >= 0 ) {
				// System.out.println( i + " " + current[i] + " ?=? " + (b-i) );
				j = i;
				parents[i] = par[++current[i++]];
				for( ; i < k; i++ ) {
					current[i] = current[i - 1];
					parents[i] = par[++current[i]];
				}
			}
		} while( i == k );
		idxOfLastBest = trueEdgeName[k - 1][child][idxOfLastBest];
		// System.out.println( max );
		return max;
	}

	/**
	 * This method decodes an index in an edge.
	 * 
	 * @param idx
	 *            the encoded index
	 * @param k
	 *            the order
	 * 
	 * @return the edge
	 */
	public int[] getEdgeFromIndex( int idx, int k ) {
		int[] erg = new int[k];
		for( int j = 0; j < k; j++ ) {
			erg[j] = ( idx % powers[1] );
			idx /= powers[1];
		}
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.graphs.tensor.Tensor#getMaximalEdgeFor(byte, int, int[])
	 */
	@Override
	public int[] getMaximalEdgeFor( byte k, int child, int... parents ) {
		return getEdgeFromIndex( k == 0 ? child : trueEdgeName[k - 1][child][getIndex( sortAndClear( child, parents, k ), k )], k + 1 );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.graphs.tensor.Tensor#getRootValue(int)
	 */
	@Override
	public double getRootValue( int child ) {
		return tensor[0][child][0];
	}

	/**
	 * Returns the edge from {@link #getBest(int, int[], byte)} in an encoded
	 * index.
	 * 
	 * @return the edge from {@link #getBest(int, int[], byte)} in an encoded
	 *         index
	 * 
	 * @see SymmetricTensor#getBest(int, int[], byte)
	 */
	public int getTrueIndexForLastGetBest() {
		return idxOfLastBest;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.graphs.tensor.Tensor#getValue(byte, int, int[])
	 */
	@Override
	public double getValue( byte k, int child, int... parents ) {
		return tensor[k][child][getIndex( sortAndClear( child, parents, k ), k )];
	}

	/**
	 * Sets the value for the edge <code>parents -&gt; child</code> to
	 * <code>Double.NEGATIVE_INFTINITY</code>.
	 * 
	 * @param k
	 *            the order
	 * @param child
	 *            the index of the child node
	 * @param parents
	 *            the indices of the parent nodes
	 */
	@Override
	public void resetValue( byte k, int child, int... parents ) {
		int i = getIndex( sortAndClear( child, parents, k ), k );
		tensor[k][child][i] = Double.NEGATIVE_INFINITY;
		trueEdgeName[k - 1][child][i] = 0;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.graphs.tensor.Tensor#setRootValue(int, double)
	 */
	@Override
	public void setRootValue( int child, double val ) {
		tensor[0][child][0] = val;
	}

	/**
	 * Sets the value if it is bigger than the current value and keeps the
	 * parents information.
	 * 
	 * @param k
	 *            the order
	 * @param val
	 *            the value to be set
	 * @param child
	 *            the index of the child node
	 * @param parents
	 *            the indices of the parent nodes
	 */
	@Override
	public void setValue( byte k, double val, int child, int... parents ) {
		int i = getIndex( sortAndClear( child, parents, k ), k );
		if( tensor[k][child][i] < val ) {
			tensor[k][child][i] = val;
			if( k > 0 ) {
				trueEdgeName[k - 1][child][i] = getAsymIndex( child, parents, k ) + powers[k] * child;
			}
		}
	}

	private static int[] sortAndClear( int child, int[] parents, int anz ) {
		int[] cleared_parents = new int[anz];
		int i = 0, stop, last = 0;
		// clear
		do {
			cleared_parents[i] = parents[i];
			if( parents[i] > child ) {
				cleared_parents[i]--;
			}
		} while( ++i < anz );
		// bubblesort
		if( anz > 1 ) {
			// System.out.println( Arrays.toString( cleared_parents ) );
			int help;
			while( last != anz ) {
				stop = last;
				last = anz;
				for( i = anz - 1; i > stop; --i ) {
					if( cleared_parents[i] < cleared_parents[i - 1] ) {
						help = cleared_parents[i];
						cleared_parents[i] = cleared_parents[i - 1];
						cleared_parents[i - 1] = help;
						last = i;
					}
				}
				// System.out.println( Arrays.toString( cleared_parents ) + " "
				// + last );
			}
			// alternative, but slow since anz is normally very small
			//Arrays.sort(cleared_parents);
		}
		return cleared_parents;
	}

	// same as: getIndices( current, 0, l )
	private int getIndex( int[] current, byte l ) {
		int erg = pascal[current[0]][l];
		for( int i = 1; i < l; i++ ) {
			erg += pascal[current[i]][l - i];
		}
		return erg;
	}

	private int getIndices( int[] current, int i, byte l ) {
		if( i == 0 ) {
			indices[0] = pascal[current[0]][l];
		} else {
			indices[i] = indices[i - 1] + pascal[current[i]][l - i];
		}
		while( ++i < l ) {
			indices[i] = indices[i - 1] + pascal[current[i]][l - i];
		}
		return indices[l - 1];
	}
}
