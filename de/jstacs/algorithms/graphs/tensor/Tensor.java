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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;

import de.jstacs.utils.SafeOutputStream;

/**
 * This is the super class for any tensor. The tensor may be symmetric or
 * asymmetric.
 * 
 * @author Jens Keilwagen
 * 
 * @see AsymmetricTensor
 * @see SymmetricTensor
 */
public abstract class Tensor {

	/**
	 * An array containing the powers for the number of nodes. This might be
	 * used in the computation of some indices.
	 */
	protected int[] powers;

	/**
	 * The number of nodes minus 1.
	 */
	protected int L;

	/**
	 * The order of the tensor.
	 */
	protected byte order;

	/**
	 * Creates a new {@link Tensor} for <code>n</code> nodes and order
	 * <code>k</code>.
	 * 
	 * @param n
	 *            the number of nodes
	 * @param k
	 *            the order
	 * 
	 * @throws IllegalArgumentException
	 *             if n &lt; 0 or k &lt; 1
	 */
	public Tensor( int n, byte k ) throws IllegalArgumentException {
		if( n < 0 ) {
			throw new IllegalArgumentException( "The number of nodes n has to be positive." );
		}
		if( k < 1 ) {
			throw new IllegalArgumentException( "The order k has to be at least 1." );
		}
		L = ( n - 1 );
		order = k;

		powers = new int[k + 2];
		powers[0] = 1;
		for( int j = 1; j <= k + 1; j++ ) {
			powers[j] = powers[j - 1] * n;
		}
	}

	/**
	 * Returns the edge
	 * <code>permute(parents[0],...,parents[k-1]) -> child</code> that maximizes
	 * the score.
	 * 
	 * @param k
	 *            the number of parents to be used
	 * @param child
	 *            the child node
	 * @param parents
	 *            the parent nodes (only the first <code>k</code> will be used)
	 * 
	 * @return the edge
	 *         <code>permute(parents[0],...,parents[k-1]) -> child</code> that
	 *         maximizes the score
	 */
	public abstract int[] getMaximalEdgeFor( byte k, int child, int... parents );

	/**
	 * Returns the number of nodes.
	 * 
	 * @return the number of nodes
	 */
	public int getNumberOfNodes() {
		return( L + 1 );
	}

	/**
	 * Returns the order.
	 * 
	 * @return the order.
	 */
	public byte getOrder() {
		return order;
	}

	/**
	 * Returns the value for <code>child</code> as root.
	 * 
	 * @param child
	 *            the name of the node
	 * 
	 * @return the value for the node
	 */
	public abstract double getRootValue( int child );

	/**
	 * Returns the value for the edge
	 * <code>parents[0],...,parents[k-1] -&gt; child</code>.
	 * 
	 * @param k
	 *            the number of parents to be used
	 * @param child
	 *            the child node
	 * @param parents
	 *            the parent nodes (only the first <code>k</code> will be used)
	 * 
	 * @return the value for the edge
	 *         <code>parents[0],...,parents[k-1] -&gt; child</code>
	 */
	public abstract double getValue( byte k, int child, int... parents );

	/**
	 * Sets the value <code>val</code> for the root node <code>child</code>.
	 * 
	 * @param child
	 *            the name of the node
	 * @param val
	 *            the value for the node
	 */
	public abstract void setRootValue( int child, double val );

	/**
	 * Sets the value for the edge
	 * <code>parents[0],...,parents[k-1] -&gt; child</code>.
	 * 
	 * @param k
	 *            the number of parents to be used
	 * @param val
	 *            the new value
	 * @param child
	 *            the child node
	 * @param parents
	 *            the parent nodes (only the first <code>k</code> will be used)
	 */
	public abstract void setValue( byte k, double val, int child, int... parents );

	/**
	 * Sets the value for the edge
	 * <code>parents[0],...,parents[k-1] -&gt; child</code> to
	 * <code>Double.NEGATIVE_INFINITY</code>.
	 * 
	 * @param k
	 *            the number of parents to be used
	 * @param child
	 *            the child node
	 * @param parents
	 *            the parent nodes (only the first <code>k</code> will be used)
	 */
	public abstract void resetValue( byte k, int child, int... parents );

	/**
	 * This method writes a {@link Tensor} in the exchange format in a specified
	 * file.
	 * 
	 * @param fname
	 *            the file name
	 * @param desc
	 *            gives you the possibility to write the description of the
	 *            tensor file
	 * 
	 * @throws IOException
	 *             if something went wrong with the file
	 * 
	 * @see Tensor#readTensorFromFile(String, boolean)
	 */
	public void writeTensorToFile( String fname, OutputStream desc ) throws IOException {
		int i, l, phase = getNumberOfNodes(), j, help, anz;
		byte k = 1;
		int[] parents;
		BufferedWriter writer = new java.io.BufferedWriter( new java.io.FileWriter( fname ) );
		writer.write( phase + "" );
		writer.newLine();
		SafeOutputStream description = SafeOutputStream.getSafeOutputStream( desc );
		description.writeln( "number of nodes" );

		writer.write( "E 0" );
		writer.newLine();
		description.writeln( "E 0" );
		for( i = 0; i < phase; i++ ) {
			writer.write( "" + getRootValue( i ) );
			writer.newLine();
			description.writeln( "v(" + i + ") = value for node " + i );
		}

		while( k <= order ) {
			writer.write( "E " + k );
			writer.newLine();
			description.writeln( "E " + k );
			parents = new int[k];
			anz = (int)Math.pow( phase, k );
			for( i = 0; i < phase; i++ ) {
				for( j = 0; j < anz; j++ ) {
					help = j;
					for( l = 0; l < k; l++ ) {
						parents[l] = help % phase;
						help /= phase;
					}
					if( check( i, parents.clone() ) ) {
						writer.write( "" + getValue( k, i, parents ) );
						writer.newLine();
						description.writeln( "v(" + i
												+ ", "
												+ Arrays.toString( parents )
												+ ") = value for edge "
												+ i
												+ " <- "
												+ toString( parents ) );
					}
				}
			}
			k++;
		}
		description.flush();
		description.close();
		writer.flush();
		writer.close();
	}

	private String toString( int[] array ) {
		if( array != null && array.length > 0 ) {
			String res = "";
			int i = 0, l = array.length - 1;
			for( ; i < l; i++ ) {
				res += array[i] + "<-";
			}
			return "[" + res + array[i] + "]";
		} else {
			return "[]";
		}
	}

	/**
	 * The opposite of the method
	 * {@link #writeTensorToFile(String, OutputStream)}.
	 * 
	 * @param fname
	 *            the file name
	 * @param asym
	 *            <code>true</code> generates an {@link AsymmetricTensor},
	 *            otherwise a {@link SymmetricTensor}
	 * 
	 * @return the tensor from the file
	 * 
	 * @throws IOException
	 *             if something went wrong with the file
	 * @throws NumberFormatException
	 *             if the file could not be parsed correctly
	 * 
	 * @see Tensor#writeTensorToFile(String, OutputStream)
	 */
	public static Tensor readTensorFromFile( String fname, boolean asym ) throws NumberFormatException, IOException {
		BufferedReader reader = new java.io.BufferedReader( new java.io.FileReader( fname ) );
		int i, n, l, help, anz = 1, parents = 0, phase = Integer.parseInt( reader.readLine() );
		ArrayList<double[][]> list = new ArrayList<double[][]>();
		double[][] values;
		int[] par;
		String s;
		while( ( s = reader.readLine() ) != null && s.equals( "E " + parents ) ) {
			values = new double[phase][anz];
			if( parents == 0 ) {
				for( n = 0; n < phase; n++ ) {
					values[n][0] = Double.parseDouble( reader.readLine() );
				}
			} else {
				par = new int[parents];
				for( n = 0; n < phase; n++ ) {
					for( i = 0; i < anz; i++ ) {
						help = i;
						for( l = 0; l < parents; l++ ) {
							par[l] = help % phase;
							help /= phase;
						}
						if( check( n, par.clone() ) ) {
							values[n][i] = Double.parseDouble( reader.readLine() );
						} else {
							values[n][i] = NOT_USED;
						}
					}
				}
			}
			anz *= phase;
			parents++;
			list.add( values );
		}
		reader.close();
		double[][][] ten = list.toArray( new double[0][0][0] );
		if( asym ) {
			return new AsymmetricTensor( ten, phase, (byte)( ten.length - 1 ) );
		} else {
			return new SymmetricTensor( ten, phase, (byte)( ten.length - 1 ) );
		}
	}

	private static final double NOT_USED = Double.NaN;

	/**
	 * Creates a three-dimensional <code>double</code> array representation of
	 * the {@link Tensor}.
	 * 
	 * @return a three-dimensional <code>double</code> array representation of
	 *         the {@link Tensor}
	 */
	public double[][][] toDouble3DArray() {
		int i, l, phase = getNumberOfNodes(), j, help, anz;
		byte k = 1;
		int[] parents;
		double[][][] erg = new double[order + 1][phase][];
		for( i = 0; i < phase; i++ ) {
			erg[0][i] = new double[]{ getRootValue( i ) };
		}

		while( k <= order ) {
			parents = new int[k];
			anz = (int)Math.pow( phase, k );
			for( i = 0; i < phase; i++ ) {
				erg[k][i] = new double[anz];
				for( j = 0; j < anz; j++ ) {
					help = j;
					for( l = 0; l < k; l++ ) {
						parents[l] = help % phase;
						help /= phase;
					}
					if( check( i, parents.clone() ) ) {
						erg[k][i][j] = getValue( k, i, parents );
					} else {
						erg[k][i][j] = NOT_USED;
					}
				}
			}
			k++;
		}
		return erg;
	}

	private static boolean check( int current, int[] array ) {
		Arrays.sort( array );
		boolean erg = array[0] != current;
		int i = 1;
		while( i < array.length && ( erg &= ( array[i] != current ) && array[i] != array[i - 1] ) ) {
			i++;
		}
		return erg;
	}

	/**
	 * Returns the index for an asymmetric tensor.
	 * 
	 * @param child
	 *            the child node
	 * @param parents
	 *            the parent notes
	 * @param k
	 *            the order
	 * 
	 * @return the index for an asymmetric tensor
	 * 
	 * @see Tensor#powers
	 */
	protected int getAsymIndex( int child, int[] parents, byte k ) {
		int erg = 0;
		for( int j = 0; j < k; j++ ) {
			erg += powers[j] * parents[j];
		}
		return erg;
	}
}
