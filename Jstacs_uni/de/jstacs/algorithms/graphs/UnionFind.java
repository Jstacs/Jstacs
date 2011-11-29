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

package de.jstacs.algorithms.graphs;

import java.util.Arrays;

/**
 * This class implements an Union-Find-algorithm with path contraction.
 * 
 * @author Jens Keilwagen
 */
public class UnionFind {

	private int[] a;

	/**
	 * Creates a new Union-Find data structure with <code>size</code> nodes.
	 * 
	 * @param size
	 *            the number of nodes
	 */
	public UnionFind( int size ) {
		a = new int[size];
		Arrays.fill( a, -1 );
	}

	/**
	 * Returns the connected components of the graph.
	 * 
	 * @return the connected components of the graph
	 */
	public int[][] getComponents() {
		int i, j, count1 = 0, count2;
		for( i = 0; i < a.length; i++ ) {
			if( a[i] < 0 ) {
				count1++;
			} else {
				find( i );
			}
		}
		int[][] erg = new int[count1][];
		int[] help;
		count1 = 0;
		i = 0;
		while( i < a.length && count1 < erg.length ) {
			if( a[i] < 0 ) {
				help = new int[a.length];
				help[0] = i;
				count2 = 1;
				for( j = 0; j < a.length; j++ ) {
					if( i == a[j] ) {
						help[count2++] = j;
					}
				}
				erg[count1] = new int[count2];
				for( j = 0; j < count2; j++ ) {
					erg[count1][j] = help[j];
				}
				Arrays.sort( erg[count1] );
				count1++;
			}
			i++;
		}
		return erg;
	}

	/**
	 * Finds the root of the tree with node <code>n</code> and does path
	 * contraction. The implemented path contraction is strict, that means it
	 * hooks every node on the way to the root directly to the root.
	 * 
	 * @param n
	 *            the node
	 * 
	 * @return the root of the tree with node <code>n</code>
	 */
	public int find( int n ) {
		int m, l = n;
		while( a[l] >= 0 ) {
			l = a[l];
		}
		// Now l == find(k)
		while( a[n] > 0 ) {
			m = a[n];
			a[n] = l;
			n = m;
		}
		// Now all nodes on the path point to l
		return l;
	}

	/**
	 * This method unions the tree that includes the node <code>k1</code> with
	 * the tree that includes node <code>k2</code>. If both nodes are already in
	 * the same tree nothing is done.
	 * 
	 * @param k1
	 *            one node
	 * @param k2
	 *            another node
	 * 
	 * @return <code>true</code> if a real union is done, <code>false</code> if
	 *         nothing is done
	 */
	public boolean union( int k1, int k2 ) {
		k1 = find( k1 );
		k2 = find( k2 );
		if( k1 != k2 ) {
			// the nodes are not in the same tree, so a real union has to be
			// done
			// the set of k1 is larger (since a[k1] = -(number of nodes in component k1))
			if( a[k1] < a[k2] ) {
				int i = k1;
				k1 = k2;
				k2 = i;
			}
			a[k2] += a[k1];
			a[k1] = k2;
			return true;
		}
		return false;
	}
}
