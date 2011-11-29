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

import java.util.ArrayList;
import java.util.Arrays;

/**
 * This class enables you to compute the maximal spanning forest for an
 * undirected, weighted graph. If a weight is
 * <code>Double.NEGATIVE_INFINITY</code> the edge will not be used.
 * 
 * @author Jens Keilwagen
 */
public class MST {

	/**
	 * Does Kruskals algorithm and finds the <b>m</b>aximal <b>s</b>panning
	 * <b>t</b>ree (MST).
	 * 
	 * @param weights
	 *            the matrix of weights, <code>weights.length</code> is the
	 *            number of nodes in the tree, <code>weights[i][j]</code> is the
	 *            weight for edge <code>(i,i+j)</code>
	 * 
	 * @return the MST of the weighted graph
	 */
	public static int[][] kruskal( double[][] weights ) {
		int iterator, n = weights.length;
		int counter[] = { 0, 0 };
		double help;
		ArrayList<Edge> list = new ArrayList<Edge>();
		do {
			iterator = n - 1 - counter[0];
			for( counter[1] = 0; counter[1] < iterator; counter[1]++ ) {
				if( Double.NEGATIVE_INFINITY != ( help = weights[counter[0]][counter[1]] ) ) {
					list.add( new Edge( counter[0], counter[0] + counter[1] + 1, help ) );
				}
			}
		} while( ++counter[0] < n - 1 );
		Edge[] e = new Edge[list.size()];
		list.toArray( e );
		Arrays.sort( e );
		UnionFind uf = new UnionFind( n );
		boolean[] in = new boolean[e.length];
		Arrays.fill( in, false );
		counter[0] = 0;
		for( iterator = e.length - 1; iterator >= 0; iterator-- ) {
			if( uf.union( e[iterator].getStartNode(), e[iterator].getEndNode() ) ) {
				in[iterator] = true;
				counter[0]++;
				//System.out.println( "(" + e[iterator].source + ", "+ e[iterator].target + ") \t" + e[iterator].weight );
			}
		}
		int[][] res = new int[counter[0]][2];
		counter[0] = 0;
		for( iterator = 0; iterator < e.length; iterator++ ) {
			if( in[iterator] ) {
				res[counter[0]][0] = e[iterator].getStartNode();
				res[counter[0]++][1] = e[iterator].getEndNode();
			}
		}
		return res;
	}
}
