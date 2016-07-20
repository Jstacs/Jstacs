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

import de.jstacs.algorithms.graphs.tensor.AsymmetricTensor;
import de.jstacs.algorithms.graphs.tensor.SymmetricTensor;
import de.jstacs.algorithms.graphs.tensor.Tensor;

/**
 * This is the main class of the graph library. Most methods are written for
 * maximization, but they can also be used for minimization if the score is
 * inverted.
 * 
 * <br>
 * <br>
 * 
 * A DAG is a <b>d</b>irected <b>a</b>cyclic <b>g</b>raph. A k-DAG is a DAG with
 * k special nodes that have 0,1,...,k-1 parents and all the other nodes have k
 * parents. A path is a sequence of nodes so that two neighboring nodes are
 * connected by an edge. A HP (Hamiltonian Path) is a path that visits all nodes
 * of an DAG exactly once and maximizes some score function. The HP(k) is the
 * Hamiltonian Path using a score function that depends on the last k visited
 * nodes.
 * 
 * @author Jens Keilwagen
 */
public class DAG {

	/**
	 * The method computes the HP(k) (see {@link DAG}). Be aware of the memory
	 * (O(2^L * L^k)) and time consumption.
	 * 
	 * @param score
	 *            the {@link Tensor} for the edge weights
	 * 
	 * @return the HP(k) (the sequence of nodes)
	 * 
	 * @see Tensor
	 * @see DAG#getStructureFromPath(int[], Tensor)
	 */
	public static int[] computeMaximalHP( Tensor score ) {
		byte k = score.getOrder(), g = (byte)( k - 1 );
		int l, L = score.getNumberOfNodes(), counter, current;
		int numOfSubsets = (int)Math.pow( 2, L ), anz = 1, i, index, oldIndex, oldAnz;
		int[] nodes, copy = new int[L], help = new int[k], perm = new int[k];

		int[] powers = new int[k + 1];
		powers[0] = 1;
		do {
			powers[anz] = powers[anz - 1] * L;
		} while( ++anz <= k );

		// init table
		double best[][];
		int[][] pre;
		if( k == 1 ) {
			best = new double[L][numOfSubsets];
			pre = new int[L][numOfSubsets];
		} else {
			best = new double[powers[k]][];
			pre = new int[powers[k]][];
			nodes = new int[k + 1];
			for( anz = 0; anz < k; anz++ ) {
				nodes[anz] = (byte)( k - anz - 1 );
			}
			do {
				System.arraycopy( nodes, 0, help, 0, k );
				Arrays.sort( help );
				anz = 1;
				i = nodes[0];
				while( anz < k && help[anz] != help[anz - 1] ) {
					i += nodes[anz] * powers[anz++];
				}
				if( anz == k ) {
					best[i] = new double[numOfSubsets];
					pre[i] = new int[numOfSubsets];
					// System.out.println( i + "\t" + toString( nodes, 0, k ) );
				}
				i = 0;
				while( nodes[i] == L - 1 ) {
					nodes[i] = 0;
					i++;
				}
				nodes[i]++;
			} while( nodes[k] == 0 );
		}

		// forward computing
		// the last symbol of a suffix (s_1,s_2,...,s_k) is s_1
		// a suffix s=(s_0,s_1,...,s_{k-1}) is encoded as index[s] =
		// \sum_{l=0}^{k-1} s_l*powers[l]
		// the index of a ancestor with suffix sa=(s_1,s_2,...,s_k) can be
		// computed as
		// => index[sa] = (index[s] - s_0)/L + s_k * powers[k-1] = index[s]/L +
		// s_k * powers[k-1]

		// compute initial values
		nodes = new int[L];
		for( counter = 0; counter < L; counter++ ) {
			nodes[counter] = counter;
		}
		help = new int[k];
		i = 0;
		while( true ) {
			System.arraycopy( nodes, 0, copy, 0, L );
			index = 0;
			anz = 0;
			for( i = 0; i < k; i++ ) {
				counter = help[i];
				perm[i] = copy[counter];
				anz += ( 1 << perm[i] );
				index += perm[i] * powers[g - i];
				while( counter < L - 1 ) {
					copy[counter] = copy[++counter];
				}
			}
			// System.out.println( anz + " " + index + "\t" + toString( perm, 0, k ) );
			best[index][anz] = getScoreForPath( score, k, k, perm );
			pre[index][anz] = -1;

			i = g;
			while( i >= 0 && help[i] == L - i - 1 ) {
				help[i--] = 0;
			}
			if( i >= 0 ) {
				help[i]++;
			} else {
				break;
			}
		}

		// compute *higher* values
		double val;
		for( l = (byte)( k + 1 ); l <= L; l++ ) {
			// sublength l => enumerate all subsets
			// System.out.println( "Teil-Laenge: " + l + " von " + L );
			anz = 0;
			for( counter = 0; counter < l; counter++ ) {
				nodes[counter] = counter;
				anz += ( 1 << counter );
			}
			while( true ) {
				// System.out.println( anz + " " + toString(nodes, 0, l) );
				// test all suffixes
				help = new int[k];
				i = 0;
				while( i >= 0 ) {
					System.arraycopy( nodes, 0, copy, 0, l );
					// the current node (last node of the partial PMM-structure)
					counter = help[0];
					current = copy[counter];
					while( counter < l - 1 ) {
						copy[counter] = copy[++counter];
					}
					index = current;
					for( i = 1; i < k; i++ ) {
						counter = help[i];
						perm[i - 1] = copy[counter];
						index += copy[counter] * powers[i];
						while( counter < l - 1 ) {
							copy[counter] = copy[++counter];
						}
					}
					// System.out.println( anz + " " + index + " " + current + " " + toString(perm, 0, k-1) );

					best[index][anz] = Double.NEGATIVE_INFINITY;
					oldIndex = index / L;
					oldAnz = anz - ( 1 << current );

					for( counter = 0; counter < l - k; counter++ ) {
						perm[g] = copy[counter];
						// System.out.println( current + " <- " + Arrays.toString( perm ) );
						val = best[oldIndex + perm[g] * powers[g]][oldAnz] + score.getValue( k, current, perm );
						if( best[index][anz] < val ) {
							best[index][anz] = val;
							pre[index][anz] = perm[g];
						}
					}

					i = g;
					while( i >= 0 && help[i] == l - i - 1 ) {
						help[i--] = 0;
					}
					if( i >= 0 ) {
						help[i]++;
					}
				}

				i = l - 1;
				counter = (byte)( L - 1 );
				while( i >= 0 && nodes[i] == counter ) {
					anz -= ( 1 << nodes[i--] );
					counter--;
				}
				if( i >= 0 ) {
					anz = anz - ( 1 << nodes[i] ) + ( 1 << ++nodes[i++] );
					while( i < l ) {
						nodes[i] = nodes[i - 1];
						anz += ( 1 << ++nodes[i++] );
					}
				} else {
					break;
				}
			}
		}

		// find max
		double max = Double.NEGATIVE_INFINITY;
		anz = numOfSubsets - 1;
		for( i = 0; i < powers[k]; i++ ) {
			if( best[i] != null && best[i][anz] > max ) {
				index = i;
				max = best[i][anz];
			}
		}

		int[] erg = new int[L];
		i = L - 1;
		while( pre[index][anz] != -1 ) {
			erg[i] = ( index % L );
			index = index / L + pre[index][anz] * powers[g];
			anz -= ( 1 << erg[i--] );
		}
		anz /= L;
		while( i >= 0 ) {
			erg[i--] = ( index % L );
			index = index / L;
		}

		return erg;
	}

	/**
	 * Computes the maximal k-DAG (see {@link DAG}), i.e. the k-DAG that
	 * maximizes the score given by a {@link Tensor}.
	 * 
	 * <br>
	 * <br>
	 * 
	 * This implementation uses a DP (<b>d</b>ynamic <b>p</b>rogramming)
	 * algorithm.
	 * 
	 * @param score
	 *            the {@link Tensor} of the score values
	 * 
	 * @return the maximal k-DAG
	 * 
	 * @see Tensor
	 */
	public static int[][] computeMaximalKDAG( Tensor score ) {
		return computeMaxKDAG( ( score instanceof SymmetricTensor )	? ( (SymmetricTensor)score )
																	: ( new SymmetricTensor( (AsymmetricTensor)score ) ) );
	}

	/**
	 * Computes the maximal k-DAG (see {@link DAG}), i.e. the k-DAG that
	 * maximizes the score given by a {@link SymmetricTensor}.
	 * 
	 * <br>
	 * <br>
	 * 
	 * This implementation uses a DP (<b>d</b>ynamic <b>p</b>rogramming)
	 * algorithm.
	 * 
	 * @param score
	 *            the {@link SymmetricTensor} of the score values
	 * 
	 * @return the maximal k-DAG
	 * 
	 * @see SymmetricTensor
	 */
	protected static int[][] computeMaxKDAG( SymmetricTensor score ) {
		// variables
		int L = score.getNumberOfNodes(), n = 0, i, l = 0;
		byte k = score.getOrder(), c;

		if( L > 31 ) {
			throw new IllegalArgumentException( "This algorithm can only handle bayesian networks with at most 31 nodes. If you are interested in bayesian networks with more nodes please try an other algorithm." );
		}

		int anz = 1, numOfSubsets = (int)Math.pow( 2, L );
		int[] nodes = new int[L + 1], par;
		nodes[0] = 0;
		double[] best = new double[numOfSubsets];
		int[] edge = new int[numOfSubsets];
		double temp;

		// forward computing of the score
		do {
			// System.out.println( anz + " " + StringParser.toString(
			// nodes,0,l+1 ) + " " );
			if( l == 0 ) {
				best[anz] = score.getRootValue( nodes[0] );
				edge[anz] = 0;
			} else {
				c = ( k < l ) ? k : (byte)l;
				par = new int[l];
				System.arraycopy( nodes, 0, par, 0, l );
				best[anz] = best[anz - ( 1 << nodes[--n] )] + score.getBest( nodes[n], par, c );
				edge[anz] = score.getTrueIndexForLastGetBest();
				while( n > 0 ) {
					par[--n] = ( nodes[n + 1] - 1 );
					temp = best[anz - ( 1 << nodes[n] )] + score.getBest( nodes[n], par, c );
					if( best[anz] < temp ) {
						best[anz] = temp;
						edge[anz] = score.getTrueIndexForLastGetBest();
					}
				}
			}
			// System.out.println( best[anz] );

			i = 0;
			while( ( anz & ( 1 << i ) ) > 0 ) {
				i++;
			}
			anz++;
			l -= i - 1;
			nodes[0] = i;
			for( n = 1; n <= l; n++ ) {
				while( ( anz & ( 1 << ++i ) ) == 0 );
				nodes[n] = i;
			}
		} while( anz < numOfSubsets );
		// if( L <=10 ) System.out.println( Arrays.toString(best) );
		// System.out.println( "score: " + best[best.length - 1] );

		// backward computing of the k-dag-structure (result)
		int[][] tempRandDeps = new int[L][];
		int[] tempEdge;
		anz--;
		l = ( L - 1 );
		do {
			i = ( k < l ) ? k : l;
			tempEdge = score.getEdgeFromIndex( edge[anz], (byte)( i + 1 ) );
			tempRandDeps[tempEdge[i]] = tempEdge;
			anz -= ( 1 << tempEdge[i] );
			l--;
		} while( l != 0 );
		tempRandDeps[tempEdge[0]] = new int[]{ tempEdge[0] };

		return tempRandDeps;
	}

	/**
	 * The method computes the HP(k) (see {@link DAG}). It tries to enumerate
	 * all permutations and returns the best. This is only possible for a small
	 * number of nodes (&lt;=10).
	 * 
	 * @param score
	 *            the {@link Tensor} for the edge weights
	 * 
	 * @return the HP(k) (the sequence of nodes)
	 * 
	 * @see DAG#getStructureFromPath(int[], Tensor)
	 */
	public static int[] enumerateHP( Tensor score ) {
		int L = score.getNumberOfNodes(), m = L, counter = 0, i;
		byte k = score.getOrder();
		m--;
		int[] copy = new int[L], index = new int[L], bestPerm = new int[L], perm = new int[L], id = new int[L];
		index[0] = -1;
		for( ; counter < L; counter++ ) {
			id[counter] = counter;
		}
		double weights, bestWeights = Double.NEGATIVE_INFINITY;

		counter = 0;
		do {
			index[counter]++;
			// compute perm
			System.arraycopy( id, 0, copy, 0, L );
			for( counter = 0; counter < L; counter++ ) {
				i = index[counter];
				perm[counter] = copy[i];
				while( i < m ) {
					copy[i++] = copy[i];
				}
			}
			// compute score for perm
			weights = getScoreForPath( score, L, k, perm );

			// System.out.println( Arrays.toString(perm) + " " + weights );

			if( weights > bestWeights ) {
				bestWeights = weights;
				System.arraycopy( perm, 0, bestPerm, 0, L );
			}

			// compute next index
			counter = m;
			while( --counter >= 0 && index[counter] == m - counter ) {
				index[counter] = 0;
			}
		} while( counter >= 0 );
		return bestPerm;
	}

	/**
	 * Returns the score for any graph.
	 * 
	 * @param t
	 *            the {@link Tensor} for the edge weights
	 * @param structure
	 *            the graph (encoded as:
	 * 
	 *            <code>(structure[i][0],...,structure[i][structure[i].length-2])-&gt;structure[i][structure[i].length-1]</code>
	 *            )
	 * 
	 * @return the score for the graph
	 * 
	 * @throws IllegalArgumentException
	 *             if the structure length and the tensor length do not match
	 */
	public static double getScore( Tensor t, int[][] structure ) throws IllegalArgumentException {
		int counter = 0, k, n = t.getNumberOfNodes();
		if( structure.length != n ) {
			throw new IllegalArgumentException( "The given structure and tensor are defined on a different number of nodes." );
		}
		double erg = 0;
		for( ; counter < n; counter++ ) {
			k = structure[counter].length - 1;
			if( k > 0 ) {
				erg += t.getValue( (byte)k, structure[counter][k], structure[counter] );
			} else if( k == 0 ) {
				erg += t.getRootValue( structure[counter][k] );
			}
		}
		return erg;
	}

	/**
	 * Returns the score for a given path <code>path</code> using the first
	 * <code>l</code> nodes and dependencies of order <code>k</code>.
	 * 
	 * @param score
	 *            the {@link Tensor} of scores
	 * @param l
	 *            the number of used nodes (from <code>perm</code>)
	 * @param k
	 *            the order
	 * @param path
	 *            the path
	 * 
	 * @return the score for the given path
	 */
	public static double getScoreForPath( Tensor score, int l, byte k, int[] path ) {
		double weight = score.getRootValue( path[0] );
		int[] parents = new int[k];
		byte i = 1;
		int j, c;
		for( ; i < k; i++ ) {
			for( c = 0, j = i - 1; j >= 0; j--, c++ ) {
				parents[c] = path[j];
			}
			weight += score.getValue( i, path[i], parents );
		}
		for( ; i < l; i++ ) {
			for( c = 0, j = i - 1; j >= i - k; j--, c++ ) {
				parents[c] = path[j];
			}
			weight += score.getValue( k, path[i], parents );
		}
		return weight;
	}

	/**
	 * Extracts the structure from a given path <code>path</code> and
	 * score-&quot;function&quot;.
	 * 
	 * @param path
	 *            the path
	 * @param score
	 *            the {@link Tensor} of scores
	 * 
	 * @return the structure of the given path (the edges)
	 */
	public static int[][] getStructureFromPath( int[] path, Tensor score ) {
		byte k = score.getOrder();
		int i = 1, l = score.getNumberOfNodes();
		int[][] structure = new int[l][];
		int[] parents = new int[k];
		int j, c;
		/*
		 * structure[0] = new byte[] { path[0] }; for( ; i < k; i++ ) { for( c = 0, j = i - 1; j >= 0; j--, c++ ) {
		 * parents[c] = path[j]; } structure[i] = score.getMaximalEdgeFor( path[i], parents, i ); } for( ; i < l; i++ ) {
		 * for( c = 0, j = i - 1; j >= i - k; j--, c++ ) { parents[c] = path[j]; } structure[i] =
		 * score.getMaximalEdgeFor( path[i], parents, k ); }
		 */
		structure[path[0]] = new int[]{ path[0] };
		for( ; i < k; i++ ) {
			for( c = 0, j = i - 1; j >= 0; j--, c++ ) {
				parents[c] = path[j];
			}
			structure[path[i]] = score.getMaximalEdgeFor( (byte)i, path[i], parents );
		}
		for( ; i < l; i++ ) {
			for( c = 0, j = i - 1; j >= i - k; j--, c++ ) {
				parents[c] = path[j];
			}
			structure[path[i]] = score.getMaximalEdgeFor( k, path[i], parents );
		}

		return structure;
	}

	/**
	 * This method returns a directed {@link String} representation of the
	 * structure that can be used in <i>Graphviz</i> to create an image.
	 * 
	 * @param structure
	 *            the structure of the graph
	 * 
	 * @return a directed {@link String} representation
	 * 
	 * @see DAG#toGraphvizFormat(int[][], String)
	 */
	public static String toDirectedGraphvizFormat( int[][] structure ) {
		return toGraphvizFormat( structure, " -> " );
	}

	/**
	 * This method returns an undirected {@link String} representation of the
	 * structure that can be used in <i>Graphviz</i> to create an image.
	 * 
	 * @param structure
	 *            the structure of the graph
	 * 
	 * @return an undirected {@link String} representation
	 * 
	 * @see DAG#toGraphvizFormat(int[][], String)
	 */
	public static String toUndirectedGraphvizFormat( int[][] structure ) {
		return toGraphvizFormat( structure, " -- " );
	}

	/**
	 * This method returns a {@link String} representation of the weighted
	 * structure that can be used in <i>Graphviz</i> to create an image.
	 * 
	 * @param structure
	 *            the structure of the graph
	 * @param arrow
	 *            the kind of arrow that is used between the nodes, e.g.
	 *            &quot;--&quot; and &quot;-&gt;&quot;
	 * @param t
	 *            the weights
	 * 
	 * @return a {@link String} representation of the weighted structure
	 */
	public static String toWeightedGraphvizFormat( int[][] structure, String arrow, Tensor t ) {
		String erg = "";
		int j, k;
		double current, sum = 0;
		for( byte i = 0; i < structure.length; i++ ) {
			// System.out.println( i );System.out.println( Arrays.toString(
			// structure[i] ) );
			j = 0;
			current = 0;
			if( structure[i].length > 2 ) {
				erg += "{" + structure[i][0];
				for( j = 1; j < structure[i].length - 1; j++ ) {
					erg += " " + structure[i][j];
				}
				erg += "} " + arrow + " ";
			} else if( structure[i].length == 2 ) {
				erg += structure[i][j++] + arrow;
			}
			erg += structure[i][j];

			k = structure[i].length - 1;
			if( k > 0 ) {
				current = t.getValue( (byte)k, structure[i][k], structure[i] );
			} else if( k == 0 ) {
				current = t.getRootValue( structure[i][k] );
			}
			erg += "\t[weight=" + current + "]\n";
			sum += current;
		}
		erg += "\nscore = " + sum;
		return erg;
	}

	/**
	 * This method returns a {@link String} representation of the structure that
	 * can be used in <i>Graphviz</i> to create an image.
	 * 
	 * @param structure
	 *            the structure of the graph
	 * @param arrow
	 *            the kind of arrow that is used between the nodes, e.g.
	 *            &quot;--&quot; and &quot;-&gt;&quot;
	 * 
	 * @return a {@link String} representation of the structure
	 */
	protected static String toGraphvizFormat( int[][] structure, String arrow ) {
		String erg = "";
		int j;
		for( int i = 0; i < structure.length; i++ ) {
			// System.out.println( i );System.out.println( Arrays.toString(
			// structure[i] ) );
			j = 0;
			if( structure[i].length > 2 ) {
				erg += "{" + structure[i][0];
				for( j = 1; j < structure[i].length - 1; j++ ) {
					erg += " " + structure[i][j];
				}
				erg += "} " + arrow + " ";
			} else if( structure[i].length == 2 ) {
				erg += structure[i][j++] + arrow;
			}
			erg += structure[i][j] + "\n";
		}
		return erg;
	}
}
