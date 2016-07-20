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
 * Class for a topological sort.
 * 
 * @author Jan Grau
 */
public class TopSort {

	/**
	 * Returns the topological order of indexes according to
	 * <code>parents2</code>. The array <code>parents2</code> at each position
	 * <code>i</code> contains the parents, i.e. incoming edges, for this
	 * position at indexes <code>0</code> to <code>parents2[i].length -
	 * 2</code> and <code>i</code> itself at position
	 * <code>parents2[i].length - 1</code>. The returned array of the
	 * topological order is organized as follows:
	 * <ul>
	 * <li> <code>parents[i][0]</code> contains the index in
	 * <code>parents2</code> with number <code>i</code> in the topological order
	 * <li> <code>parents[j][1]</code> contains the number of index
	 * <code>j</code> in the topological order
	 * </ul>
	 * 
	 * @param parents2
	 *            the array of parents
	 * 
	 * @return the topological order
	 */
	public static int[][] getTopologicalOrder( int[][] parents2 ) {
		int[][] parents = new int[parents2.length][];
		for( int i = 0; i < parents.length; i++ ) {
			parents[i] = parents2[i].clone();
		}
		int[][] order = new int[parents.length][2];

		int curr = 0;
		while( curr < order.length ) {
			for( int i = 0; i < parents.length; i++ ) {
				if( parents[i][0] == i ) {

					parents[i][0] = -1;
					order[curr][0] = i;
					order[i][1] = curr;
					curr++;

					for( int j = 0; j < parents.length; j++ ) {
						if( i != j ) {
							checkloop: for( int k = 0; k < parents[j].length - 1; k++ ) {
								if( parents[j][k] == i ) {
									for( int l = k + 1; l < parents[j].length; l++ ) {
										parents[j][l - 1] = parents[j][l];
									}
									break checkloop;
								}
							}
						}
					}
				}
			}
		}
		return order;
	}

	/**
	 * Computes a topological ordering for a given graph.
	 * 
	 * @param aM
	 *            a two-dimensional array describing the graph: <br>
	 *            <ol>
	 *            <li><code>aM.length</code> is equal no number of nodes in the
	 *            graph</li>
	 *            <li><code>aM[i]</code> contains information (parents) of node
	 *            <code>i</code></li>
	 *            <li><code>aM[i]={3,2,i}</code> means, that node <code>i</code>
	 *            has parents 3 and 2</li>
	 *            <li><code>aM[i]={i}</code> means, that node <code>i</code> has
	 *            no parents</li>
	 *            <li>the <code>i</code> at the end of each <code>aM[i]</code>
	 *            is mandatory</li>
	 *            </ol>
	 * 
	 * @return <code>byte</code>-array of length
	 *         &quot;numberOfNodesInGraph&quot; that contains at pos 0 the id of
	 *         the first node in the ordering, at pos 1 the id of the second
	 *         node in the order and so on !!!<br>
	 *         If the graph contained a cyclus or other problems occurred the
	 *         returned array is of length zero.
	 */
	static public byte[] getTopologicalOrder2( byte[][] aM ) {

		//bestimme Vorgaengeranzahlen
		byte[] vorG = new byte[aM.length];
		for( int i = 0; i < vorG.length; vorG[i] = (byte)( aM[i++].length - 1 ) );

		//initiiere Liste aller Knoten, die keine Vorgaenger haben
		byte[] listOfNodesWithZeroPredecessor = new byte[aM.length];
		byte[] numOfNodesWithZeroPredecessor = { 0 };
		for( byte i = 0; i < vorG.length; i++ ) {
			if( vorG[i] == 0 ) {
				listOfNodesWithZeroPredecessor[numOfNodesWithZeroPredecessor[0]] = i;
				numOfNodesWithZeroPredecessor[0]++;
			}
		}

		//wenn es keine Knoten ohne Vorgaenger gibt, gibts auch keine topo. Sortierung
		if( numOfNodesWithZeroPredecessor[0] == 0 ) {
			return new byte[0];
		}

		//finde Element mit Vorgaengeranzahl=0 
		//-> falls es keines gibt aber noch nicht alle behandelt 
		//-> Zyklus in Graph enthalten 
		//-> keine topol. Sortierung moeglich 
		//-> gib Array der Laenge Null zurueck
		byte[] topoSort = new byte[aM.length];
		Arrays.fill( topoSort, (byte)-1 );

		for( int i = 0; i < aM.length; i++ ) {

			topoSort[i] = getNextNodeWithZeroPredecessor( listOfNodesWithZeroPredecessor, numOfNodesWithZeroPredecessor, vorG, aM );

			if( i < ( aM.length - 1 ) && numOfNodesWithZeroPredecessor[0] == 0 ) {
				return new byte[0];
			}
		}

		return topoSort;
	}

	private static byte getNextNodeWithZeroPredecessor( byte[] LONWZP, byte[] NONWZP, byte[] VG, byte[][] aM ) {

		//einen Knoten nehmen, der keine Vorgaenger mehr hat
		byte ret = LONWZP[NONWZP[0] - 1];
		NONWZP[0]--;
		VG[ret] = -1;

		//diesen Knoten aus Graph entfernen und neue Knoten ohne Vorgaenger in Liste aufnehmen
		//Knoten, die nicht mehr im Graphen sind, haben in vorG (VG) eine -1 zu stehen
		for( int i = 0; i < aM.length; i++ ) {

			if( VG[i] == -1 ) {
				continue;
			}

			for( int j = 0; j < aM[i].length - 1; j++ ) {
				if( aM[i][j] == ret ) {
					VG[i]--;
				}
			}

			if( VG[i] == 0 ) {
				VG[i] = -1;
				LONWZP[NONWZP[0]] = (byte)i;
				NONWZP[0]++;
			}
		}

		return ret;
	}

}
