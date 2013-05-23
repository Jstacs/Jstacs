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
 * This class implements the algorithm of Chu_Liu_Edmonds to determine an
 * optimal branching (optimal directed graph of order 1).
 * 
 * @author Andre Gohr
 */
public class Chu_Liu_Edmonds {

	/**
	 * Compute the branching yielding the maximum sum of weights.
	 */
	final public static byte MAXIMALBRANCHING = 0;

	/**
	 * Compute the branching yielding the minimum sum of weights.
	 */
	final public static byte MINIMALBRANCHING = 1;

	private final static double[] INFTY = { Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY };

	private final static byte[] NULLBYTEARRAY = new byte[0];

	private Chu_Liu_Edmonds() {
		super();
	}

	/**
	 * Returns an optimal branching.
	 * 
	 * @param graph
	 *            <ul>
	 *            <li>of dimension
	 *            <code>[number_of_nodes][number_of_nodes]</code>;
	 *            <li>e.g. <code>[1][2]</code> contains the weight of edge
	 *            (1,2), if node 2 is parent of node 1
	 *            <li>the weight of non-existing edges should be
	 *            <code>Double.POSITIVE_INFINITY</code>/
	 *            <code>Double.NEGATIVE_INFINITY</code> if a minimal/maximal
	 *            branching is desired
	 *            </ul>
	 * @param rootWeights
	 *            of dimension <code>[number_of_nodes][1]</code> and contains
	 *            for each node the additional costs, if this node becomes a
	 *            root of the optimal branching. Whether the returned optimal
	 *            branching is a optimal directed spanning tree has to be
	 *            checked by the client.
	 * @param type
	 *            either <code>Chu_Liu_Edmonds.MAXIMALBRANCHING</code> or
	 *            <code>Chu_Liu_Edmonds.MINIMALBRANCHING</code>
	 * 
	 * @return an optimal branching (probable a forest) decoded using a
	 *         <code>byte[][]</code>. The first dimension is defined by the
	 *         number of nodes N. If node i has parent <code>j ->
	 *         byte[i]</code>={j,i}; if node i has no parent -> <code>byte[i]</code>={i}. The last position
	 * 		   always contains the child node. The positions before contains the parent node if any.
	 * 
	 * @throws IllegalArgumentException
	 *             if the type was chosen wrongly
	 * @throws Exception
	 *             if something went wrong
	 */
	static public byte[][] getOptimalBranching( double[][] graph, double[][] rootWeights, byte type ) throws IllegalArgumentException,
			Exception {
		if( !( type == MAXIMALBRANCHING ^ type == MINIMALBRANCHING ) ) {
			throw new IllegalArgumentException( "Error in Chu_Liu_Edmonds.getOptimalBranching(): " + "given type has to be either Chu_Liu_Edmonds.MAXIMALBRANCHING or Chu_Liu_Edmonds.MINIMALBRANCHING" );
		}

		byte root = -1;
		double[][] optimalResult = null;
		// die Wurzel spielt keine Rolle
		if( rootWeights == null || rootWeights.length == 0 ) {
			// graph wird waehrend Abarbeitung veraendert
			optimalResult = Chu_Liu_Edmonds_Algo( copyGraph( graph, type ), root, type );
		} else {
			double[][] result;
			double tempScore, optimalScore = INFTY[type];
			for( byte i = 0; i < rootWeights.length; i++ ) {
				root = i;
				result = Chu_Liu_Edmonds_Algo( copyGraph( graph, type ), root, type );

				tempScore = 0;
				for( byte l = 0; l < result.length; l++ ) {
					for( byte j = 0; j < result[l].length; j++ ) {
						// posteriorM[i][j] war immer i<-j (j ist parent von i
						// -> kann also auch als j->i aufgefasst werden)
						if( !Double.isInfinite( result[l][j] ) ) {
							tempScore += result[l][j];
							break;
						}
					}
				}
				if( ( type == MAXIMALBRANCHING && tempScore + rootWeights[root][0] > optimalScore ) || ( type == MINIMALBRANCHING && tempScore + rootWeights[root][0] < optimalScore ) ) {
					optimalScore = tempScore + rootWeights[root][0];
					optimalResult = result;
				}
			}
		}

		boolean foundDep;
		byte[][] randVarDeps = new byte[graph.length][];

		for( byte l = 0; l < optimalResult.length; l++ ) {
			foundDep = false;
			for( byte j = 0; j < optimalResult[l].length; j++ ) {
				// posteriorM[i][j] war immer i<-j (j ist parent von i -> kann
				// also auch als j->i aufgefasst werden)
				if( !Double.isInfinite( optimalResult[l][j] ) ) {
					byte[] temp = new byte[2];
					temp[0] = j;
					temp[1] = l;
					randVarDeps[l] = temp;
					foundDep = true;
					break;
				}
			}
			if( !foundDep ) {
				byte[] temp = new byte[1];
				temp[0] = l;
				randVarDeps[l] = temp;
			}
		}
		return randVarDeps;
	}

	private static double[][] Chu_Liu_Edmonds_Algo( double[][] graph, byte root, byte type ) throws Exception {

		/*
		 * System.out.println("an Chu_Liu_Edmonds uebergebene Matrix:");
		 * printMatrix(graph);
		 */

		byte[] parents = new byte[graph.length];
		Arrays.fill( parents, (byte)-1 );
		byte tempParent;
		byte[] tempCycle;

		// finde optimale einlaufende Kanten fuer jeden Knoten
		for( byte i = 0; i < graph.length; i++ ) {
			if( i != root ) {
				tempParent = getMaxParent( i, graph, type );
				parents[i] = tempParent;
				// System.out.println("if-case: parents["+i+"]="+tempParent);
			} else {
				parents[i] = -1;
				// System.out.println("else-case: parents["+i+"]=-1");
			}

		}

		// die entstandenen Zyklen finden
		byte[][] cycli = getCycles( parents );
		if( cycli.length == 0 ) {
			tempCycle = NULLBYTEARRAY;
		} else {
			tempCycle = cycli[0];
		}

		// System.out.println("found Cycle: "+cycli.length);
		// for(int i=0; i<cycli.length;i++){
		// System.out.println("cyclus "+i);
		// printCycle(cycli[i]);
		// }
		// System.exit(0);

		// falls ein Zyklus gefunden wurde
		if( tempCycle.length != 0 ) {

			// System.out.println("Zyklus gefunden");

			// falls gesamter Graph aus einem grossen Zyklus besteht oder falls
			// eine Wurzel gewaehlt wurde und der Graph nur noch
			// aus 2 Knoten besteht (Wurzel und noch ein anderer) -> dann kann
			// Graph keinen Zyklus enthalten
			if( tempCycle.length == graph.length ) {

				// System.out.println("Zyklus umfasst gesamten Graphen");

				// nur unoptimalste Kante aus Zyklus loeschen
				byte unoptimalRandVar = (byte)getMostUnoptimalParams( graph, tempCycle, parents, type )[0];
				parents[unoptimalRandVar] = -1;

				// aus dem gegebenen Graphen die minimalste Zykluskante und alle
				// anderen nicht am dMST beteiligten Kanten
				// loeschen
				reduceGraph( graph, parents, type );

				return graph;
			}// gesamter Graph ist grosser Zyklus

			if( tempCycle.length == 1 ) {
				throw new Exception( "Error in Chu_Liu_Edmonds: Cycle has length 1" );
			}

			/** ****************** */
			/** Zykluskontraktion* */
			/** ****************** */

			// System.out.println("Zyklus umfasst nur einen Teil des gesamten
			// Graphen");
			Arrays.sort( tempCycle );

			/*
			 * System.out.println("Zyklus:"); printCycle(tempCycle);
			 */

			// unoptimalste Zykluskante finden
			double[] temp = getMostUnoptimalParams( graph, tempCycle, parents, type );
			double cycleMostUnoptimalScore = temp[2];
			byte mostUnoptimalRandVar = (byte)temp[0], mostUnoptimalParent = (byte)temp[1];

			// Matrix umschreiben, so dass ZyklusKnoten ganz hinten und ganz
			// unten
			// dies geschieht nur virtuell
			byte tempCycleCounter = 0, newRandVarOrderCounter = 0, newPosOfRoot = -1;
			byte[] newRandVarOrder = new byte[graph.length];
			for( byte m = 0; m < newRandVarOrder.length; m++ ) {
				if( tempCycleCounter >= tempCycle.length || m != tempCycle[tempCycleCounter] ) {
					newRandVarOrder[newRandVarOrderCounter++] = m;
					if( m == root ) {
						newPosOfRoot = (byte)( newRandVarOrderCounter - 1 );
					}
				} else {
					newRandVarOrder[newRandVarOrder.length - 1 - tempCycleCounter] = m;
					tempCycleCounter++;
				}
			}

			double[][] reducedMatrix = new double[graph.length - tempCycle.length + 1][graph.length - tempCycle.length + 1];
			byte[] maxHeadFromIngoingEdge = new byte[reducedMatrix.length - 1];
			double tempScore;
			// die maximalen ausgehenden und einlaufenden Kanten in Zyklus und
			// aus Zyklus bestimmen
			// und alle anderen Knaten in new Matrix auf NaN setzten
			// die maximalen ausgehenden Kanten finden und setzten
			int cycleRegion = newRandVarOrder.length - tempCycle.length;
			for( byte m = 0; m < newRandVarOrder.length; m++ ) {
				for( byte l = 0; l < newRandVarOrder.length; l++ ) {
					// der Kern
					if( m < cycleRegion && l < cycleRegion ) {
						// die Diagonale
						if( m == l ) {
							reducedMatrix[m][l] = INFTY[type];
						} else {
							// der Rest
							reducedMatrix[m][l] = graph[newRandVarOrder[m]][newRandVarOrder[l]];
						}
					} else {
						// die maximal/minimal aus Zyklus auslaufenden Kanten
						if( m < cycleRegion && l >= cycleRegion ) {
							if( l == cycleRegion ) {
								reducedMatrix[m][cycleRegion] = graph[newRandVarOrder[m]][newRandVarOrder[l]];
							} else {
								switch( type ) {
									case ( MAXIMALBRANCHING ): {
										if( reducedMatrix[m][cycleRegion] < graph[newRandVarOrder[m]][newRandVarOrder[l]] ) {
											reducedMatrix[m][cycleRegion] = graph[newRandVarOrder[m]][newRandVarOrder[l]];
										}
										break;
									}
									case ( MINIMALBRANCHING ): {
										if( reducedMatrix[m][cycleRegion] > graph[newRandVarOrder[m]][newRandVarOrder[l]] ) {
											reducedMatrix[m][cycleRegion] = graph[newRandVarOrder[m]][newRandVarOrder[l]];
										}
										break;
									}
								}
							}
						} else {
							// die maximal/minimal in Zyklus einlaufenden Kanten
							// mit neuen Gewichten
							if( m >= cycleRegion && l < cycleRegion ) {
								if( m == cycleRegion ) {
									tempScore = graph[newRandVarOrder[m]][newRandVarOrder[l]] + cycleMostUnoptimalScore
												- graph[newRandVarOrder[m]][parents[newRandVarOrder[m]]];
									reducedMatrix[cycleRegion][l] = tempScore;
									maxHeadFromIngoingEdge[l] = m;
								} else {
									tempScore = graph[newRandVarOrder[m]][newRandVarOrder[l]] + cycleMostUnoptimalScore
												- graph[newRandVarOrder[m]][parents[newRandVarOrder[m]]];

									switch( type ) {
										case ( MAXIMALBRANCHING ): {
											if( reducedMatrix[cycleRegion][l] < tempScore ) {
												reducedMatrix[cycleRegion][l] = tempScore;
												maxHeadFromIngoingEdge[l] = m;
											}
											break;
										}
										case ( MINIMALBRANCHING ): {
											if( reducedMatrix[cycleRegion][l] > tempScore ) {
												reducedMatrix[cycleRegion][l] = tempScore;
												maxHeadFromIngoingEdge[l] = m;
											}
											break;
										}
									}
								}
							} else {
								// die nicht optimalen Kanten innerhalb des
								// Zyklus loeschen
								if( newRandVarOrder[l] != parents[newRandVarOrder[m]] ) {
									graph[newRandVarOrder[m]][newRandVarOrder[l]] = INFTY[type];
								}
							}
						}
					}
				}
			}

			reducedMatrix[reducedMatrix.length - 1][reducedMatrix.length - 1] = INFTY[type];

			/*
			 * System.out.println("reduzieter Matrix:");
			 * printMatrix(reducedMatrix);
			 */

			// System.exit(0);
			// System.out.println("\n ***Aufruf Chu_Liu_Edmonds***\n");
			double[][] newWeights = Chu_Liu_Edmonds_Algo( reducedMatrix, newPosOfRoot, type );

			/*
			 * System.out.println("zurueckgegebene Matrix (enthaelt partiellen
			 * optimales branching"); printMatrix(newWeights);
			 */

			/** ******************** */
			/** Zyklusdekontraktion* */
			/** ******************** */

			// Knotenkontraktion rueckgaengig machen
			boolean exsistsIngoingEdge = false;
			byte headNodeOfIngoingEdge = -1;
			for( byte m = 0; m < newRandVarOrder.length; m++ ) {
				for( byte l = 0; l < newRandVarOrder.length; l++ ) {
					// Kern der zurueckerhaltenen Matrix kopieren
					if( m < cycleRegion && l < cycleRegion ) {
						graph[newRandVarOrder[m]][newRandVarOrder[l]] = reducedMatrix[m][l];
					} else {
						// die aus Zyklus herauslaufenden Kanten setzten
						if( m < cycleRegion && l >= cycleRegion ) {
							if( Double.isInfinite( reducedMatrix[m][cycleRegion] ) || graph[newRandVarOrder[m]][newRandVarOrder[l]] != reducedMatrix[m][cycleRegion] ) {
								graph[newRandVarOrder[m]][newRandVarOrder[l]] = INFTY[type];
							} else {
								graph[newRandVarOrder[m]][newRandVarOrder[l]] = reducedMatrix[m][cycleRegion];
							}
						} else {
							// die in Zyklus hineinlaufenden Kanten setzten
							if( m >= cycleRegion && l < cycleRegion ) {
								if( Double.isInfinite( reducedMatrix[cycleRegion][l] ) ) {
									graph[newRandVarOrder[m]][newRandVarOrder[l]] = INFTY[type];
								} else {
									exsistsIngoingEdge = true;
									if( m != maxHeadFromIngoingEdge[l] ) {
										graph[newRandVarOrder[m]][newRandVarOrder[l]] = INFTY[type];
									} else {
										headNodeOfIngoingEdge = newRandVarOrder[m];
									}
								}
							}
						}
					}
				}
			}

			// fuer die neue in Zyklus einlaufende Kante muss eine
			// konkurrierdene zyklusinterne Kante geloescht werden
			if( exsistsIngoingEdge ) {
				graph[headNodeOfIngoingEdge][parents[headNodeOfIngoingEdge]] = INFTY[type];
			} else {
				// falls es keine in den Zyklus einlaufende Kante von ausserhalb
				// des Zyklus gibt
				// muss noch die minimalste Kante des Zyklus geloescht werden
				graph[mostUnoptimalRandVar][mostUnoptimalParent] = INFTY[type];
			}

			return graph;
		}// falls ein Zyklus gefunden wurde

		// falls kein Zyklus gefunden wurde (entweder wurde ein optimaler dMST
		// gefunden oder ein
		// optimales branching (moeglich, falls per Graphdefinition der Graph
		// nicht zusammenhaengend ist)
		// deswegen wurde bei Auswahl der je maximalen Kanten kein
		// zusammenhaengendes Etwas gefunden
		// aus dem gegebenen Graphen die minimalste Zykluskante und alle anderen
		// nicht am dMST oder am optimalen branching beteiligten Kanten
		// loeschen
		reduceGraph( graph, parents, type );

		return graph;
	}

	private static byte[][] getCycles( byte[] parents ) {

		boolean[] traversed = new boolean[parents.length];
		boolean[] tempWalk = new boolean[parents.length];
		Arrays.fill( traversed, false );
		int cycleLength;
		byte tempRandVar;
		java.util.LinkedList cycli = new java.util.LinkedList();

		for( byte randVar = 0; randVar < parents.length; randVar++ ) {
			if( !traversed[randVar] ) {
				// System.out.println("test Var:"+randVar);
				tempRandVar = randVar;
				Arrays.fill( tempWalk, false );

				while( tempRandVar != -1 && !tempWalk[tempRandVar] && !traversed[tempRandVar] ) {
					tempWalk[tempRandVar] = true;
					// System.out.print(" ->"+tempRandVar);
					tempRandVar = parents[tempRandVar];
				}

				// keinen Zyklus gefunden
				if( tempRandVar == -1 || traversed[tempRandVar] ) {
					// System.out.println("keinen Z gefunden");
				}
				// einen Zyklus gefunden
				else {
					// System.out.println("Z gefunden");
					cycleLength = 0;
					byte cycleElement = tempRandVar;
					do {
						tempRandVar = parents[tempRandVar];
						cycleLength++;
					} while( tempRandVar != cycleElement );

					byte[] cyclus = new byte[cycleLength];
					tempRandVar = cycleElement;
					for( byte i = 0; i < cyclus.length; i++ ) {
						cyclus[i] = tempRandVar;
						tempRandVar = parents[tempRandVar];
					}
					cycli.add( cyclus );
				}

				for( int i = 0; i < traversed.length; traversed[i] = ( traversed[i] | tempWalk[i++] ) );

			}
		}

		if( cycli.size() == 0 ) {
			return new byte[0][0];
		} else {
			int size = cycli.size();
			byte[][] ret = new byte[size][];
			for( int i = 0; i < size; i++ ) {
				ret[i] = (byte[])cycli.removeFirst();
			}
			return ret;
		}

	}

	private static byte getMaxParent( byte randVar, double[][] weights, byte type ) {

		byte optimalParent = -1;
		double optimalScore = INFTY[type];

		switch( type ) {
			case ( MINIMALBRANCHING ): {

				for( byte i = 0; i < weights[randVar].length; i++ ) {
					if( i != randVar && weights[randVar][i] < optimalScore ) {
						optimalScore = weights[randVar][i];
						optimalParent = i;
					}
				}
				break;
			}

			case ( MAXIMALBRANCHING ): {

				for( byte i = 0; i < weights[randVar].length; i++ ) {
					if( i != randVar && weights[randVar][i] > optimalScore ) {
						optimalScore = weights[randVar][i];
						optimalParent = i;
					}
				}
				break;
			}
		}
		return optimalParent;
	}

	private static void reduceGraph( double[][] graph, byte[] parents, byte type ) {

		// aus dem gegebenen Graphen die minimalste Zykluskante und alle anderen
		// nicht am dMST beteiligten Kanten
		// loeschen
		for( byte m = 0; m < graph.length; m++ ) {
			for( byte l = 0; l < graph.length; l++ ) {
				if( m != l && l != parents[m] ) {
					graph[m][l] = INFTY[type];
				}
			}
		}
	}

	// {mostUnoptimalRandVar,mostUnoptimalParent,mostUnoptimalEdgeWeight}
	private static double[] getMostUnoptimalParams( double[][] graph, byte[] tempCycle, byte[] parents, byte type ) {

		double mostUnoptimalScore = INFTY[1 - type];
		double[] ret = new double[3];
		// unoptimalste Kante finden
		for( byte m = 0; m < tempCycle.length; m++ ) {
			switch( type ) {
				case ( MAXIMALBRANCHING ): {
					if( graph[tempCycle[m]][parents[tempCycle[m]]] < mostUnoptimalScore ) {
						mostUnoptimalScore = graph[tempCycle[m]][parents[tempCycle[m]]];
						ret[0] = tempCycle[m];
						ret[1] = parents[tempCycle[m]];
					}
					break;
				}
				case ( MINIMALBRANCHING ): {
					if( graph[tempCycle[m]][parents[tempCycle[m]]] > mostUnoptimalScore ) {
						mostUnoptimalScore = graph[tempCycle[m]][parents[tempCycle[m]]];
						ret[0] = tempCycle[m];
						ret[1] = parents[tempCycle[m]];
					}
					break;
				}
			}
		}

		ret[2] = mostUnoptimalScore;
		return ret;
	}

	private static double[][] copyGraph( double[][] graph, byte type ) {

		double[][] ret = new double[graph.length][graph.length];

		for( int i = 0; i < graph.length; i++ ) {
			for( int j = 0; j < graph.length; j++ ) {
				if( i == j ) {
					ret[i][j] = INFTY[type];
				} else {
					ret[i][j] = graph[i][j];
				}
			}
		}

		return ret;
	}

	private static void printMatrix( double[][] m ) {

		System.out.println( "matrixausgabe:***********************" );
		for( byte i = 0; i < m.length; i++ ) {
			for( byte j = 0; j < m[i].length; j++ ) {
				System.out.print( m[i][j] + "\t" );
			}
			System.out.println();
		}
		System.out.println( "ende*********************************" );
	}

	private static void printCycle( byte[] c ) {

		System.out.println( "start Zyklus:****************" );
		for( int i = 0; i < c.length; i++ ) {
			System.out.print( c[i] + "\t" );
		}
		System.out.println( "\nende Zyklus******************" );
	}

}
