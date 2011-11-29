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

/**
 * This class is a representation of a weighted edge.
 * 
 * @author Jens Keilwagen
 */
public class Edge implements Comparable<Edge>, Cloneable {

	/**
	 * The source node.
	 */
	protected int source;

	/**
	 * The target node.
	 */
	protected int target;

	/**
	 * The weight of the edge.
	 */
	protected double weight;

	private Edge( Edge e ) {
		source = e.source;
		target = e.target;
		weight = e.weight;
	}

	/**
	 * Creates a new weighted edge.
	 * 
	 * @param s
	 *            the start node
	 * @param t
	 *            the target (=end) node
	 * @param w
	 *            the weight of the edge
	 * 
	 * @throws IllegalArgumentException
	 *             if the names of the nodes are less than 0
	 */
	public Edge( int s, int t, double w ) throws IllegalArgumentException {
		if( s < 0 || t < 0 ) {
			throw new IllegalArgumentException( "The nodes have to be integer, that are bigger than or equal to zero." );
		}
		source = s;
		target = t;
		weight = w;
	}

	/**
	 * Returns the start node of the edge.
	 * 
	 * @return the start node of the edge
	 */
	public int getStartNode() {
		return source;
	}

	/**
	 * Returns the end node of the edge.
	 * 
	 * @return the end node of the edge
	 */
	public int getEndNode() {
		return target;
	}

	/**
	 * Returns the weight of the edge.
	 * 
	 * @return the weight of the edge
	 */
	public double getWeight() {
		return weight;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "(" + source + ", " + target + "; " + weight + ")";
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo( Edge e ) throws ClassCastException {
		if( weight < e.weight ) {
			return -1;
		} else if( weight == e.weight ) {
			return 0;
		} else {
			return 1;
		}
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	public Edge clone() {
		return new Edge( this );
	}
}
