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

package de.jstacs.clustering.distances;

/**
 * This abstract class defined a DistanceMetric (which may be used for clustering) on a generic type <code>T</code>.
 * @author Jan Grau
 *
 * @param <T> the generic type on which the metric is defined
 */
public abstract class DistanceMetric<T> {

	/**
	 * Returns the distance according to the metric of the two supplied objects.
	 * @param o1 the first object
	 * @param o2 the second object
	 * @return the distance
	 * @throws Exception if the distance could not be computed
	 */
	public abstract double getDistance(T o1, T o2) throws Exception;
	
	/**
	 * Returns the matrix of all pairwise distance of the supplied objects, where rows and colums are indexed in the order
	 * of the supplied objects.
	 * @param metric the metric
	 * @param objects the objects
	 * @param <T> the generic type on which the metric is defined
	 * @return the pairwise distances
	 * @throws Exception if the distance could not be computed for one pair of objects
	 */
	public static <T> double[][] getPairwiseDistanceMatrix(DistanceMetric<T> metric, T... objects) throws Exception {
		double[][] matrix = new double[objects.length][];
		for(int i=0;i<matrix.length;i++){
			matrix[i] = new double[i];
			for(int j=0;j<i;j++){
				matrix[i][j] = metric.getDistance( objects[i], objects[j] );
			}
		}
		return matrix;
	}
	
	
}
