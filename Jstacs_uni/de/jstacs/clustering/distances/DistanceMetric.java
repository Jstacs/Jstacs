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
