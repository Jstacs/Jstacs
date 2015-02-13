package de.jstacs.clustering.distances;


public abstract class DistanceMetric<T> {

	public abstract double getDistance(T o1, T o2) throws Exception;
	
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
