package de.jstacs.utils;

import java.util.Comparator;

/**
 * This class implements a {@link Comparator} of double arrays.
 * It just uses the double value at a specific index of the array to compare two arrays.
 * This index can be set in the constructor.
 * 
 * @author Jens Keilwagen
 */
public class DoubleArrayComparator implements Comparator<double[]> {
	
	private int index;

	/**
	 * Creates a new instance with a specific index for comparing double arrays.
	 * 
	 * @param idx the index of the element which should be used to compare two arrays
	 */
	public DoubleArrayComparator( int idx ) {
		index = idx;
	}
	
	@Override
	public int compare(double[] o1, double[] o2) {
		return (int) Math.signum( o1[index] - o2[index] );
	}		
}