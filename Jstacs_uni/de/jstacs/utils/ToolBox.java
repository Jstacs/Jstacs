/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.utils;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Random;

/**
 * This class is a collection of methods which might be useful for the programmer.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class ToolBox {

	/**
	 * Handling of tied ranks in {@link ToolBox#rank(double[], TiedRanks)}.
	 * 
	 * @author Jan Grau
	 *
	 */
	public enum TiedRanks{
		/**
		 * Identical values obtain different ranks. The ranks are assigned in the order of the input array.
		 */
		IN_ORDER,
		/**
		 * Identical values obtain identical ranks. If multiple values obtain identical ranks, the following rank is incremented by the multiplicity.
		 */
		SPORTS,
		/**
		 * Identical values obtain identical ranks. If multiple values obtain identical ranks, the following rank is incremented by one.
		 */
		CONTIGUOUS
	}
	
	
	/**
	 * This method converts a {@link HashSet} in a {@link Hashtable} with unique indices starting at 0.
	 * The indices are derived from the order of the {@link Iterator} of the {@link HashSet}.
	 * 
	 * @param <K> the type of the keys
	 * 
	 * @param set the set of keys
	 * 
	 * @return a {@link Hashtable} with keys and unique indices starting at 0
	 */
	public static <K> Hashtable<K,Integer> parseHashSet2IndexHashtable( HashSet<K> set ){
		Hashtable<K, Integer> hash = new Hashtable<K, Integer>();
		if( set != null ) {
			Iterator<K> it = set.iterator();
			int i = 0;
			
			while( it.hasNext() ) {
				hash.put( it.next(), i );
				i++;
			}
		}
		return hash;
	}
	
	/**
	 * This method returns the maximum of the elements of an <code>array</code>.
	 * 
	 * @param array
	 *            the array of values
	 * 
	 * @return the maximum
	 */
	public static double max( double... array ) {
		return max( 0, array.length, array );
	}
	
	/**
	 * This method returns the minimum of the elements of an <code>array</code>.
	 * 
	 * @param array
	 *            the array of values
	 * 
	 * @return the minimum
	 */
	public static double min( double... array ) {
		return min( 0, array.length, array );
	}
	
	/**
	 * This method returns the maximum of the elements of an <code>array</code>
	 * between <code>start</code> and <code>end</code>.
	 * 
	 * @param start
	 *            start position (inclusive)
	 * @param end
	 *            end position (exclusive)
	 * @param array
	 *            the array of values
	 * 
	 * @return the maximum
	 */
	public static double max( int start, int end, double[] array ) {
		if( end <= start ) {
			throw new IllegalArgumentException();
		} else {
			double max = array[start];
			start++;
			while( start < end ) {
				if( max < array[start] ) {
					max = array[start];
				}
				start++;
			}
			return max;
		}
	}
	
	/**
	 * This method returns the minimum of the elements of an <code>array</code>
	 * between <code>start</code> and <code>end</code>.
	 * 
	 * @param start
	 *            start position (inclusive)
	 * @param end
	 *            end position (exclusive)
	 * @param array
	 *            the array of values
	 * 
	 * @return the minimum
	 */
	public static double min( int start, int end, double[] array ) {
		if( end <= start ) {
			throw new IllegalArgumentException("The start index has to be smaller than the end index.");
		} else {
			double min = array[start];
			start++;
			while( start < end ) {
				if( min > array[start] ) {
					min = array[start];
				}
				start++;
			}
			return min;
		}
	}
	
	/**
	 * Returns the index with maximal value in a <code>double</code> array.
	 * 
	 * @param w
	 *            the given <code>double</code> array
	 * 
	 * @return the index
	 */
	public static final int getMaxIndex( double[] w ) {
		return getMaxIndex( 0, w.length, w );
	}

	/**
	 * Returns the index with maximal value in a <code>double</code> array.
	 * 
	 * @param w
	 *            the given <code>double</code> array
	 * @param start
	 *            start position (inclusive)
	 * @param end
	 *            end position (exclusive)
	 * 
	 * @return the index
	 */
	public static final int getMaxIndex( int start, int end, double[] w ) {
		if(start>=end){
			return -1;
		}
		int max = start, i = start+1;
		for( ; i < end; i++ ) {
			if( w[i] > w[max] ) {
				max = i;
			}
		}
		return max;
	}
	
	/**
	 * This method returns the median of the elements of <code>array</code>.
	 * 
	 * @param array
	 *            the array of values
	 * 
	 * @return the median
	 */
	public static double median( double... array ) {
		return median( 0, array.length, array );
	}
	
	/**
	 * This method returns the median of the elements of an <code>array</code>
	 * between <code>start</code> and <code>end</code>.
	 * 
	 * @param start
	 *            start position (inclusive)
	 * @param end
	 *            end position (exclusive)
	 * @param array
	 *            the array of values
	 * 
	 * @return the median
	 */
	public static double median( int start, int end, double[] array ) {
		if( end < start ) {
			throw new IllegalArgumentException("The start index has to be smaller than the end index.");
		} else if(end == start){
			return Double.NaN;
		} else {
			double[] ar2 = new double[end-start];
			System.arraycopy( array, start, ar2, 0, ar2.length );
			Arrays.sort( ar2 );
			int idx = ar2.length/2;
			if(ar2.length % 2 == 0){
				return (ar2[idx-1]+ar2[idx])/2.0;
			}else{
				return ar2[idx];
			}
		}
	}
	
	/**
	 * This method returns the mean of the elements of an <code>array</code>
	 * between <code>start</code> and <code>end</code>.
	 * 
	 * @param start
	 *            start position (inclusive)
	 * @param end
	 *            end position (exclusive)
	 * @param array
	 *            the array of values
	 * @return the mean
	 */
	public static double mean( int start, int end, double[] array ) {
		if( start >= end ) {
			throw new IllegalArgumentException( "The start index has to be smaller than the end index." );
		}
		double sum = 0;
		for( int k = start; k < end; k++ ) {
			sum += array[k];
		}
		return sum / (double)( end - start );
	}
	
	/**
	 * This method returns the standard deviation of the elements of an <code>array</code>
	 * between <code>start</code> and <code>end</code>.
	 * 
	 * @param start
	 *            start position (inclusive)
	 * @param end
	 *            end position (exclusive)
	 * @param array
	 *            the array of values
	 * @return the standard deviation
	 */
	public static double sd( int start, int end, double[] array ) {
		double mean = mean( start, end, array );
		double var = 0, h;
		for( int k = start; k < end; k++ ) {
			h = mean - array[k];
			var += h*h;
		}
		return Math.sqrt( var / (double)( end - start ) );
	}
	
	/**
	 * Returns the <code>percent</code> percentile of the array, i.e.,
	 * returns the element at <code>percent*(array.length)</code> of the sorted array.
	 * 
	 * @param array the array
	 * @param percent the percentage
	 * @return the percentile value
	 */
	public static double percentile( double[] array, double percent ) {
		return percentile( 0, array.length, array, percent );
	}
	
	/**
	 * Returns the <code>percent</code> percentile of the array between start and end, i.e.,
	 * returns the element at <code>percent*(end-start)</code> of the sorted sub-array.
	 * 
	 * @param start the first index to consider
	 * @param end the end index (exclusive)
	 * @param array the array
	 * @param percent the percentage
	 * @return the percentile value
	 */
	public static double percentile( int start, int end, double[] array, double percent ) {
		if( end <= start || percent < 0 || percent > 1 ) {
			throw new IllegalArgumentException();
		} else {
			double[] ar2 = new double[end-start];
			System.arraycopy( array, start, ar2, 0, ar2.length );
			Arrays.sort( ar2 );
			int idx = (int)Math.ceil( percent*ar2.length );
			return ar2[idx];
		}
	}
	
	/**
	 * Returns the element in the sorted variant of <code>array</code> (considered in descending order)
	 * that yields <code>percent</code> of the cumulative value of array elements.
	 * 
	 * @param array the array
	 * @param percent the percentage
	 * @return the element
	 */
	public static double weightPercentile(double[] array, double percent){
		double t = ToolBox.sum( array )*percent;
		double[] ar2 = new double[array.length];
		System.arraycopy( array, 0, ar2, 0, array.length );
		Arrays.sort( ar2 );
		int j=ar2.length-1;
		while( j >= 0 && (t - ar2[j]) > 0 ){
			t -= ar2[j];
			j--;
		}
		return ar2[j];
	}
	
	/**
	 * Ranks the values in <code>values</code>, where the greatest value obtains the lowest rank.
	 * The boolean <code>sameRank</code> allows to decide whether tied values should obtain the same rank.
	 * 
	 * @param values the values
	 * @param sameRank a switch whether tied values obtain the same rank.
	 * @return the ranks
	 */
	public static final int[] rank( double[] values, boolean sameRank ){
		return rank(values, sameRank ? TiedRanks.CONTIGUOUS : TiedRanks.IN_ORDER);
	}
	
	/**
	 * This method computes the order of the elements to obtain a sorted array.
	 * 
	 * @param values the array
	 * @param descending a switch whether the array should be sorted ascending or descending
	 * 
	 * @return the array containing the indices of the elements to obtain a sorted array
	 */
	public static final int[] order(double[] values, boolean descending){//TODO FIXME desc? asc?
		ComparableElement<Integer, Double>[] comp = new ComparableElement[values.length];
		for(int i=0;i<values.length;i++){
			comp[i] = new ComparableElement<Integer, Double>(i, descending ? -values[i] : values[i] );
		}
		
		Arrays.sort(comp);
		int[] res = new int[values.length];
		for(int i=0;i<res.length;i++){
			res[ i ] = comp[i].getElement();
		}
		return res;
	}
	
	/**
	 * Ranks the values in <code>values</code>, where the greatest value obtains the lowest rank.
	 * The enum <code>ties</code> allows to choose how tied ranks are handled.
	 * 
	 * @param values the values
	 * @param ties a switch how to handle tied ranks
	 * @return the ranks
	 * @see TiedRanks
	 */
	public static final int[] rank( double[] values, TiedRanks ties ){
		
		int n = values.length;
		double[][] help= new double[n][2];
		for(int i=0;i<n;i++){
			help[i][0] = values[i];
			help[i][1] = i;
		}
		Arrays.sort( help, new DoubleArrayComparator(0) );
		int[] ranks = new int[n];
		n--;
		int rank = 0;
		int temp = 0;
		for(int i=n;i>=0;i--){
			ranks[(int)Math.round(help[i][1])] = rank; 
			if( ties == TiedRanks.IN_ORDER || (i > 0 && help[i][0] != help[i-1][0]) ){
				rank+= temp + 1;
				temp = 0;
			}else if(ties == TiedRanks.CONTIGUOUS){
				temp = 0;
			}else if(ties == TiedRanks.SPORTS){
				temp++;
			}else{
				throw new RuntimeException( ties.name()+" not supported." );
			}
		}
		
		return ranks;
	}
	
	/**
	 * The method computes the Spearman correlation of two vectors.
	 * 	
	 * @param v1 the first vector
	 * @param v2 the second vector
	 * 
	 * @return the Spearman correlation of the two vectors
	 * 
	 * @throws Exception if the vectors have different length
	 */
	public static double spearmanCorrelation(double[] v1, double[] v2) throws Exception{
		if(v1.length != v2.length){
			throw new Exception("Number of values in both vector differ.");
		}
		
		int[] rankTruth = ToolBox.rank( v2, true );
		int[] rankPred = ToolBox.rank( v1, true );
		
		double sumTruth = 0;
		double sumPred = 0;
		double sqTruth = 0;
		double sqPred = 0;
		double cross = 0;
		double n = rankTruth.length;
		
		for(int i=0;i<rankTruth.length;i++){
			sumTruth += rankTruth[i];
			sumPred += rankPred[i];
			sqTruth += rankTruth[i]*rankTruth[i];
			sqPred += rankPred[i]*rankPred[i];
			cross += rankTruth[i]*rankPred[i];
		}
		
		return (cross - sumTruth*sumPred/n)/( Math.sqrt( sqTruth - sumTruth*sumTruth/n )*Math.sqrt( sqPred - sumPred*sumPred/n) );
	}
	
	/**
	 * The method computes the Pearson correlation of two vectors.
	 * 	
	 * @param v1 the first vector
	 * @param v2 the second vector
	 * 
	 * @return the Pearson correlation of the two vectors
	 * 
	 * @throws Exception if the vectors have different length
	 */
	public static double pearsonCorrelation(double[] v1, double[] v2) throws Exception{
		if(v1.length != v2.length){
			throw new Exception("Number of values in both vector differ.");
		}
		
		double sumTruth = 0;
		double sumPred = 0;
		double sqTruth = 0;
		double sqPred = 0;
		double cross = 0;
		double n = v2.length;
		
		for(int i=0;i<v2.length;i++){
			sumTruth += v2[i];
			sumPred += v1[i];
			sqTruth += v2[i]*v2[i];
			sqPred += v1[i]*v1[i];
			cross += v2[i]*v1[i];
		}
		
		return (cross - sumTruth*sumPred/n)/( Math.sqrt( sqTruth - sumTruth*sumTruth/n )*Math.sqrt( sqPred - sumPred*sumPred/n) );
	}

	/**
	 * The method computes the Pearson correlation of two vectors beginning at specific offsets.
	 * 	
	 * @param v1 the first vector
	 * @param v2 the second vector
	 * @param off1 the offset for the first vector
	 * @param off2 the offset for the second vector
	 * 
	 * @return the Pearson correlation of the two vectors
	 * 
	 * @throws Exception if the vectors have different length
	 */
	public static double pearsonCorrelation(double[] v1, double[] v2, int off1, int off2) throws Exception{
		int l1 = v1.length-off1;
		int l2 = v2.length-off2;
		if(l1 < 1 || l2 < 1){
			return 0;
		}
		return pearsonCorrelation( v1, v2, off1, off2, Math.min( l1, l2 ) );
	}
	
	/**
	 * The method computes the Pearson correlation of two vectors beginning at specific offsets and using a given number of entries.
	 * 	
	 * @param v1 the first vector
	 * @param v2 the second vector
	 * @param off1 the offset for the first vector
	 * @param off2 the offset for the second vector
	 * @param length the number of elements to be used for computing the correlation
	 * 
	 * @return the Pearson correlation of the two vectors
	 * 
	 * @throws Exception if the vectors have different length
	 */
	public static double pearsonCorrelation(double[] v1, double[] v2, int off1, int off2, int length) throws Exception{
		if(length < 1 || v1.length < off1+length || v2.length < off2+length){
			throw new Exception("Wrong length");
		}
		
		double sumTruth = 0;
		double sumPred = 0;
		double sqTruth = 0;
		double sqPred = 0;
		double cross = 0;
		double n = length;
		
		for(int i=0;i<length;i++){
			sumTruth += v2[off2+i];
			sumPred += v1[off1+i];
			sqTruth += v2[off2+i]*v2[off2+i];
			sqPred += v1[off1+i]*v1[off1+i];
			cross += v2[off2+i]*v1[off1+i];
		}
		
		return (cross - sumTruth*sumPred/n)/( Math.sqrt( sqTruth - sumTruth*sumTruth/n )*Math.sqrt( sqPred - sumPred*sumPred/n) );
	}
	
	/**
	 * Computes the sum of the values in <code>array</code>
	 * @param array the array
	 * @return the sum of the elements in <code>array</code>
	 */
	public static double sum(double... array){
		return sum(0,array.length,array);
	}
	
	/**
	 * Computes the sum of the values in <code>array</code> starting at
	 * <code>start</code> until <code>end</code>.
	 * @param start the index of the first element in the sum (inclusive)
	 * @param end the end index (exclusive)
	 * @param array the array
	 * @return the sum of the elements in <code>array</code> from <code>start</code> to <code>end</code>
	 */
	public static double sum( int start, int end, double[] array ) {
		double sum = 0.0;
		for(int i=start;i<end;i++){
			sum += array[i];
		}
		return sum;
	}
	
	/**
	 * Transpose a <code>double</code> matrix.
	 * 
	 * @param sp the array
	 * @return the transposed array
	 * 
	 * @throws Exception if <code>sp</code> is no matrix
	 */
	public static double[][] transpose( double[][] sp ) throws Exception {
		double[][] tr = new double[sp[0].length][sp.length];
		for(int i=0;i<sp.length;i++){
			if(sp[i].length != tr.length){
				throw new Exception("Not a matrix.");
			}
			for(int j=0;j<sp[i].length;j++){
				tr[j][i] = sp[i][j];
			}
		}
		return tr;
	}
	
	
	/**
	 * Creates an array of hue values that can be used for the representation
	 * of probabilities. The first four hues are always equal to those of the sequence-logo colors.
	 * @param number the number of colors
	 * @return the hue values
	 */
	public static double[] getUniqueHueValues(int number){
		double[] vals = new double[number];
		for(int i=0;i<vals.length;i++){
			vals[i] = i * 360.0 / number;
			if(vals[i] == 120 || vals[i] == 240 || vals[i] == 40 || vals[i] == 0){
				vals[i] = Double.NaN;
			}
		}
		double[] vals2 = new double[number];
		vals2[0] = 120.0/360.0;
		vals2[1] = 240.0/360.0;
		vals2[2] = 40.0/360.0;
		vals2[3] = 0;
		int curr = (int) Math.round( number/2.0 );
		int off = (int) Math.round( number/4.0 );
		for(int i=4;i<vals2.length;i++){
			if(!Double.isNaN( vals[curr] ) ){
				vals2[i] = vals[curr]/360.0;
				vals[curr] = Double.NaN;
			}else{
				i--;
			}
			if(curr + off < number){
				curr += off;
			}else{
				curr = off;
				off = (int) Math.round( off/2.0 );
			}
		}
		return vals2;
	}

	/**
	 * This method implements a sort algorithm on the array <code>arrayToBeSorted</code>.
	 * However, the entries of the array <code>alongWith</code> are reordered in the same way as<code>arrayToBeSorted</code>.
	 * 
	 * @param arrayToBeSorted the array to be sorted
	 * @param alongWith a series of arrays that is in the same way reordered as <code>arrayToBeSorted</code>.
	 * 
	 * @see Arrays#sort(Object[], Comparator)
	 * @see DoubleArrayComparator
	 */
	public static void sortAlongWith( double[] arrayToBeSorted, double[]... alongWith ) {
		boolean simple = alongWith == null;
		if( !simple ) {
			int i = 0; 
			while( i < alongWith.length && alongWith[i] == null ) {
				i++;
			}
			simple = alongWith.length == i;
		}
		if( simple ) {
			Arrays.sort( arrayToBeSorted );
		} else {
			double[][] h = new double[arrayToBeSorted.length][1+alongWith.length];
			for( int i = 0; i < h.length; i++ ) {
				h[i][0] = arrayToBeSorted[i];
				for( int j = 1; j <= alongWith.length; j++ ) {
					h[i][j] = alongWith[j-1][i];
				}
			}
			Arrays.sort( h, new DoubleArrayComparator(0) );
			for( int i = 0; i < h.length; i++ ) {
				arrayToBeSorted[i] = h[i][0];
				for( int j = 1; j <= alongWith.length; j++ ) {
					alongWith[j-1][i] = h[i][j];
				}
			}
		}		
	}
	
	/**
	 * This methods returns a {@link String} representation of a double array using the specified {@link NumberFormat}.
	 * 
	 * @param array the array to be represented
	 * @param nf the {@link NumberFormat} to be used
	 * 
	 * @return a {@link String} representation
	 */
	public static String toString( double[] array, NumberFormat nf ) {
		if( array == null ) {
			return "null";
		}
		StringBuffer sb = new StringBuffer();
		sb.append("[");
		for( int i = 0; i < array.length; i++ ) {
			sb.append( (i==0?"":", ") + nf.format(array[i]) );
		}
		sb.append("]");
		return sb.toString();
	}
	
	/**
	 * Creates and returns a permutation of the values in <code>values</code>.
	 * @param values the values
	 * @return the permutation
	 */
	public static double[] permute(double[] values){
		return permute(values,null);
	}
	
	/**
	 * Creates a new permutation of the values in <code>values</code>. If the array supplied as
	 * <code>permutation</code> has the same dimension as <code>values</code> and is not <code>null</code>,
	 * the values are returned in this array (and, additionally, as return value), which may save memory
	 * in multiple successive calls to this method.
	 * @param values the values
	 * @param permutation the permutation (may be <code>null</code>)
	 * @return the permutation
	 */
	public static double[] permute(double[] values, double[] permutation){
		Random r = new Random();
		if(permutation == null || permutation.length != values.length){
			permutation = new double[values.length];
		}
		for( int i = 0; i < permutation.length; i++ ) {
			permutation[i] = values[i];
		}
		int temp;
		double temp2;
		for( int i = permutation.length - 1; i >= 0; i-- ) {
			temp = r.nextInt( i + 1 );
			temp2 = permutation[i];
			permutation[i] = permutation[temp];
			permutation[temp] = temp2;
		}
		return permutation;
	}

	/**
	 * Counts the number of <code>true</code> values in <code>bools</code> (similar to sum on booleans in R).
	 * @param bools the boolean values
	 * @return the number of <code>true</code>s
	 */
	public static int sum( boolean[] bools ) {
		int num = 0;
		for(int i=0;i<bools.length;i++){
			if(bools[i]){
				num++;
			}
		}
		return num;
	}

	/**
	 * Returns the index with minimum value in a <code>double</code> array.
	 * 
	 * @param w
	 *            the given <code>double</code> array
	 * 
	 * @return the index
	 */
	public static int getMinIndex(double[] w){
		return getMinIndex(0,w.length,w);
	}
	
	/**
	 * Returns the index with minimum value in a <code>double</code> array.
	 * 
	 * @param w
	 *            the given <code>double</code> array
	 * @param start
	 *            start position (inclusive)
	 * @param end
	 *            end position (exclusive)
	 * 
	 * @return the index
	 */
	public static int getMinIndex(int start, int end, double[] w) {
		if(start>=end){
			return -1;
		}
		int min = start, i = start+1;
		for( ; i < end; i++ ) {
			if( w[i] < w[min] ) {
				min = i;
			}
		}
		return min;
	}
	
	
}
