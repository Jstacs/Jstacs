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

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;

/**
 * This class is a collection of methods which might be useful for the programmer.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class ToolBox {

	
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
			throw new IllegalArgumentException();
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
		int max = start, i = start+1;
		for( ; i < end; i++ ) {
			if( w[i] > w[max] ) {
				max = i;
			}
		}
		return max;
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
		
		int n = values.length;
		double[][] help= new double[n][2];
		for(int i=0;i<n;i++){
			help[i][0] = values[i];
			help[i][1] = i;
		}
		Arrays.sort( help, new DoubleArrayComparator() );
		int[] ranks = new int[n];
		n--;
		int rank = 0;
		for(int i=n;i>=0;i--){
			ranks[(int)Math.round(help[i][1])] = rank; 
			if( !sameRank || (i > 0 && help[i][0] != help[i-1][0]) ){
				rank++;
			}
		}
		/*
		double[] v = values.clone();
		int[] ranks = new int[v.length];
		int rank = -1;
		double last = Double.NEGATIVE_INFINITY;
		int maxIdx = 0;
		for(int i=0;i<v.length;i++){
			maxIdx = getMaxIndex( v );
			if(v[maxIdx] != last){
				rank++;
			}
			ranks[maxIdx] = rank;
			last = v[maxIdx];
			v[maxIdx] = Double.NEGATIVE_INFINITY;
		}
		*/
		return ranks;
	}
	
	
	private static class DoubleArrayComparator implements Comparator<double[]> {
		public int compare(double[] o1, double[] o2) {
			return (int) Math.signum( o1[0] - o2[0] );
		}
		
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
	
}
