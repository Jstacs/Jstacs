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

package de.jstacs.classifiers.utils;

/**
 * This class can be used to compute any p-value from a given statistic.
 * Therefore the array of scores has to be sorted.
 * 
 * @author Jens Keilwagen
 * 
 * @see java.util.Arrays#sort(double[])
 */
public class PValueComputation {

	/**
	 * This method searches for the insertion point of the score in a given
	 * sorted array of scores and returns the p-value for this score.
	 * 
	 * @param sortedScores
	 *            the array of sorted scores
	 * @param myScore
	 *            the score that is searched in the array and that corresponds
	 *            to the p-value
	 * 
	 * @return the p-value for the score
	 * 
	 * @see java.util.Arrays#sort(double[])
	 * @see PValueComputation#getIndex(double[], double)
	 */
	public final static double getPValue( double[] sortedScores, double myScore ) {
		return 1d - getIndex( sortedScores, myScore, 0 ) / (double)( sortedScores.length );
	}

	/**
	 * This method searches for the insertion point of the score in a given
	 * sorted array of scores from index <code>start</code> and returns the
	 * p-value for this score.
	 * 
	 * @param sortedScores
	 *            the array of sorted scores
	 * @param myScore
	 *            the score that is searched in the array and that corresponds
	 *            to the p-value
	 * @param start
	 *            the start position in the array where to start the search
	 *            (accelerates the search)
	 * 
	 * @return the p-value for the score
	 * 
	 * @see java.util.Arrays#sort(double[])
	 * @see PValueComputation#getIndex(double[], double, int)
	 */
	public final static double getPValue( double[] sortedScores, double myScore, int start ) {
		return 1d - getIndex( sortedScores, myScore, start ) / (double)( sortedScores.length );
	}

	/**
	 * This method searches in <code>sortedScores</code> for the index
	 * <code>i</code> so that
	 * <code>sortedScores[i-1] &lt; myScore &lt;= sortedScores[i]</code>. This
	 * index could be seen as insertion point for <code>myScore</code>.
	 * 
	 * @param sortedScores
	 *            the array of sorted scores
	 * @param myScore
	 *            the current score
	 * 
	 * @return the index where <code>myScore</code> should be inserted
	 * 
	 * @see java.util.Arrays#sort(double[])
	 * @see PValueComputation#getIndex(double[], double, int)
	 */
	public final static int getIndex( double[] sortedScores, double myScore ) {
		return getIndex( sortedScores, myScore, 0 );
	}

	/**
	 * This method searches in <code>sortedScores</code> beginning at
	 * <code>start</code> for the index <code>i</code> so that
	 * <code>sortedScores[i-1] &lt; myScore &lt;= sortedScores[i]</code>. This
	 * index could be seen as insertion point for <code>myScore</code>.
	 * 
	 * @param sortedScores
	 *            the array of sorted scores
	 * @param myScore
	 *            the current score
	 * @param start
	 *            the start index (inclusive)
	 * 
	 * @return the index where <code>myScore</code> should be inserted
	 * 
	 * @see java.util.Arrays#sort(double[])
	 */
	public final static int getIndex( double[] sortedScores, double myScore, int start ) {
		if( start >= sortedScores.length ) {
			return sortedScores.length;
		} else if( start < 0 ) {
			return 0;
		} else if( myScore <= sortedScores[start] ) {
			return start;
		} else {
			int bound, end = sortedScores.length;
			// binary search between start and end
			// sortedScores[start] < myScore <= "sortedScores[end]"
			do {
				bound = start + ( end - start ) / 2;
				if( sortedScores[bound] < myScore ) {
					start = bound;
				} else {
					end = bound;
				}
			} while( end - start > 1 );
			return end;
		}
	}

	/**
	 * This method finds the first index that has a significant score.
	 * 
	 * @param sortedScores
	 *            the array of sorted scores
	 * @param pValue
	 *            the p-value that is used to determine what is significant and
	 *            what not (the smaller, the more significant)
	 * 
	 * @return the index of the first significant score
	 * 
	 * @see java.util.Arrays#sort(double[])
	 */
	public final static int getBorder( double[] sortedScores, double pValue ) {
		int index = (int)Math.ceil( sortedScores.length * ( 1d - pValue ) );
		if( index == 0 || index == sortedScores.length || sortedScores[index - 1] < sortedScores[index] ) {
			return index;
		} else {
			int bound, start = index, end = sortedScores.length;
			// sortedScores[index] == sortedScores[start] < "sortedScores[end]"
			while( end - start > 1 ) {
				bound = start + ( end - start ) / 2;
				if( sortedScores[start] < sortedScores[bound] ) {
					end = bound;
				} else { // sortedScores[start] == sortedScores[bound]
					start = bound;
				}
			}
			return end;
		}
	}

	/**
	 * This method returns the threshold that determines if an observed score is
	 * significant. The observed score is significant if it is bigger than the
	 * threshold.
	 * 
	 * @param sortedScores
	 *            the array of sorted scores
	 * @param signIndex
	 *            the index of the first significant score
	 * 
	 * @return the threshold
	 */
	public final static double getThreshold( double[] sortedScores, int signIndex ) {
		if( signIndex < 0 ) {
			throw new IndexOutOfBoundsException();
		} else if( signIndex >= sortedScores.length ) {
			return sortedScores[sortedScores.length-1] * ( 1d + 1E-15 * Math.signum( sortedScores[sortedScores.length-1] ) )
				+ Double.MIN_VALUE;
		} else {
			return sortedScores[signIndex];
		}
	}
}
