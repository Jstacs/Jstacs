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

package de.jstacs.utils;

/**
 * This class can be used for normalisation of any <code>double</code> array or
 * a part of a <code>double</code> array.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class Normalisation {

	/**
	 * Returns the logarithm of the sum of values <code>val[i]</code> given as
	 * <code>lnVal[i] = Math.log( val[i] )</code>.
	 * 
	 * @param lnVal
	 *            the logs of the values, i.e.
	 *            <code>lnVal[i] = Math.log( val[i] )</code>
	 * 
	 * @return the logarithm of the sum of values
	 * 			{@latex.inline $\\log(\\sum_i \\mathrm{val}[i])$}
	 * 
	 * @see Normalisation#getLogSum(int, int, double...)
	 * 
	 */
	public static double getLogSum( double... lnVal ) {
		return getLogSum( 0, lnVal.length, lnVal );
	}

	/**
	 * Returns the logarithm of the sum of values <code>v[i]</code> given as
	 * <code>lnVal[i] = Math.log( val[i] )</code> between a start and end index.
	 * 
	 * @param start
	 *            the first index in <code>lnVal</code> considered for the sum
	 * @param end
	 *            the index after the last index considered for the sum
	 * @param lnVal
	 *            the logs of the values, i.e.
	 *            <code>lnVal[i] = Math.log( val[i] )</code>
	 * 
	 * @return the logarithm of the sum of values between the start and end
	 *         index {@latex.inline $\\log(\\sum_{i=start}^{end - 1} \\mathrm{val}[i])$}
	 * 
	 */
	public static double getLogSum( int start, int end, double... lnVal ) {
		/*
		//default
		double logSum = lnVal[0];
		if( lnVal.length > 1 )
		{
			for( int i = 1; i < lnVal.length; i++ )
			{
				logSum = getLogSumAB( logSum, lnVal[i] );
			}
		}
		return logSum;
		*/
		// same proceeding as in normalisation-methods
		double offset = Double.NEGATIVE_INFINITY, sum = 0;
		int i;
		for( i = start; i < end; i++ ) {
			offset = Math.max( offset, lnVal[i] );

		}
		if( Double.isInfinite( offset ) ) {
			return Double.NEGATIVE_INFINITY;
		} else {
			for( i = start; i < end; i++ ) {
				sum += Math.exp( lnVal[i] - offset );

			}
			return offset + Math.log( sum );
		}

	}

	/**
	 * The method does a log-sum-normalisation on the array <code>d</code>, where
	 * the values of <code>d</code> are assumed to be logarithmised.
	 * Let {@latex.inline $d_i = \\log(v_i) \\Leftrightarrow v_i = \\exp(d_i)$} and {@latex.inline $s := \\sum_{i=0}^{\\mathrm{length}(d)-1} v_i$}.
	 * Then after log-sum-normalisation, the array <code>d</code> contains the normalized original values, i.e., {@latex.inline $d_i$} is set to
	 * {@latex.inline $d_i := \\frac{v_i}{s}$ }. The method returns the log-sum of the values, {@latex.inline $\\log(s)$ }.
	 * 
	 * @param d
	 *            the array with the logarithmised values that should be
	 *            normalised
	 * 
	 * @return the logarithm of the sum of the values
	 *         {@latex.inline $\\log(\\sum_{i=0}^{\\mathrm{length}(d)-1} v_i)$ }
	 * 
	 * @see Normalisation#logSumNormalisation(double[], int, int, double[], int)
	 */
	public static double logSumNormalisation( double[] d ) {
		return logSumNormalisation( d, 0, d.length, d, 0 );
	}

	/**
	 * The method does a log-sum-normalisation on the values of the array <code>d</code> between start
	 * index <code>startD</code> and end index <code>endD</code>, where
	 * the values of <code>d</code> are assumed to be logarithmised.
	 * Let {@latex.inline $d_i = \\log(v_i) \\Leftrightarrow v_i = \\exp(d_i)$} and {@latex.inline $s := \\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i$}.
	 * Then after log-sum-normalisation, the part of the array <code>d</code> between start
	 * index <code>startD</code> and end index <code>endD</code> contains the normalized original values, i.e., {@latex.inline $\\forall i=\\mathrm{startD},\\ldots,\\mathrm{endD}\\ d_i$} is set to
	 * {@latex.inline $d_i := \\frac{v_i}{s}$ }. The method returns the log-sum of the values, {@latex.inline $\\log(s)$ }.
	 * 
	 * @param d
	 *            the array with the logarithms of the values that should be
	 *            normalised
	 * @param startD
	 *            the first index in <code>d</code> considered for the
	 *            log-sum-normalisation
	 * @param endD
	 *            the index after the last index in <code>d</code> considered
	 *            for the log-sum-normalisation
	 * 
	 * @return the logarithm of the sum of the values between
	 *         <code>startD</code> and <code>endD</code>
	 *         {@latex.inline $\\log(\\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} val[i])$ }
	 * 
	 * @see Normalisation#logSumNormalisation(double[], int, int, double[], int)
	 */
	public static double logSumNormalisation( double[] d, int startD, int endD ) {
		return logSumNormalisation( d, startD, endD, d, startD );
	}

	/**
	 * The method does a log-sum-normalisation on the values of the array <code>d</code> between start
	 * index <code>startD</code> and end index <code>endD</code>, where
	 * the values of <code>d</code> are assumed to be logarithmised. In addition to <code>d</code> another array of values 
	 * <code>secondValues</code> is considered for the normalization constant, but not normalized itself.
	 * 
	 * Let {@latex.inline $d_i = \\log(v_i) \\Leftrightarrow v_i = \\exp(d_i)$} and {@latex.inline $s := \\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i + \\sum_{i=0}^{\\mathrm{length(secondValues)}-1} \\exp(\\mathrm{secondValues}_i)$}.
	 * Then after log-sum-normalisation, the part of the array <code>d</code> starting at
	 * index <code>startD</code> contains the normalized original values, i.e., {@latex.inline $\\forall i=\\mathrm{startD},\\ldots,\\mathrm{endD}-1\\ d_i$} is set to
	 * {@latex.inline $d_i := \\frac{v_i}{s}$ }. The method returns the log-sum of the values, {@latex.inline $\\log(s)$ }.
	 * The method overwrites the values of <code>d</code>
	 * starting at position <code>startD</code>!. 
	 * <code>secondValues</code> will be changed during
	 * log-sum-normalisation and will not be written to <code>d</code>.
	 * 
	 * @param d
	 *            the array with the logarithmised values that should be
	 *            normalised
	 * @param startD
	 *            the first index in <code>d</code> considered for the
	 *            log-sum-normalisation
	 * @param endD
	 *            the index after the last index in <code>d</code> considered
	 *            for the log-sum-normalisation
	 * @param secondValues
	 *            second array with additional values, the whole array is
	 *            considered for the log-sum-normalisation
	 * 
	 * @return the logarithm of the sum of the values of <code>d</code> between
	 *         <code>startD</code> and <code>endD</code> and the values of
	 *         <code>secondValue</code><br>
	 *         {@latex.inline $\\log(\\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i + \\sum_{i=0}^{\\mathrm{length(secondValues)}-1} \\exp(\\mathrm{secondValues}_i))$ }
	 */
	public static double logSumNormalisation( double[] d, int startD, int endD, double[] secondValues ) {
		return logSumNormalisation( d, startD, endD, secondValues, d, startD );
	}

	/**
	 * The method does a log-sum-normalisation on the values of the array <code>d</code> between start
	 * index <code>startD</code> and end index <code>endD</code>, where
	 * the values of <code>d</code> are assumed to be logarithmised.
	 * Let {@latex.inline $d_i = \\log(v_i) \\Leftrightarrow v_i = \\exp(d_i)$} and {@latex.inline $s := \\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i$}.
	 * Then after log-sum-normalisation, the part of the array <code>dest</code> starting at
	 * index <code>startDest</code> contains the normalized original values, i.e., {@latex.inline $\\forall i=\\mathrm{startDest},\\ldots,\\mathrm{startDest}+(\\mathrm{endD}-\\mathrm{startD})-1\\ \\mathrm{dest}_i$} is set to
	 * {@latex.inline $\\mathrm{dest}_i := \\frac{v_j}{s}$ } where {@latex.inline $j = i - \\mathrm{startDest} + \\mathrm{startD}$}. The method returns the log-sum of the values, {@latex.inline $\\log(s)$ }.
	 * The method writes the result of <code>d</code> in <code>dest</code>
	 * starting at position <code>startDest</code> while <code>d</code> remains
	 * unchanged. <code>secondValues</code> will be changed during
	 * log-sum-normalisation and will not be written to <code>dest</code>.
	 * 
	 * @param d
	 *            the array with the logarithmised values that should be
	 *            normalised
	 * @param startD
	 *            the first index in <code>d</code> considered for the
	 *            log-sum-normalisation
	 * @param endD
	 *            the index after the last index in <code>d</code> considered
	 *            for the log-sum-normalisation
	 * @param dest
	 *            the destination array for the normalised values
	 * @param startDest
	 *            the start index of the destination array
	 * 
	 * @return the logarithm of the sum of the values of <code>d</code> between
	 *         <code>startD</code> and <code>endD</code>
	 *         {@latex.inline $\\log(\\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i )$ }
	 */
	public static double logSumNormalisation( double[] d, int startD, int endD, double[] dest, int startDest ) {
		return logSumNormalisation( d, startD, endD, null, dest, startDest );
	}

	/**
	 * The method does a log-sum-normalisation on the values of the array <code>d</code> between start
	 * index <code>startD</code> and end index <code>endD</code>, where
	 * the values of <code>d</code> are assumed to be logarithmised. In addition to <code>d</code> another array of values 
	 * <code>secondValues</code> is considered for the normalization constant, but not normalized itself.
	 * 
	 * Let {@latex.inline $d_i = \\log(v_i) \\Leftrightarrow v_i = \\exp(d_i)$} and {@latex.inline $s := \\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i + \\sum_{i=0}^{\\mathrm{length(secondValues)}-1} \\exp(\\mathrm{secondValues}_i)$}.
	 * Then after log-sum-normalisation, the part of the array <code>dest</code> starting at
	 * index <code>startDest</code> contains the normalized original values, i.e., {@latex.inline $\\forall i=\\mathrm{startDest},\\ldots,\\mathrm{startDest}+(\\mathrm{endD}-\\mathrm{startD})-1\\ \\mathrm{dest}_i$} is set to
	 * {@latex.inline $\\mathrm{dest}_i := \\frac{v_j}{s}$ } where {@latex.inline $j = i - \\mathrm{startDest} + \\mathrm{startD}$}. The method returns the log-sum of the values, {@latex.inline $\\log(s)$ }.
	 * The method writes the result of <code>d</code> in <code>dest</code>
	 * starting at position <code>startDest</code> while <code>d</code> remains
	 * unchanged. <code>secondValues</code> will be changed during
	 * log-sum-normalisation and will not be written to <code>dest</code>.
	 * 
	 * @param d
	 *            the array with the logarithmised values that should be
	 *            normalised
	 * @param startD
	 *            the first index in <code>d</code> considered for the
	 *            log-sum-normalisation
	 * @param endD
	 *            the index after the last index in <code>d</code> considered
	 *            for the log-sum-normalisation
	 * @param secondValues
	 *            second array with additional values, the whole array is
	 *            considered for the log-sum-normalisation
	 * @param dest
	 *            the destination array for the normalised values
	 * @param startDest
	 *            the start index of the destination array
	 * 
	 * @return the logarithm of the sum of the values of <code>d</code> between
	 *         <code>startD</code> and <code>endD</code> and the values of
	 *         <code>secondValue</code><br>
	 *         {@latex.inline $\\log(\\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i + \\sum_{i=0}^{\\mathrm{length(secondValues)}-1} \\exp(\\mathrm{secondValues}_i))$ }
	 */
	public static double logSumNormalisation( double[] d, int startD, int endD, double[] secondValues, double[] dest, int startDest ) {
		double offset = Double.NEGATIVE_INFINITY;
		int i = startD;
		for( ; i < endD; i++ ) {
			offset = Math.max( offset, d[i] );
		}
		if( secondValues != null ) {
			for( i = 0; i < secondValues.length; i++ ) {
				offset = Math.max( offset, secondValues[i] );
			}
		}
		return logSumNormalisation( d, startD, endD, offset, secondValues, dest, startDest );
	}

	
	/**
	 * The method does a log-sum-normalisation on the array <code>d</code>, where
	 * the values of <code>d</code> are assumed to be logarithmised.
	 * Let {@latex.inline $d_i = \\log(v_i) \\Leftrightarrow v_i = \\exp(d_i)$} and {@latex.inline $s := \\sum_{i=0}^{\\mathrm{length}(d)-1} v_i$}.
	 * Then after log-sum-normalisation, the array <code>d</code> contains the normalized original values, i.e., {@latex.inline $d_i$} is set to
	 * {@latex.inline $d_i := \\frac{v_i}{s}$ }. The method returns the log-sum of the values, {@latex.inline $\\log(s)$ }.
	 * 
	 * @param d
	 *            the array with the logarithmised values that should be
	 *            normalised
	 * @param offset
	 *            the offset on the log-values which is used to get more accurate results in
	 *            the normalization. Typically, this is set to the maximum of the log-values.
	 * 
	 * @return the logarithm of the sum of the values
	 *         {@latex.inline $\\log(\\sum_{i=0}^{\\mathrm{length}(d)-1} v_i)$ }
	 * 
	 * @see Normalisation#logSumNormalisation(double[], int, int, double[], int)
	 */
	public static double logSumNormalisation( double[] d, double offset ) {
		return logSumNormalisation( d, 0, d.length, offset, d, 0 );
	}

	/**
	 * The method does a log-sum-normalisation on the values of the array <code>d</code> between start
	 * index <code>startD</code> and end index <code>endD</code>, where
	 * the values of <code>d</code> are assumed to be logarithmised.
	 * Let {@latex.inline $d_i = \\log(v_i) \\Leftrightarrow v_i = \\exp(d_i)$} and {@latex.inline $s := \\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i$}.
	 * Then after log-sum-normalisation, the part of the array <code>dest</code> starting at
	 * index <code>startDest</code> contains the normalized original values, i.e., {@latex.inline $\\forall i=\\mathrm{startDest},\\ldots,\\mathrm{startDest}+(\\mathrm{endD}-\\mathrm{startD})-1\\ \\mathrm{dest}_i$} is set to
	 * {@latex.inline $\\mathrm{dest}_i := \\frac{v_j}{s}$ } where {@latex.inline $j = i - \\mathrm{startDest} + \\mathrm{startD}$}. The method returns the log-sum of the values, {@latex.inline $\\log(s)$ }.
	 * The method writes the result of <code>d</code> in <code>dest</code>
	 * starting at position <code>startDest</code> while <code>d</code> remains
	 * unchanged. <code>secondValues</code> will be changed during
	 * log-sum-normalisation and will not be written to <code>dest</code>.
	 * 
	 * @param d
	 *            the array with the logarithmised values that should be
	 *            normalised
	 * @param startD
	 *            the first index in <code>d</code> considered for the
	 *            log-sum-normalisation
	 * @param endD
	 *            the index after the last index in <code>d</code> considered
	 *            for the log-sum-normalisation
	 * @param offset
	 *            the offset on the log-values which is used to get more accurate results in
	 *            the normalization. Typically, this is set to the maximum of the log-values.
	 * @param dest
	 *            the destination array for the normalised values
	 * @param startDest
	 *            the start index of the destination array
	 * 
	 * @return the logarithm of the sum of the values of <code>d</code> between
	 *         <code>startD</code> and <code>endD</code>
	 *         {@latex.inline $\\log(\\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i )$ }
	 */
	public static double logSumNormalisation( double[] d, int startD, int endD, double offset, double[] dest, int startDest ) {
		return logSumNormalisation( d, startD, endD, offset, null, dest, startDest );
	}

	/**
	 * The method does a log-sum-normalisation on the values of the array <code>d</code> between start
	 * index <code>startD</code> and end index <code>endD</code>, where
	 * the values of <code>d</code> are assumed to be logarithmised. In addition to <code>d</code> another array of values 
	 * <code>secondValues</code> is considered for the normalization constant, but not normalized itself.
	 * 
	 * Let {@latex.inline $d_i = \\log(v_i) \\Leftrightarrow v_i = \\exp(d_i)$} and {@latex.inline $s := \\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i + \\sum_{i=0}^{\\mathrm{length(secondValues)}-1} \\exp(\\mathrm{secondValues}_i)$}.
	 * Then after log-sum-normalisation, the part of the array <code>dest</code> starting at
	 * index <code>startDest</code> contains the normalized original values, i.e., {@latex.inline $\\forall i=\\mathrm{startDest},\\ldots,\\mathrm{startDest}+(\\mathrm{endD}-\\mathrm{startD})-1\\ \\mathrm{dest}_i$} is set to
	 * {@latex.inline $\\mathrm{dest}_i := \\frac{v_j}{s}$ } where {@latex.inline $j = i - \\mathrm{startDest} + \\mathrm{startD}$}. The method returns the log-sum of the values, {@latex.inline $\\log(s)$ }.
	 * The method writes the result of <code>d</code> in <code>dest</code>
	 * starting at position <code>startDest</code> while <code>d</code> remains
	 * unchanged. <code>secondValues</code> will be changed during
	 * log-sum-normalisation and will not be written to <code>dest</code>.
	 * 
	 * @param d
	 *            the array with the logarithmised values that should be
	 *            normalised
	 * @param startD
	 *            the first index in <code>d</code> considered for the
	 *            log-sum-normalisation
	 * @param endD
	 *            the index after the last index in <code>d</code> considered
	 *            for the log-sum-normalisation
	 * @param offset
	 *            the offset on the log-values which is used to get more accurate results in
	 *            the normalization. Typically, this is set to the maximum of the log-values.
	 * @param secondValues
	 *            second array with additional values, the whole array is
	 *            considered for the log-sum-normalisation
	 * @param dest
	 *            the destination array for the normalised values
	 * @param startDest
	 *            the start index of the destination array
	 * 
	 * @return the logarithm of the sum of the values of <code>d</code> between
	 *         <code>startD</code> and <code>endD</code> and the values of
	 *         <code>secondValue</code><br>
	 *         {@latex.inline $\\log(\\sum_{i=\\mathrm{startD}}^{\\mathrm{endD}-1} v_i + \\sum_{i=0}^{\\mathrm{length(secondValues)}-1} \\exp(\\mathrm{secondValues}_i))$ }
	 */
	public static double logSumNormalisation( double[] d, int startD, int endD, double offset, double[] secondValues, double[] dest,
			int startDest ) {
		double sum = 0;
		int i = 0, l = endD - startD;
		for( ; i < l; i++ ) {
			dest[startDest + i] = Math.exp( d[startD + i] - offset );

			sum += dest[startDest + i];

		}
		if( secondValues != null ) {
			for( i = 0; i < secondValues.length; i++ ) {
				secondValues[i] = Math.exp( secondValues[i] - offset );
				sum += secondValues[i];
			}
		}
		if( sum != 1d ) {
			normalisation( dest, sum, startDest, startDest + l );
			if( secondValues != null ) {
				normalisation( secondValues, sum );
			}
		}

		return offset + Math.log( sum );
	}

	/**
	 * The method does a sum-normalisation on <code>d</code>, i.e. divides all values
	 * in <code>d</code> by the sum over all values in <code>d</code> and returns the
	 * sum of the values.
	 * 
	 * @param d
	 *            the array with the values that should be normalised
	 * 
	 * @return the sum of the values of <code>d</code>
	 *         {@latex.inline $\\sum_{i=0}^{\\mathrm{length}(d)-1} d[i]$}
	 * 
	 * @see Normalisation#sumNormalisation(double[], double[], int)
	 */
	public static double sumNormalisation( double[] d ) {
		return sumNormalisation( d, d, 0 );
	}

	/**
	 * The method does a sum-normalisation on <code>d</code>, i.e. divides all values
	 * in <code>d</code> by the sum over all values in <code>d</code>
	 * and writes the result to <code>dest</code> starting at position
	 * <code>start</code> while <code>d</code> remains unchanged. The sum of the
	 * values of <code>d</code> will be returned.
	 * 
	 * @param d
	 *            the array with the values that should be normalised
	 * @param dest
	 *            the destination array for the normalised values
	 * @param start
	 *            the start index of the destination array
	 * 
	 * @return the sum of the values of <code>d</code>
	 *         {@latex.inline $\\sum_{i=0}^{\\mathrm{length}(d)-1} d[i] $}
	 */
	public static double sumNormalisation( double[] d, double[] dest, int start ) {
		int i;
		double sum = d[0];
		for( i = 1; i < d.length; i++ ) {
			sum += d[i];
		}
		normalisation( d, sum, dest, start );
		return sum;
	}

	/**
	 * The method does a normalisation on <code>d</code> using the value
	 * <code>v</code> for normalisation.
	 * 
	 * @param d
	 *            the array with the values that should be normalised
	 * @param v
	 *            the value for the normalisation
	 * 
	 * @see Normalisation#normalisation(double[], double, double[], int)
	 */
	public static void normalisation( double[] d, double v ) {
		normalisation( d, v, d, 0 );
	}

	/**
	 * The method does a normalisation on <code>d</code> writing the result to
	 * <code>dest</code> starting at position <code>start</code> while
	 * <code>d</code> remains unchanged. The value <code>v</code> is used for
	 * the normalisation.
	 * 
	 * @param d
	 *            the array with the values that should be normalised
	 * @param v
	 *            the value for normalisation
	 * @param dest
	 *            the destination array for the normalised values
	 * @param start
	 *            the start index of the destination array
	 */
	public static void normalisation( double[] d, double v, double[] dest, int start ) {
		for( int i = 0; i < d.length; i++, start++ ) {
			dest[start] = d[i] / v;
		}
	}

	/**
	 * The method does a sum normalisation on <code>d</code> between start index
	 * <code>start</code> and end index <code>end</code> using the value
	 * <code>v</code> for the normalisation.
	 * 
	 * @param d
	 *            the array with the values that should be normalised
	 * @param v
	 *            the value for normalisation
	 * @param start
	 *            the first index in <code>d</code> considered for the
	 *            log-sum-normalisation
	 * @param end
	 *            the index after the last index in <code>d</code> considered
	 *            for the log-sum-normalisation
	 */
	public static void normalisation( double[] d, double v, int start, int end ) {
		for( ; start < end; start++ ) {
			d[start] /= v;
		}
	}
}
