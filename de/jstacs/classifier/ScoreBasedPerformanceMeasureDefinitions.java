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

package de.jstacs.classifier;

import java.util.AbstractList;
import java.util.Arrays;

import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.RangeParameter;

/**
 * This class contains the methods that are needed to evaluate a score based
 * 2-class-classifier.
 * 
 * @author Jens Keilwagen, Jan Grau
 * 
 * @see AbstractScoreBasedClassifier
 * @see AbstractClassifier#evaluate(MeasureParameters, boolean, de.jstacs.data.Sample...)
 * @see AbstractClassifier#evaluateAll(MeasureParameters, boolean, de.jstacs.data.Sample...)
 * @see MeasureParameters
 */
public class ScoreBasedPerformanceMeasureDefinitions {

	/**
	 * This method allows to compute a partial ROC curve.
	 * 
	 * @param scoresClass0
	 *            the scores of class 0 (have to be sorted)
	 * @param scoresClass1
	 *            the scores of class 1 (have to be sorted)
	 * @param specs
	 *            these specificities are used as grid for the partial curve
	 * 
	 * @return a matrix
	 *         <ul>
	 *         <li> <code>m[0][i] i</code>-th specificity
	 *         <li> <code>m[1][i] i</code>-th sensitivity
	 *         <li> <code>m[2][i] i</code>-th threshold
	 *         </ul>
	 * 
	 * @throws ParameterException
	 *             if the method {@link RangeParameter#next()} throws one
	 */
	public static double[][] getPartialROC( double[] scoresClass0, double[] scoresClass1, RangeParameter specs ) throws ParameterException {
		int i = 0, fn = 0, n = specs.getNumberOfValues();
		double[][] val = new double[3][n];

		specs.resetToFirst();

		for( i = 0; i < n; i++ ) {
			val[0][i] = (Double)specs.getValue();
			specs.next();
		}
		Arrays.sort( val[0] );

		for( i = 0; i < n; i++ ) {

			val[2][i] = scoresClass1[(int)Math.ceil( val[0][i] * ( scoresClass1.length - 1 ) )];
			while( fn < scoresClass0.length && scoresClass0[fn] <= val[2][i] ) {
				fn++;
			}
			val[1][i] = ( (double)( scoresClass0.length - fn ) ) / (double)scoresClass0.length;

		}

		return val;
	}

	/**
	 * This method computes the area under the precision recall curve. It also
	 * enables the user to compute the complete curve so that it can be used for
	 * plotting. If you like to get the curve, you should use an empty
	 * {@link AbstractList} for <code>list</code>, otherwise use <code>null</code>.
	 * 
	 * @param scoresClass0
	 *            the scores of class 0 (have to be sorted)
	 * @param scoresClass1
	 *            the scores of class 1 (have to be sorted)
	 * @param list
	 *            <code>null</code> or an empty list
	 * 
	 * @return the area under the precision recall curve
	 */
	public static double getAUC_PR( double[] scoresClass0, double[] scoresClass1, AbstractList<double[]> list ) {
		int i_old = 0, j_old = 0, i = 0, j = 0, d = scoresClass1.length, m = scoresClass0.length;
		double helpJ, propTerm;
		double erg = 0, help1 = 0, help2 = 0;

		// find correct start point
		while( ( j < d ) && ( scoresClass0[i] > scoresClass1[j] ) ) {
			j++;
		}
		//i is zero, so p[0] = 1 ...
		double[] p = new double[]{ ( m - i ) / (double)m, ( m - i ) / ( (double)( m - i + d - j ) ) };
		if( list != null ) {
			list.add( p.clone() );
		}

		// which class defines the threshold
		boolean unique, fromMotif = false;
		if( j < d && scoresClass0[i] == scoresClass1[j] ) {
			unique = false;
		} else {
			unique = true;
			fromMotif = true;
		}
		while( i < m && j < d ) {
			i_old = i;
			j_old = j;
			// find next possible threshold
			if( unique ) {
				if( fromMotif ) {
					while( i < m && scoresClass1[j] > scoresClass0[i] ) {
						i++;
					}
				} else {
					while( j < d && scoresClass0[i] > scoresClass1[j] ) {
						j++;
					}
				}
			} else {
				while( i + 1 < m && scoresClass0[i] == scoresClass0[i + 1] ) {
					i++;
				}
				while( j + 1 < d && scoresClass1[j] == scoresClass1[j + 1] ) {
					j++;
				}
				i++;
				j++;
			}

			// now we have (i_old,j_old) and (i,j) with i_old <= i and j_old <= j
			if( i == i_old ) {
				//only the number of fp changes
				p[1] = (double)( m - i ) / (double)( m - i + d - j );
				//add point
				if( list != null ) {
					list.add( p.clone() );
				}
			} else {
				if( i < m || j < d ) {
					// the interpolation between point A (i,j) and B (i_old, j_old)

					// FP_A => d-j, FP_B => d-j_old
					// ==> FP_B-FP_A = j - j_old

					// TP_A => m-i, TP_B => m-i_old
					// ==> TP_B - TP_A = i - i_old 

					propTerm = ( j - j_old ) / (double)( i - i_old );
					for( i_old++, helpJ = j_old + propTerm; i_old <= i; i_old++ ) {
						// compute sn and ppv
						help1 = (double)( m - i_old ) / (double)m;
						help2 = (double)( m - i_old ) / (double)( m - i_old + d - helpJ );
						helpJ += propTerm;

						// compute AUC
						erg += ( p[1] + help2 ) / 2d * ( p[0] - help1 );
						p[0] = help1;
						p[1] = help2;

						// add point
						if( list != null ) {
							list.add( p.clone() );
						}
					}
				} else {
					erg += p[1] * p[0];
					p[0] = 0;

					if( list != null ) {
						list.add( p.clone() );
					}
				}
			}

			if( i < m && j < d ) {
				//next
				if( scoresClass0[i] == scoresClass1[j] ) {
					unique = false;
				} else {
					unique = true;
					if( scoresClass0[i] < scoresClass1[j] ) {
						fromMotif = true;
					} else {
						fromMotif = false;
					}
				}
			}
		}

		// left side of the plot
		if( i < m ) {
			help1 = 0;
			//compute AUC
			erg += p[1] * ( p[0] - help1 );
			p[0] = help1;

			// add point
			if( list != null ) {
				list.add( p.clone() );
			}
		}

		return erg;
	}

	/**
	 * This method computes the area under the receiver operator characteristics
	 * curve. It also enables the user to compute the complete curve so that it
	 * can be used for plotting. If you like to get the curve, you should use an
	 * empty {@link AbstractList} for <code>list</code>, otherwise use
	 * <code>null</code>.
	 * 
	 * @param scoresClass0
	 *            the scores of class 0 (have to be sorted)
	 * @param scoresClass1
	 *            the scores of class 1 (have to be sorted)
	 * @param list
	 *            <code>null</code> or an empty list
	 * 
	 * @return the area under the receiver operator characteristics curve
	 */
	public static double getAUC_ROC( double[] scoresClass0, double[] scoresClass1, AbstractList<double[]> list ) {
		int i = 0, j = 0, d = scoresClass1.length, m = scoresClass0.length;
		double erg = 0, help1, help2;
		double[] p = new double[]{ 1, 1 };

		if( list != null ) {
			list.add( p.clone() );
		}
		//which class defines the threshold
		boolean unique, fromMotif = false;
		if( scoresClass0[i] == scoresClass1[j] ) {
			unique = false;
		} else {
			unique = true;
			if( scoresClass0[i] < scoresClass1[j] ) {
				fromMotif = true;
			} else {
				fromMotif = false;
			}
		}
		while( i < m && j < d ) {
			// find next possible threshold
			// discard values that are not interesting
			if( unique ) {
				if( fromMotif ) {
					while( i < m && scoresClass0[i] < scoresClass1[j] ) {
						i++;
					}
				} else {
					while( j < d && scoresClass0[i] > scoresClass1[j] ) {
						j++;
					}
				}
			} else {
				while( i + 1 < m && scoresClass0[i] == scoresClass0[i + 1] ) {
					i++;
				}
				while( j + 1 < d && scoresClass1[j] == scoresClass1[j + 1] ) {
					j++;
				}
				i++;
				j++;
			}

			// erg += height * width
			help1 = (double)( d - j ) / (double)d;
			help2 = (double)( m - i ) / (double)m;
			erg += ( p[1] + help2 ) / 2d * ( p[0] - help1 );
			p[0] = help1;
			p[1] = help2;
			if( list != null ) {
				list.add( p.clone() );
			}

			if( i < m && j < d ) {
				//next
				if( scoresClass0[i] == scoresClass1[j] ) {
					unique = false;
				} else {
					unique = true;
					if( scoresClass0[i] < scoresClass1[j] ) {
						fromMotif = true;
					} else {
						fromMotif = false;
					}
				}
			}
		}
		if( list != null ) {
			list.add( new double[]{ 0, 0 } );
		}
		return erg;
	}

	// double since otherwise an overflow is easily possible
	private static double getCorrelationCoefficient( double tp, double fp, double fn, double tn ) {
		return ( tp * tn - fn * fp ) / Math.sqrt( ( tp + fn ) * ( tn + fp ) * ( tp + fp ) * ( tn + fn ) );
	}

	/**
	 * This method computes the false positive rate (FPR) for a given
	 * sensitivity.
	 * 
	 * @param scoresClass0
	 *            the scores of class 0 (have to be sorted)
	 * @param scoresClass1
	 *            the scores of class 1 (have to be sorted)
	 * @param sensitivity
	 *            the fixed sensitivity which is used to measure the FPR
	 * 
	 * @return the {@link ThresholdMeasurePair} containing the threshold and the
	 *         FPR
	 * 
	 * @throws IllegalArgumentException
	 *             if the sensitivity is not in [0,1]
	 * 
	 * @see ThresholdMeasurePair
	 */
	public static ThresholdMeasurePair getFPRForSensitivity( double[] scoresClass0, double[] scoresClass1, double sensitivity ) throws IllegalArgumentException {
		if( !( 0 <= sensitivity && sensitivity <= 1 ) ) {
			throw new IllegalArgumentException( "The value of percentage has to be in [0,1]." );
		}
		int d = scoresClass1.length, i = d - 1;
		double threshold = scoresClass0[(int)Math.ceil( ( 1 - sensitivity ) * ( scoresClass0.length - 1 ) )];
		while( i >= 0 && scoresClass1[i] >= threshold ) {
			i--;
		}
		return new ThresholdMeasurePair( threshold, (double)( d - 1 - i ) / (double)d );
	}

	/**
	 * This method computes the positive predictive value (PPV) for a given
	 * sensitivity.
	 * 
	 * @param scoresClass0
	 *            the scores of class 0 (have to be sorted)
	 * @param scoresClass1
	 *            the scores of class 1 (have to be sorted)
	 * @param sensitivity
	 *            the fixed sensitivity which is used to measure the PPV
	 * 
	 * @return the {@link ThresholdMeasurePair} containing the threshold and the
	 *         PPV
	 * 
	 * @throws IllegalArgumentException
	 *             if the sensitivity is not in [0,1]
	 * 
	 * @see ThresholdMeasurePair
	 */
	public static ThresholdMeasurePair getPPVForSensitivity( double[] scoresClass0, double[] scoresClass1, double sensitivity ) throws IllegalArgumentException {
		// at least (!!!) <code>sensitivity</code> % of the motifs are correct classified
		if( !( 0 <= sensitivity && sensitivity <= 1 ) ) {
			throw new IllegalArgumentException( "The value of percentage has to be in [0,1]." );
		}
		int d = scoresClass1.length, j = (int)Math.ceil( ( 1 - sensitivity ) * ( scoresClass0.length - 1 ) ), i = d - 1;
		double threshold = scoresClass0[j];
		while( j >= 0 && scoresClass0[j] == threshold ) {
			j--;
		}
		// => (j+1) false negatives
		// => (scoresClass0.length-1-j) true positives
		j = scoresClass0.length - 1 - j; // true positives
		while( i >= 0 && scoresClass1[i] >= threshold ) {
			i--;
		}
		// => (i+1) true negatives
		// => (d-1-i) false positives
		i = d - 1 - i; // false positives
		return new ThresholdMeasurePair( threshold, ( (double)j ) / (double)( i + j ) );
	}

	/**
	 * This method computes the maximal correlation coefficient (CC_max).
	 * 
	 * @param scoresClass0
	 *            the scores of class 0 (have to be sorted)
	 * @param scoresClass1
	 *            the scores of class 1 (have to be sorted)
	 * 
	 * @return the {@link ThresholdMeasurePair} containing the threshold and the
	 *         CC_max
	 * 
	 * @throws Exception
	 *             if the CC_max can not be computed
	 * 
	 * @see ThresholdMeasurePair
	 */
	public static ThresholdMeasurePair getMaxOfCC( double[] scoresClass0, double[] scoresClass1 ) throws Exception {
		int i = 1, j = 1, d = scoresClass1.length, m = scoresClass0.length;
		double[] current = new double[2];
		double[] erg = new double[]{ Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };

		//find second smallest value => first threshold current[0]

		if( scoresClass0[0] == scoresClass1[0] ) {
			while( j < d && scoresClass1[j - 1] == scoresClass1[j] ) {
				j++;
			}
			while( i < m && scoresClass0[i - 1] == scoresClass0[i] ) {
				i++;
			}
			if( i != m && j != d ) {
				current[0] = Math.min( scoresClass1[j], scoresClass0[i] );
				// decoys_lr[j-1] = motifs_lr[i-1] < c <= Math.min( decoys_lr[j], motifs_lr[i] )
			} else if( i == m && j == d ) {
				throw new IllegalArgumentException( "All likelihood-ratios had the same values for both samples." );
			} else if( i == m ) // j != d
			{
				current[0] = scoresClass1[j];
				// motifs_lr[i-1] = decoys_lr[j-1] < c = decoys_lr[j], i == m
			} else
			// if( j == d && i != m )
			{
				current[0] = scoresClass0[i];
				// decoys_lr[j-1] = motifs_lr[i-1] < c = motifs_lr[i] , j == d
			}
		} else {
			if( scoresClass0[0] < scoresClass1[0] ) {
				while( i < m && scoresClass0[i - 1] == scoresClass0[i] ) {
					i++;
				}
				j = 0;
				if( i == m ) {
					current[0] = scoresClass1[j];
					// motifs_lr[i-1] < c == decoys_lr[j], i == m
				} else {
					if( scoresClass0[i] <= scoresClass1[j] ) {
						current[0] = scoresClass0[i];
						// motifs_lr[i-1] < c == motifs_lr[i] <= decoys_lr[j]
					} else {
						current[0] = scoresClass1[j];
						// motifs_lr[i-1] < c == decoys_lr[j] < motifs_lr[i]
					}
				}
			} else // ( motifs_lr[0] > decoys_lr[0] )
			{
				while( j < d && scoresClass1[j - 1] == scoresClass1[j] ) {
					j++;
				}
				i = 0;
				if( j == d ) {
					current[0] = scoresClass0[i];
					//decoys_lr[j-1] < c == motis_lr[i], j == d
				} else {
					if( scoresClass1[j] <= scoresClass0[i] ) {
						current[0] = scoresClass1[j];
						// decoys_lr[j-1] < c = decoys_lr[j] <= motifs_lr[i]
					} else {
						current[0] = scoresClass0[i];
						// motifs_lr[i-1] < c = motifs_lr[i] < decoys_lr[j]
					}
				}
			}
		}
		// forall 0 <= l < i <= k < m:  motifs_lr[l] < c <= motifs_lr[k]
		// forall 0 <= l < j <= k < d:  decoys_lr[l] < c <= decoys_lr[k]

		while( i < m || j < d ) {
			// compute CC
			current[1] = getCorrelationCoefficient( m - i, d - j, i, j );
			if( current[1] > erg[1] ) {
				erg[0] = current[0];
				erg[1] = current[1];
			}

			//find next threshold
			while( j < d && scoresClass1[j] == current[0] ) {
				j++;
			}
			while( i < m && scoresClass0[i] == current[0] ) {
				i++;
			}
			if( i < m && j < d ) {
				if( scoresClass0[i] <= scoresClass1[j] ) {
					current[0] = scoresClass0[i];
				} else // if( motifs_lr[ i ] > decoys_lr[ j ] )
				{
					current[0] = scoresClass1[j];
				}
			} else if( i < m ) // j == d
			{
				current[0] = scoresClass0[i];
			} else if( j < d ) // i == m
			{
				current[0] = scoresClass1[j];
			}
			// else // i == m && j == d => the outer loop will break
		}
		return new ThresholdMeasurePair( erg[0], erg[1] );
	}

	/**
	 * This method computes the sensitivity for a given specificity.
	 * 
	 * @param scoresClass0
	 *            the scores of class 0 (have to be sorted)
	 * @param scoresClass1
	 *            the scores of class 1 (have to be sorted)
	 * @param specificity
	 *            the specificity that is used to measure the sensitivity
	 * 
	 * @return the {@link ThresholdMeasurePair} containing the threshold and the
	 *         sensitivity
	 * 
	 * @throws IllegalArgumentException
	 *             if the specificity is not in [0,1]
	 * 
	 * @see ThresholdMeasurePair
	 */
	public static ThresholdMeasurePair getSensitivityForSpecificity( double[] scoresClass0, double[] scoresClass1, double specificity ) throws IllegalArgumentException {
		if( !( 0 < specificity && specificity < 1 ) ) {
			throw new IllegalArgumentException( "The value of specificity has to be in (0,1)." );
		}
		int i = 0, m = scoresClass0.length;
		double threshold = scoresClass1[(int)Math.ceil( specificity * ( scoresClass1.length - 1 ) )];
		while( i < m && scoresClass0[i] <= threshold ) {
			i++;
		}
		return new ThresholdMeasurePair( threshold, (double)( m - i ) / (double)m );
	}

	/**
	 * This method computes the classification rate for two classes.
	 * 
	 * @param scoresClass0
	 *            the scores of class 0 (have to be sorted)
	 * @param scoresClass1
	 *            the scores of class 1 (have to be sorted)
	 * 
	 * @return the classification rate
	 */
	public static double getClassificationRateFor2Classes( double[] scoresClass0, double[] scoresClass1 ) {
		// TODO change implementation, currently naive
		int i = 0, m = scoresClass0.length;
		while( i < m && scoresClass0[i] < 0 ) {
			i++;
		}

		int d = scoresClass1.length, j = d - 1;
		while( j >= 0 && scoresClass1[j] >= 0 ) {
			j--;
		}

		return ( ( m - i ) + j + 1 ) / (double)( m + d );
	}

	/**
	 * This class is used as a container that allows to store a threshold and
	 * the result of a measure together.
	 * 
	 * @author Jan Grau
	 */
	public static class ThresholdMeasurePair {

		private double threshold;

		private double measure;

		/**
		 * Creates a filled instance of a {@link ThresholdMeasurePair}.
		 * 
		 * @param threshold
		 *            the value of the threshold
		 * @param measure
		 *            the value of the measure
		 */
		public ThresholdMeasurePair( double threshold, double measure ) {
			this.threshold = threshold;
			this.measure = measure;
		}

		/**
		 * This method returns the value of threshold.
		 * 
		 * @return the value of threshold
		 */
		public double getThreshold() {
			return threshold;
		}

		/**
		 * This method returns the value of the measure.
		 * 
		 * @return the value of the measure
		 */
		public double getMeasure() {
			return measure;
		}
		
		public String toString() {
			return "measure: " + measure + "\tthreshold: " + threshold;
		}
	}

}
