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
package de.jstacs.classifiers.performanceMeasures;

import java.util.ArrayList;

import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.io.NonParsableException;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.ToolBox;

/**
 * This class implements the precision-recall curve and its area under the curve.
 * The precision-recall curve is the plot of precision (also called positive predictive value, {@latex.inline $\\frac{TP}{TP+FP}$}) 
 * against recall (also called sensitivity, {@latex.inline $\\frac{TP}{TP+FN}$}) for all possible classification thresholds.
 * 
 * If you are only interested in the area under this curve, you can use {@link AucPR} instead.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class PRCurve extends AbstractTwoClassPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link PRCurve}.
	 */
	public PRCurve() {
		super();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link PRCurve} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link PRCurve} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public PRCurve( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	/**
	 * The name of the performance measure return by {@link #getName()} 
	 */
	public static final String NAME = "Precision-Recall curve";
	
	@Override
	public String getName() {
		return NAME;
	}
	
	@Override
	public ResultSet compute( double[] sortedScoresClass0, double[] weightClass0, double[] sortedScoresClass1, double[] weightClass1 ) {
		boolean goadrichAndDavis = simpleWeights( weightClass0 ) && simpleWeights( weightClass1 );
		
		ArrayList<double[]> list = null;
		if( !(this instanceof NumericalPerformanceMeasure) ){
			list = new ArrayList<double[]>();
		}
		
		int k, i_old = 0, j_old = 0, i = 0, j = 0, d = sortedScoresClass1.length, m = sortedScoresClass0.length;
		double helpJ, propTerm;
		double help1 = 0, help2 = 0;
		double aucGD = 0, aucIntegral = 0, pos, neg, fn = 0, tn = 0, tn_old, fn_old;
		if( weightClass0 == null ) {
			pos = m;
		} else {
			pos = ToolBox.sum( weightClass0 );
		}
		if( weightClass1 == null ) {
			neg = d;
		} else {
			neg = ToolBox.sum( weightClass1 );
		}

		// find correct start point
		while( ( j < d ) && ( sortedScoresClass0[i] > sortedScoresClass1[j] ) ) {
			tn += getWeight( weightClass1, j );
			j++;
		}
		//i is zero, so p[0] = 1 ...
		double[] p = new double[]{ ( pos-fn ) / pos, ( pos-fn ) / ( pos-fn + neg-tn ) };
		if( list != null ) {
			list.add( p.clone() );
		}

		// Do positive and negative have different scores 
		boolean unique = !( j < d && sortedScoresClass0[i] == sortedScoresClass1[j] );
		// which class defines the threshold
		boolean fromMotif = unique;
		while( i < m && j < d ) {
			i_old = i;
			j_old = j;
			tn_old = tn;
			fn_old = fn;
			
			// find next possible threshold
			if( !unique || fromMotif ) {
				while( i + 1 < m && sortedScoresClass0[i] == sortedScoresClass0[i + 1] ) {
					fn += getWeight( weightClass0, i );
					i++;
				}
				fn += getWeight( weightClass0, i );
				i++;
			}
			if( !unique || !fromMotif ) {
				while( j + 1 < d && sortedScoresClass1[j] == sortedScoresClass1[j + 1] ) {
					tn += getWeight( weightClass1, j );
					j++;
				}
				tn += getWeight( weightClass1, j );
				j++;
			}
			if( i < m && j < d ) {
				//next
				if( sortedScoresClass0[i] == sortedScoresClass1[j] ) {
					unique = false;
				} else {
					unique = true;
					if( sortedScoresClass0[i] < sortedScoresClass1[j] ) {
						fromMotif = true;
					} else {
						fromMotif = false;
					}
				}
			}

			// now we have (i_old,j_old) and (i,j) with i_old <= i and j_old <= j
			if( fn == fn_old ) {//i == i_old ) {
				//only the number of fp changes
				p[1] = ( pos-fn ) / ( pos-fn + neg-tn );
				//add point
				if( list != null ) {
					list.add( p.clone() );
				}
			} else {
				// the interpolation between point A (i,j) and B (i_old, j_old)
				
				//point A =^= new point
				//point B =^= old point
				
				double pB = p[0];
				double pA = ( pos-fn ) / pos;
				
				//unweighted
				// FP_A => d-j, FP_B => d-j_old
				// ==> FP_B-FP_A = j - j_old

				// TP_A => m-i, TP_B => m-i_old
				// ==> TP_B - TP_A = i - i_old
								
				if( goadrichAndDavis ) {
					if( i < m || j < d ) {
		
						propTerm = ( j - j_old ) / (double)( i - i_old );
						double h1 = p[0], h2 = p[1];
						int c = i_old+1;
						for( helpJ = j_old + propTerm; c <= i; c++ ) {
							// compute sn and ppv
							help1 = (double)( m - c ) / (double)m;
							help2 = (double)( m - c ) / (double)( m - c + d - helpJ );
							helpJ += propTerm;
		
							// compute AUC
							aucGD += ( h2 + help2 ) / 2d * ( h1 - help1 );
							h1 = help1;
							h2 = help2;
						}
					} else {
						//i != i_old, i == m && j == d
						aucGD += p[1] * p[0];
					}
				}
				
				//integral
				
				//weighted
				// FP_A => neg-tn, FP_B => neg-tn_old
				// ==> FP_B-FP_A = tn - tn_old

				// TP_A => pos-fn, TP_B => pos-fn_old
				// ==> TP_B - TP_A = fn - fn_old

				double h = (tn-tn_old) / (fn-fn_old);
				double a = 1 + h;
				double b = (neg-tn - h * (pos-fn) ) / pos;
				
				//System.out.println( a + "\t" + b + "\t" + h);
				if( b != 0 ) {
					aucIntegral += (pB - pA - b/a*(Math.log(a*pB+b) - Math.log(a*pA+b)))/a;
				} else {
					aucIntegral += (pB - pA) / a;
				}

				//curve
				if( list != null ) {
					propTerm = Math.min( (fn - fn_old) / (i - i_old), minStepSize );
					h *= propTerm ;
					double helpI = fn_old + propTerm;
					helpJ = tn_old + h;
					k=1;
					while( helpI < fn ) {
						// compute sn and ppv
						p[0] = ( pos - helpI ) / pos;
						p[1] = ( pos - helpI ) / ( pos - helpI + neg - helpJ );
						
						// add point
						list.add( p.clone() );
						
						k++;
						helpJ = tn_old + k*h;
						helpI = fn_old + k*propTerm;
					}
				}
				if( pA != p[0] ) {
					p[0] = pA;
					double temp = ( pos-fn ) / ( pos-fn + neg-tn );
					if(!Double.isNaN( temp )){
						p[1] = temp;
					}
					if( list != null ) {
						list.add( p.clone() );
					}
				}
			}
			
			/*
			out( sortedScoresClass0, i );
			out( sortedScoresClass1, j );
			out( sortedScoresClass0, i_old );
			out( sortedScoresClass1, j_old );
			System.out.println( aucGD + "\t" + aucIntegral );
			*/
		}

		// left side of the plot
		if( i < m ) {
			help1 = 0;
			if( goadrichAndDavis ) {
				//compute AUC
				aucGD += p[1] * ( p[0] - help1 );
			}
			
			aucIntegral += p[1] * ( p[0] - help1 );
			
			// add point
			p[0] = help1;
			if( list != null ) {
				list.add( p.clone() );
			}
		}
		NumericalResult auc1 = new NumericalResult( "AUC-PR (Davis and Goadrich)", getName(), aucGD ); 
		NumericalResult auc2 = new NumericalResult( "AUC-PR (Integral)", getName(), aucIntegral );
		if(list == null){
			if( goadrichAndDavis ) {
				return new NumericalResultSet( new NumericalResult[]{auc1, auc2} );
			} else {
				return new NumericalResultSet( auc2 );
			}
		}else{
			Result curve = new AbstractScoreBasedClassifier.DoubleTableResult(getName(), getName(), list);
			if( goadrichAndDavis ) {
				return new ResultSet( new Result[]{
	               auc1,
	               auc2,
	               curve
				} );
			} else {
				return new ResultSet( new Result[]{
	               auc2,
	               curve
				} );
			}
		}
	}
		
	public static double minStepSize=1;
	
	private void out( double[] array, int index ) {
		System.out.print( index + " (" + (index < array.length ? array[index] : "###" ) + ")\t" );
	}
}
