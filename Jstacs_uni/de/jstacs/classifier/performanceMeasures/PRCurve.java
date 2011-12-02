package de.jstacs.classifier.performanceMeasures;

import java.util.ArrayList;

import de.jstacs.NonParsableException;
import de.jstacs.classifier.AbstractScoreBasedClassifier;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;


public class PRCurve extends TwoClassAbstractPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link PRCurve}.
	 */
	public PRCurve() throws Exception {
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
	
	public static final String NAME = "Precision-Recall curve";
	
	@Override
	public String getName() {
		return NAME;
	}

	@Override
	public ResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
		
		ArrayList<double[]> list = null;
		if( !(this instanceof NumericalPerformanceMeasure) ){
			list = new ArrayList<double[]>();
		}
		
		int i_old = 0, j_old = 0, i = 0, j = 0, d = sortedScoresClass1.length, m = sortedScoresClass1.length;
		double helpJ, propTerm;
		double erg = 0, help1 = 0, help2 = 0;

		// find correct start point
		while( ( j < d ) && ( sortedScoresClass1[i] > sortedScoresClass1[j] ) ) {
			j++;
		}
		//i is zero, so p[0] = 1 ...
		double[] p = new double[]{ ( m - i ) / (double)m, ( m - i ) / ( (double)( m - i + d - j ) ) };
		if( list != null ) {
			list.add( p.clone() );
		}

		// which class defines the threshold
		boolean unique, fromMotif = false;
		if( j < d && sortedScoresClass1[i] == sortedScoresClass1[j] ) {
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
					while( i < m && sortedScoresClass1[j] > sortedScoresClass1[i] ) {
						i++;
					}
				} else {
					while( j < d && sortedScoresClass1[i] > sortedScoresClass1[j] ) {
						j++;
					}
				}
			} else {
				while( i + 1 < m && sortedScoresClass1[i] == sortedScoresClass1[i + 1] ) {
					i++;
				}
				while( j + 1 < d && sortedScoresClass1[j] == sortedScoresClass1[j + 1] ) {
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
				if( sortedScoresClass1[i] == sortedScoresClass1[j] ) {
					unique = false;
				} else {
					unique = true;
					if( sortedScoresClass1[i] < sortedScoresClass1[j] ) {
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
		NumericalResult auc = new NumericalResult( "AUC-PR", "Area under the " + getName(), erg );
		if(list == null){
			return new NumericalResultSet( auc );
		}else{
			return new ResultSet( new Result[]{
			                                   auc,
			                                   new AbstractScoreBasedClassifier.DoubleTableResult("PR curve", getName(), list)
			} );
		}
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList();
	}
}
