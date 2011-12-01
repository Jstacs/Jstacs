package de.jstacs.classifier.measures;

import java.util.ArrayList;

import de.jstacs.NonParsableException;
import de.jstacs.classifier.AbstractScoreBasedClassifier;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;


public class PRCurve extends CurveMeasure {

	public PRCurve() throws Exception {
		this(false);
	}

	public PRCurve( boolean getCurve ) throws Exception {
		super( getCurve );
	}

	public PRCurve( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Precision-Recall curve";
	}

	@Override
	public ResultSet compute( double[] scoresClass0, double[] scoresClass1 ) {
		
		ArrayList<double[]> list = null;
		if((Boolean)getParameterAt( 0 ).getValue()){
			list = new ArrayList<double[]>();
		}
		
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
		if(list == null){
			return new NumericalResultSet( new NumericalResult( "AUC-PR", "Area under the PR curve", erg ) );
		}else{
			return new ResultSet( new Result[]{
			                                   new NumericalResult( "AUC-PR", "Area under the PR curve", erg ),
			                                   new AbstractScoreBasedClassifier.DoubleTableResult("PR curve", "Precision recall curve", list)
			} );
		}
	}

}
