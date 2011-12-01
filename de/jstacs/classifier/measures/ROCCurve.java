package de.jstacs.classifier.measures;

import java.util.ArrayList;

import de.jstacs.NonParsableException;
import de.jstacs.classifier.AbstractScoreBasedClassifier;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;


public class ROCCurve extends CurveMeasure {

	public ROCCurve() throws Exception {
		this(false);
	}

	public ROCCurve( boolean getCurve ) throws Exception {
		super( getCurve );
	}

	public ROCCurve( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Receiver Operating Characteristic";
	}

	@Override
	public ResultSet compute( double[] scoresClass0, double[] scoresClass1 ) {
		
		ArrayList<double[]> list = null;
		if((Boolean)getParameterAt( 0 ).getValue()){
			list = new ArrayList<double[]>();
		}
		
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
		if(list == null){
			return new NumericalResultSet( new NumericalResult( "AUC-ROC", "Area under the ROC curve", erg ) );
		}else{
			return new ResultSet( new Result[]{
			                                   new NumericalResult( "AUC-ROC", "Area under the ROC curve", erg ),
			                                   new AbstractScoreBasedClassifier.DoubleTableResult("ROC curve", "Receiver operating characteristic curve", list)
			} );
		}
	}

}
