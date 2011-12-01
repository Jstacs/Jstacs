package de.jstacs.classifier.performanceMeasures;

import java.util.ArrayList;

import de.jstacs.NonParsableException;
import de.jstacs.classifier.AbstractScoreBasedClassifier;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;


public class ROCCurve extends TwoClassAbstractPerformanceMeasure {

	public ROCCurve() throws Exception {
		super();
	}

	public ROCCurve( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	public static final String NAME = "Receiver Operating Characteristic";
	
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
		
		int i = 0, j = 0, d = sortedScoresClass1.length, m = sortedScoresClass1.length;
		double erg = 0, help1, help2;
		double[] p = new double[]{ 1, 1 };

		if( list != null ) {
			list.add( p.clone() );
		}
		//which class defines the threshold
		boolean unique, fromMotif = false;
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
		while( i < m && j < d ) {
			// find next possible threshold
			// discard values that are not interesting
			if( unique ) {
				if( fromMotif ) {
					while( i < m && sortedScoresClass1[i] < sortedScoresClass1[j] ) {
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
		if( list != null ) {
			list.add( new double[]{ 0, 0 } );
		}
		NumericalResult auc = new NumericalResult( "AUC-ROC", "Area under the " + getName() + " curve", erg );
		if(list == null){
			return new NumericalResultSet( auc );
		}else{
			return new ResultSet( new Result[]{
			                                   auc,
			                                   new AbstractScoreBasedClassifier.DoubleTableResult("ROC curve", getName() + " curve", list)
			} );
		}
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList();
	}
}
