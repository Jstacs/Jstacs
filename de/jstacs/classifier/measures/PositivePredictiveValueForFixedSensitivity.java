package de.jstacs.classifier.measures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.ResultSet;


public class PositivePredictiveValueForFixedSensitivity extends TwoClassAbstractMeasure {

	public PositivePredictiveValueForFixedSensitivity() {
	}
	
	public PositivePredictiveValueForFixedSensitivity(double sensitivity) throws Exception {
		loadParameters();
		getParameterAt( 0 ).setValue( sensitivity );
	}

	public PositivePredictiveValueForFixedSensitivity( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Positive predictive value for fixed sensitivity";
	}

	@Override
	public ResultSet compute( double[] scoresClass0, double[] scoresClass1 ) {
		double sensitivity = (Double)getParameterAt( 0 ).getValue();
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
		return new NumericalResultSet(new NumericalResult[]{
		                                                    new NumericalResult( "Threshold", "Threshold for the positive predictive value", threshold ),
		                                                    new NumericalResult( "Maximum correlation coefficient", "The positive predictive value for a fixed sensitivity of "+sensitivity, ( (double)j ) / (double)( i + j ) )
		});
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList( 1 );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Sensitivity", "The fixed sensitivity for the positive predictive value.", true ) );
	}

}
