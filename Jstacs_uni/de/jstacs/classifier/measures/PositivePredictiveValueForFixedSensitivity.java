package de.jstacs.classifier.measures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;


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
	public NumericalResultSet compute( double[] scoresClass0, double[] scoresClass1 ) {
		double sensitivity = (Double)getParameterAt( 0 ).getValue();
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
				new NumericalResult( getName(), "The " + getName().toLowerCase() + " of "+sensitivity, ( (double)j ) / (double)( i + j ) ),
				new NumericalResult( "Threshold", "Threshold for the " + getName().toLowerCase() + " of "+sensitivity, threshold ),
		});
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList( 1 );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Sensitivity", "The fixed sensitivity for the positive predictive value.", true, new NumberValidator<Double>(0d,1d) ) );
	}

}
