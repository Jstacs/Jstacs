package de.jstacs.classifier.performanceMeasures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;


public class FalsePositiveRateForFixedSensitivity extends TwoClassAbstractPerformanceMeasure implements NumericalPerformanceMeasure {

	public FalsePositiveRateForFixedSensitivity() {
	}
	
	public FalsePositiveRateForFixedSensitivity(double sensitivity) throws Exception {
		loadParameters();
		getParameterAt( 0 ).setValue( sensitivity );
	}

	public FalsePositiveRateForFixedSensitivity( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "False positive rate for fixed sensitivity";
	}

	@Override
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
		double sensitivity = (Double)getParameterAt( 0 ).getValue();
		int d = sortedScoresClass1.length, i = d - 1;
		double threshold = sortedScoresClass1[(int)Math.ceil( ( 1 - sensitivity ) * ( sortedScoresClass1.length - 1 ) )];
		while( i >= 0 && sortedScoresClass1[i] >= threshold ) {
			i--;
		}

		return new NumericalResultSet( new NumericalResult[]{
				new NumericalResult( getName(), "The " + getName().toLowerCase() + " of "+sensitivity, (double)( d - 1 - i ) / (double)d ),
				new NumericalResult( "Threshold", "Threshold for the " + getName().toLowerCase() + " of "+sensitivity, threshold )
		});
	}
	
	public NumericalResultSet compute( double[][][] classSpecificScores ) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList( 1 );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Sensitivity", "The fixed sensitivity for the false positive rate.", true, new NumberValidator<Double>(0d,1d) ) );
	}

}
