package de.jstacs.classifier.measures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.ResultSet;


public class FalsePositiveRateForFixedSensitivity extends TwoClassAbstractMeasure {

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
	public ResultSet compute( double[] scoresClass0, double[] scoresClass1 ) {
		double sensitivity = (Double)getParameterAt( 0 ).getValue();
		if( !( 0 <= sensitivity && sensitivity <= 1 ) ) {
			throw new IllegalArgumentException( "The value of percentage has to be in [0,1]." );
		}
		int d = scoresClass1.length, i = d - 1;
		double threshold = scoresClass0[(int)Math.ceil( ( 1 - sensitivity ) * ( scoresClass0.length - 1 ) )];
		while( i >= 0 && scoresClass1[i] >= threshold ) {
			i--;
		}

		return new NumericalResultSet(new NumericalResult[]{
		                                                    new NumericalResult( "Threshold", "Threshold for the false positive rate", threshold ),
		                                                    new NumericalResult( "Maximum correlation coefficient", "The false positive rate for a fixed sensitivity of "+sensitivity, (double)( d - 1 - i ) / (double)d )
		});
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList( 1 );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Sensitivity", "The fixed sensitivity for the false positive rate.", true ) );
	}

}
