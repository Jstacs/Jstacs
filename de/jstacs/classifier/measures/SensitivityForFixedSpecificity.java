package de.jstacs.classifier.measures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.classifier.ScoreBasedPerformanceMeasureDefinitions.ThresholdMeasurePair;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.ResultSet;


public class SensitivityForFixedSpecificity extends TwoClassAbstractMeasure {

	public SensitivityForFixedSpecificity() {
	}
	
	public SensitivityForFixedSpecificity(double specificity) throws Exception {
		loadParameters();
		getParameterAt( 0 ).setValue( specificity );
	}

	public SensitivityForFixedSpecificity( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Sensitivity for a fixed specificity";
	}

	@Override
	public ResultSet compute( double[] scoresClass0, double[] scoresClass1 ) {
		double specificity = (Double)getParameterAt( 0 ).getValue();
		if( !( 0 < specificity && specificity < 1 ) ) {
			throw new IllegalArgumentException( "The value of specificity has to be in (0,1)." );
		}
		int i = 0, m = scoresClass0.length;
		double threshold = scoresClass1[(int)Math.ceil( specificity * ( scoresClass1.length - 1 ) )];
		while( i < m && scoresClass0[i] <= threshold ) {
			i++;
		}
		
		return new NumericalResultSet(new NumericalResult[]{
		                                                    new NumericalResult( "Threshold", "Threshold for the sensitivity", threshold ),
		                                                    new NumericalResult( "Sensitivity", "The sensitivity for a specificity of "+specificity, (double)( m - i ) / (double)m  )
		});
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList( 1 );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Specificity", "The fixed specificity for the sensitivity.", true ) );
	}

}
