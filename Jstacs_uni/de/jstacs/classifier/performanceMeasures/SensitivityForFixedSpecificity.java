package de.jstacs.classifier.performanceMeasures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;


public class SensitivityForFixedSpecificity extends TwoClassAbstractPerformanceMeasure implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link SensitivityForFixedSpecificity} with empty parameter values.
	 */
	public SensitivityForFixedSpecificity() {
	}
	
	public SensitivityForFixedSpecificity(double specificity) throws Exception {
		loadParameters();
		getParameterAt( 0 ).setValue( specificity );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link SensitivityForFixedSpecificity} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SensitivityForFixedSpecificity} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public SensitivityForFixedSpecificity( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Sensitivity for a fixed specificity";
	}

	@Override
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
		double specificity = (Double)getParameterAt( 0 ).getValue();
		int i = 0, m = sortedScoresClass1.length;
		double threshold = sortedScoresClass1[(int)Math.ceil( specificity * ( sortedScoresClass1.length - 1 ) )];
		while( i < m && sortedScoresClass1[i] <= threshold ) {
			i++;
		}
		
		return new NumericalResultSet(new NumericalResult[]{
				new NumericalResult( "Sensitivity", "The "+getName().toLowerCase() +" of "+specificity, (double)( m - i ) / (double)m  ),
				new NumericalResult( "Threshold", "Threshold for the "+getName().toLowerCase() +" of "+specificity, threshold )
		});
	}

	public NumericalResultSet compute( double[][][] classSpecificScores ) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}
	
	@Override
	protected void loadParameters() throws Exception {
		initParameterList( 1 );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Specificity", "The fixed specificity for the sensitivity.", true, new NumberValidator<Double>(0d,1d) ) );
	}
}
