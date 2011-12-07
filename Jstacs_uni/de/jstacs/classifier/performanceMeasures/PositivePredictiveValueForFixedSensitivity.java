package de.jstacs.classifier.performanceMeasures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;

/**
 * This class implements the positive predictive value for a fixed sensitivity.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class PositivePredictiveValueForFixedSensitivity extends TwoClassAbstractPerformanceMeasure implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link PositivePredictiveValueForFixedSensitivity} with empty parameter values.
	 */
	public PositivePredictiveValueForFixedSensitivity() {}
	
	/**
	 * Constructs a new instance of the performance measure {@link PositivePredictiveValueForFixedSensitivity} with given <code>sensitivity</code>.
	 * 
	 * @param sensitivity the sensitivity for which the positive predictive value should be computed
	 * 
	 * @throws Exception if the internal parameters can not be created or the value can not be set
	 */
	public PositivePredictiveValueForFixedSensitivity(double sensitivity) throws Exception {
		loadParameters();
		getParameterAt( 0 ).setValue( sensitivity );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link PositivePredictiveValueForFixedSensitivity} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link PositivePredictiveValueForFixedSensitivity} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public PositivePredictiveValueForFixedSensitivity( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Positive predictive value for fixed sensitivity";
	}

	@Override
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
		double sensitivity = (Double)getParameterAt( 0 ).getValue();
		int d = sortedScoresClass1.length, j = (int)Math.ceil( ( 1 - sensitivity ) * ( sortedScoresClass0.length - 1 ) ), i = d - 1;
		double threshold = sortedScoresClass0[j];
		while( j >= 0 && sortedScoresClass0[j] == threshold ) {
			j--;
		}
		// => (j+1) false negatives
		// => (scoresClass0.length-1-j) true positives
		j = sortedScoresClass0.length - 1 - j; // true positives
		while( i >= 0 && sortedScoresClass1[i] >= threshold ) {
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
	
	public NumericalResultSet compute( double[][][] classSpecificScores ) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList( 1 );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Sensitivity", "The fixed sensitivity for the positive predictive value.", true, new NumberValidator<Double>(0d,1d) ) );
	}

}
