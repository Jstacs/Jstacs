package de.jstacs.classifiers.performanceMeasures;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;

/**
 * This class implements the false positive rate for a fixed sensitivity.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class FalsePositiveRateForFixedSensitivity extends TwoClassAbstractPerformanceMeasure implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link FalsePositiveRateForFixedSensitivity} with empty parameter values.
	 */
	public FalsePositiveRateForFixedSensitivity() {
		super();
		try{
			parameters.add( new SimpleParameter( DataType.DOUBLE, "Sensitivity", "The fixed sensitivity for the false positive rate.", true, new NumberValidator<Double>(0d,1d),0.95 ) );
		}catch(ParameterException doesnothappen){}
	}
	
	/**
	 * Constructs a new instance of the performance measure {@link FalsePositiveRateForFixedSensitivity} with given <code>sensitivity</code>.
	 * 
	 * @param sensitivity the sensitivity for which the false positive rate should be computed
	 * 
	 * @throws Exception if the internal parameters can not be created or the value can not be set
	 */
	public FalsePositiveRateForFixedSensitivity(double sensitivity) throws Exception {
		this();
		getParameterAt( 0 ).setValue( sensitivity );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link FalsePositiveRateForFixedSensitivity} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link FalsePositiveRateForFixedSensitivity} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
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
		double threshold = sortedScoresClass0[(int)Math.ceil( ( 1 - sensitivity ) * ( sortedScoresClass0.length - 1 ) )];
		while( i >= 0 && sortedScoresClass1[i] >= threshold ) {
			i--;
		}

		return new NumericalResultSet( new NumericalResult[]{
				new NumericalResult( getName() + " of "+sensitivity, "", (double)( d - 1 - i ) / (double)d ),
				new NumericalResult( "Threshold for the " + getName().toLowerCase() + " of "+sensitivity, "", threshold )
		});
	}
	
	public NumericalResultSet compute( double[][][] classSpecificScores ) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}


}
