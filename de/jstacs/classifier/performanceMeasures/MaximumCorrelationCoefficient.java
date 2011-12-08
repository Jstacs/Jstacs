package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;

/**
 * This class implements the maximum of the correlation coefficient.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class MaximumCorrelationCoefficient extends MaximumNumericalTwoClassMeasure implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link MaximumCorrelationCoefficient}.
	 */
	public MaximumCorrelationCoefficient() {}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link MaximumCorrelationCoefficient} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MaximumCorrelationCoefficient} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public MaximumCorrelationCoefficient( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "correlation coefficient";
	}

	// double since otherwise an overflow is easily possible
	protected double getMeasure( double tp, double fp, double fn, double tn ) {
		return ( tp * tn - fn * fp ) / Math.sqrt( ( tp + fn ) * ( tn + fp ) * ( tp + fp ) * ( tn + fn ) );
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList();
	}
}
