package de.jstacs.classifier.performanceMeasures;

import de.jstacs.io.NonParsableException;

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
	protected String getMeasureName() {
		return "correlation coefficient";
	}
	
	@Override
	protected String getSpecificName() {
		return getName();
	}

	// double since otherwise an overflow is easily possible
	protected double getMeasure( double tp, double fp, double fn, double tn ) {
		return ( tp * tn - fn * fp ) / Math.sqrt( ( tp + fn ) * ( tn + fp ) * ( tp + fp ) * ( tn + fn ) );
	}
}
