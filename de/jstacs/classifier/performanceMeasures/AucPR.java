package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;
import de.jstacs.results.NumericalResultSet;

/**
 * This class implements the area under curve of the precision-recall curve.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class AucPR extends PRCurve implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link AucPR}.
	 */
	public AucPR() {
		super();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link AucPR} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AucPR} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public AucPR(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	public NumericalResultSet compute(double[] sortedScoresClass0, double[] sortedScoresClass1) {
		return (NumericalResultSet) super.compute( sortedScoresClass0, sortedScoresClass1 );
	}

	public NumericalResultSet compute( double[][][] classSpecificScores ) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}
	
	@Override
	public String getName() {
		return "Area under curve (" + NAME + ")";
	}
}
