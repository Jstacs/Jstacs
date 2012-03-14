package de.jstacs.classifiers.performanceMeasures;

import de.jstacs.io.NonParsableException;
import de.jstacs.results.NumericalResultSet;

/**
 * 
 * @author Jens Keilwagen
 */
public abstract class AbstractNumericalTwoClassPerformanceMeasure extends AbstractTwoClassPerformanceMeasure implements NumericalPerformanceMeasure {

	protected AbstractNumericalTwoClassPerformanceMeasure() {
	}

	protected AbstractNumericalTwoClassPerformanceMeasure(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	public final NumericalResultSet compute( double[][][] classSpecificScores ) {
		return compute( classSpecificScores, null );
	}
	
	public final NumericalResultSet compute( double[][][] classSpecificScores, double[][] weights ) {
		return (NumericalResultSet) super.compute( classSpecificScores, weights );
	}

	@Override
	public final NumericalResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
		return (NumericalResultSet) compute( sortedScoresClass0, null, sortedScoresClass1, null );
	}
	
	public abstract NumericalResultSet compute(double[] sortedScoresClass0, double[] weightClass0, double[] sortedScoresClass1, double[] weightsClass1 );
}
