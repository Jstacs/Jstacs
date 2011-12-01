package de.jstacs.classifier.performanceMeasures;

import de.jstacs.results.NumericalResultSet;

public interface NumericalPerformanceMeasure {

	public NumericalResultSet compute(double[] classificationScoresFg, double[] classificationScoresBg);

	public NumericalResultSet compute(double[][][] classSpecificScores);
}
