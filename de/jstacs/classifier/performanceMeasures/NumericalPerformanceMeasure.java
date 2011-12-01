package de.jstacs.classifier.performanceMeasures;

import de.jstacs.results.NumericalResultSet;

public interface NumericalPerformanceMeasure {

	public NumericalResultSet compute(double[] sortedScoresClass0, double[] sortedScoresClass1);

	public NumericalResultSet compute(double[][][] classSpecificScores);
}
