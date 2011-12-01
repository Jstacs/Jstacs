package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;
import de.jstacs.results.NumericalResultSet;

public class AucROC extends ROCCurve implements NumericalPerformanceMeasure {

	public AucROC() throws Exception {
		super();
	}

	public AucROC(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
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
