package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;

public class NumericalPerformanceMeasureParameters extends PerformanceMeasureParameters {

	public NumericalPerformanceMeasureParameters( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}
	
	public NumericalPerformanceMeasureParameters( int numClasses ) throws Exception {
		super( numClasses, AbstractPerformanceMeasure.getCollectionOfAllMeasures( numClasses, true ), new AbstractPerformanceMeasure[0] );
	}
}
