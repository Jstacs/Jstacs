package de.jstacs.classifier.performanceMeasures;

import de.jstacs.io.NonParsableException;

/**
 * This class implements a container for {@link NumericalPerformanceMeasure}s that can be used, for instance, in an repeated assessment,
 * (cf. {@link de.jstacs.classifier.assessment.ClassifierAssessment}).
 * 
 * @author Jens Keilwagen
 */
public class NumericalPerformanceMeasureParameterSet extends PerformanceMeasureParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link NumericalPerformanceMeasureParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link NumericalPerformanceMeasureParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public NumericalPerformanceMeasureParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	/**
	 * Constructs a new {@link NumericalPerformanceMeasureParameterSet} that can be used for classifiers that
	 * handle the given number of classes.
	 * 
	 * @param numClasses the number of classes
	 *  
	 * @throws Exception if something went wrong
	 * 
	 * @see de.jstacs.classifier.AbstractClassifier#getNumberOfClasses()
	 * @see AbstractPerformanceMeasure#getCollectionOfAllMeasures(int, boolean)
	 * @see PerformanceMeasureParameterSet#PerformanceMeasureParameterSet(int, de.jstacs.parameters.SelectionParameter, AbstractPerformanceMeasure...)
	 */
	public NumericalPerformanceMeasureParameterSet( int numClasses ) throws Exception {
		super( numClasses, AbstractPerformanceMeasure.getCollectionOfAllMeasures( numClasses, true ), new AbstractPerformanceMeasure[0] );
	}
	
	/**
	 * Constructs a new {@link NumericalPerformanceMeasureParameterSet} that can be used for binary classifiers.
	 * 
	 * @throws Exception if something went wrong
	 * 
	 * @see #NumericalPerformanceMeasureParameterSet(int)
	 */
	public NumericalPerformanceMeasureParameterSet() throws Exception {
		this( 2 );
	}
}
