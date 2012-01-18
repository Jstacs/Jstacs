package de.jstacs.classifiers.performanceMeasures;

import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.SubclassFinder;

/**
 * This class is the abstract super class of any performance measure used to evaluate
 * an {@link de.jstacs.classifiers.AbstractClassifier}. It is recommended to use the method
 * {@link de.jstacs.classifiers.AbstractClassifier#evaluate(PerformanceMeasureParameterSet, boolean, de.jstacs.data.DataSet...)}
 * for evaluating the performance of any classifier.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see de.jstacs.classifiers.AbstractClassifier
 * @see PerformanceMeasureParameterSet
 */
public abstract class AbstractPerformanceMeasure extends ParameterSet {
	
	/**
	 * Constructs a new {@link AbstractPerformanceMeasure} with empty parameter values.
	 */
	protected AbstractPerformanceMeasure(){}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link AbstractPerformanceMeasure} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AbstractPerformanceMeasure} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	protected AbstractPerformanceMeasure(StringBuffer xml) throws NonParsableException{
		super(xml);
	}
	
	/**
	 * The method returns the name of the performance measure.
	 * 
	 * @return the name of the performance measure
	 */
	public abstract String getName();
	
	/**
	 * This method allows to compute the performance measure of given sorted score ratios.
	 * 
	 * <b>This method can only be used for binary classifiers.</b>
	 * 
	 * @param sortedScoresClass0 the sorted score ratios of class 0
	 * @param sortedScoresClass1 the sorted score ratios of class 1
	 *  
	 * @return a result set containing the results of the performance measure
	 * 
	 * @see java.util.Arrays#sort(double[])
	 */
	public abstract ResultSet compute(double[] sortedScoresClass0, double[] sortedScoresClass1);

	/**
	 * This method allows to compute the performance measure of given class specific scores.
	 * 
	 * @param classSpecificScores the scores; first dimension = data sets, second dimension = sequences of the data set, third dimension classes of the classifier
	 *  
	 * @return a result set containing the results of the performance measure
	 */
	public abstract ResultSet compute(double[][][] classSpecificScores);
	
	/**
	 * This method returns the allowed number of classes. For many performance measures this
	 * number is fixed, e.g. for AUC-ROC the number is 2. If the number is not fixed the
	 * method returns 0, e.g. for the classification rate.
	 * 
	 * @return the allowed number of classes
	 * 
	 * @see de.jstacs.classifiers.AbstractClassifier#getNumberOfClasses()
	 */
	public abstract int getAllowedNumberOfClasses();
	
	/**
	 * This method creates an instance of an {@link SelectionParameter} that can be used to create
	 * an instance of {@link PerformanceMeasureParameterSet} or {@link NumericalPerformanceMeasureParameterSet}.
	 * 
	 * @param numClasses the number of classes
	 * @param numerical 
	 * 		a switch indicating whether all performance measures or only those implementing
	 * 		{@link NumericalPerformanceMeasure} shall be contained in the returned
	 * 		{@link SelectionParameter} 
	 * @return a {@link SelectionParameter} that can be used to create an instance of {@link PerformanceMeasureParameterSet} or {@link NumericalPerformanceMeasureParameterSet}
	 * 
	 * @throws Exception if something went wrong, e.g. missing empty constructor of any performance measure.
	 * 
	 * @see de.jstacs.classifiers.AbstractClassifier#getNumberOfClasses()
	 */
	public static SelectionParameter getCollectionOfAllMeasures(int numClasses, boolean numerical) throws Exception {
		LinkedList<Class<? extends AbstractPerformanceMeasure>> list = SubclassFinder.findInstantiableSubclasses( AbstractPerformanceMeasure.class, AbstractPerformanceMeasure.class.getPackage().getName() );
		Iterator<Class<? extends AbstractPerformanceMeasure>> it = list.iterator();
		LinkedList<AbstractPerformanceMeasure> found = new LinkedList<AbstractPerformanceMeasure>();
		while(it.hasNext()){
			Class<? extends AbstractPerformanceMeasure> cl = it.next();
			if( !numerical || NumericalPerformanceMeasure.class.isAssignableFrom(cl) ) {
				try{
					AbstractPerformanceMeasure mea =  cl.getConstructor().newInstance();
					if(mea.getAllowedNumberOfClasses() == 0 || mea.getAllowedNumberOfClasses() == numClasses){
						found.add( mea );
					}
				}catch(NoSuchMethodException e){
					
				}
			}
		}
		
		return new SelectionParameter( "Performance Measures", "Performance measures that can be computed for "+(numClasses == 0 ? "any number of" : numClasses)+" classes.", true, found.toArray( new AbstractPerformanceMeasure[0] ) );
	}	
}
