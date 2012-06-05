/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 *
 * For more information on Jstacs, visit http://www.jstacs.de
 */
package de.jstacs.classifiers.performanceMeasures;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.SubclassFinder;
import de.jstacs.utils.ToolBox;

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
	 * @see #compute(double[], double[], double[], double[])
	 */
	public ResultSet compute(double[] sortedScoresClass0, double[] sortedScoresClass1) {
		return compute( sortedScoresClass0, null, sortedScoresClass1, null) ;
	}

	/**
	 * This method allows to compute the performance measure of given class specific scores.
	 * 
	 * @param classSpecificScores the scores; first dimension = data sets, second dimension = sequences of the data set, third dimension classes of the classifier
	 *  
	 * @return a result set containing the results of the performance measure
	 * 
	 * @see #compute(double[][][], double[][])
	 */
	public ResultSet compute(double[][][] classSpecificScores) {
		return compute( classSpecificScores, null );
	}
	
	/**
	 * This method allows to compute the performance measure of given sorted score ratios.
	 * 
	 * <b>This method can only be used for binary classifiers.</b>
	 * 
	 * @param sortedScoresClass0 the sorted score ratios of class 0
	 * @param weightsClass0 the weights of the sequences of class 0 sorted along with the scores <code>sortedScoresClass0</code>
	 * @param sortedScoresClass1 the sorted score ratios of class 1
	 * @param weightsClass1 the weights of the sequences of class 1 sorted along with the scores <code>sortedScoresClass1</code>
	 *  
	 * @return a result set containing the results of the performance measure
	 * 
	 * @see de.jstacs.utils.ToolBox#sortAlongWith(double[], double[])
	 */
	public abstract ResultSet compute(double[] sortedScoresClass0, double[] weightsClass0, double[] sortedScoresClass1, double[] weightsClass1 );

	/**
	 * This method allows to compute the performance measure of given class specific scores.
	 * 
	 * @param classSpecificScores the scores; first dimension = data sets, second dimension = sequences of the data set, third dimension classes of the classifier
	 * @param weights the weights for all sequence in all data sets
	 *  
	 * @return a result set containing the results of the performance measure
	 */
	public abstract ResultSet compute(double[][][] classSpecificScores, double[][] weights );
	
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
	
	protected static double findThreshold(double[] sortedReferenceScores, double[] sortedMeasureScores, double[] referenceWeights, double percentage, boolean atLeast){
		
		double sum = referenceWeights == null ? sortedReferenceScores.length : ToolBox.sum( 0,sortedReferenceScores.length, referenceWeights );
		double hypThreshWeight = sum*percentage;
		
		double curr = 0;
		int i=0;
		while( ( curr ) < hypThreshWeight && i < sortedReferenceScores.length ){
			curr += getWeight( referenceWeights, i );
			i++;
		}

		//skip all items with the same score
		if( atLeast ) {
			while( i+1 < sortedReferenceScores.length && sortedReferenceScores[i] == sortedReferenceScores[i+1] ) {
				i++;
			}
		} else {
			while( i-1 >= 0 && sortedReferenceScores[i] == sortedReferenceScores[i-1] ) {
				curr -= getWeight( referenceWeights, i );
				i--;
			}
			if( curr != hypThreshWeight ){
				// We did not exactly meet the percentage and want that percentage at most.
				// Since we are a bit above the percentage now, we decrease i.
				i--;
			}
		}
		
		double min;
		if( i == sortedReferenceScores.length ){
			//we have all scores below threshold, 
			//so we need at least the last score (plus epsilon)
			min = sortedReferenceScores[ sortedReferenceScores.length-1 ];
		}else{
			//we use the score that does almost meet the weight (!atLeast) or is a bit above (atLeast)
			min = sortedReferenceScores[i];
		}
		if(!atLeast){
			//if we want at most that percentage, we can return what we have to far
			return min;
		}else{
			//otherwise we need to find the optimal threshold above min
			
			//we search for the score above min in sortedMeasureScores
			int j = findSplitIndex( sortedMeasureScores, min );
			while( j < sortedMeasureScores.length && sortedMeasureScores[j] == min ){
				j++;
			}
			double max;
			if(j == sortedMeasureScores.length){
				//if all scores in sortedMeasureScores are already below min
				// we can savely return a threshold above this value
				max =  min * 1.01;//TODO what else?
			}else{
				//otherwise the threshold should be below the next-greater value,
				//which we store to max
				max = sortedMeasureScores[j];
			}
			//now we need to find the first value above min in sortedReferenceScores
			//since this might be (substantially) lower than max
			while(i < sortedReferenceScores.length && sortedReferenceScores[i] == min){
				i++;
			}
			//if we find such a value
			if(i < sortedReferenceScores.length){
				//we set the new max to the smaller of that value
				//and the previous max
				max = Math.min( max, sortedReferenceScores[i] );
			}
			//we return the middle value between
			// - the value that is an epsilon below the required threshold and
			// - the value that is the smallest value greater than min
			//   in any of the arrays sortedReferenceScore and sortedMeasureScores
			return (min+max)/2.0;
		}
	}
	
	/**
	 * Returns the index in <code>sortedScores</code> with value greater or equal to <code>t</code>.
	 * 
	 * @param sortedScores the array of sorted values
	 * @param t the threshold
	 * @return the index in <code>sortedScores</code> with value greater or equal to <code>t</code>
	 */
	protected static int findSplitIndex( double[] sortedScores, double t ) {
		int i = Arrays.binarySearch( sortedScores, t );
		if( i >= 0 ) {
			//find first occurrence
			while( i >= 0 && sortedScores[i] == t ) {
				i--;
			}
			i++;
		} else {
			//compute insertion point
			//binary search returns i = (-(insertion point) - 1);
			i = -(i+1);
		}
		return i;
	}
	
	protected final static double getWeight( double[] weight, int index ) {
		if( weight == null ) {
			return 1;
		} else {
			return weight[index];
		}
	}
	
	public static boolean simpleWeights( double[] weight ) {
		if( weight != null ) {
			for( int i = 0; i < weight.length; i++ ) {
				if( weight[i] != 1 ) {
					return false;
				}
			}
		}
		return true;
	}
}
