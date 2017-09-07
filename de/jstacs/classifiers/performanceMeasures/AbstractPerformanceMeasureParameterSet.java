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

import java.util.LinkedList;

import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.SimpleParameterSet;


/**
 * This class implements a container of {@link PerformanceMeasure}s that can be used
 * in {@link de.jstacs.classifiers.AbstractClassifier#evaluate(AbstractPerformanceMeasureParameterSet, boolean, de.jstacs.data.DataSet...)}.
 * 
 * @author Jens Keilwagen, Jan Grau
 * @param <T> the class of the {@link PerformanceMeasure}, typically {@link NumericalPerformanceMeasure} or {@link PerformanceMeasure}
 */
public class AbstractPerformanceMeasureParameterSet<T extends PerformanceMeasure> extends ExpandableParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link AbstractPerformanceMeasureParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AbstractPerformanceMeasureParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public AbstractPerformanceMeasureParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	/**
	 * Returns the number of classes the {@link PerformanceMeasure}s in <code>measures</code>
	 * can be applied to.
	 * @param measures the measures
	 * @return the number of classes
	 */
	protected static int getNumberOfClasses( PerformanceMeasure[] measures ) {
		int res = 0;
		for( int i = 0; i < measures.length; i++ ) {
			int n = measures[i].getAllowedNumberOfClasses();
			if( res == 0 ) {
				if( n != 0 ) {
					res = n;
				}
			} else {
				if( n != res ) {
					throw new IllegalArgumentException( "The performance measures are defined for different number of classes." );
				}
			}
		}
		return res==0?2:res;
	}

	/**
	 * Creates a new {@link AbstractPerformanceMeasureParameterSet} for the given number of classes and measures using only numerical
	 * performance measures or not.
	 * 
	 * @param numClasses the number of classes
	 * @param numerical the information if all measures shall be {@link NumericalPerformanceMeasure}s
	 * @param measures the measures
	 * @throws Exception if something went wrong
	 */
	protected AbstractPerformanceMeasureParameterSet( int numClasses, boolean numerical, T[] measures ) throws Exception {
		this( numClasses, AbstractPerformanceMeasure.getCollectionOfAllMeasures( numClasses, numerical ), measures );
	}

	private static ParameterSet[] getParameterSets(int numClasses, SelectionParameter selection, PerformanceMeasure... measures) throws Exception {
		ParameterSet template = new SimpleParameterSet( selection );
		int n = measures == null || measures.length== 0 ? 0 : measures.length;
		ParameterSet[] pars = new ParameterSet[Math.max(1,n)];
		pars[0] = template;
		
		for(int i=0;i<n;i++){
			if(measures[i].getAllowedNumberOfClasses() != 0 && measures[i].getAllowedNumberOfClasses() != numClasses){
				throw new Exception("Provided measure "+measures[i].getName()+" not allowed for "+numClasses+" classes.");
			}
			if( i != 0 ) {
				pars[i] = template.clone();
			}
			((SelectionParameter)pars[i].getParameterAt( 0 )).setValue( measures[i] );
		}
		return pars;
	}
	
	/**
	 * This constructor creates an instance with a given template <code>selection</code> that can be used for classifiers handling a given number of classes.
	 * Additional it allows to set some measure initially.
	 * 
	 * @param numClasses the number of classes
	 * @param selection the template that can be used to add an select performance measures
	 * @param measures the initially set measures
	 * 
	 * @throws Exception if the <code>measures</code> could not be set (e.g. number of classes differs, ...) 
	 */
	protected AbstractPerformanceMeasureParameterSet( int numClasses, SelectionParameter selection, T... measures ) throws Exception {
		super( getParameterSets(numClasses, selection, measures), "Performance measures", "Performance measures for evaluating a classifier of "+numClasses+" classes" );
	}	
	
	/**
	 * Creates a filled {@link NumericalPerformanceMeasureParameterSet} that can be used in
	 * {@link de.jstacs.classifiers.AbstractClassifier#evaluate(AbstractPerformanceMeasureParameterSet, boolean, de.jstacs.data.DataSet...)}
	 * or in a {@link de.jstacs.classifiers.assessment.ClassifierAssessment}.
	 * 
	 * @return a filled {@link NumericalPerformanceMeasureParameterSet}
	 * 
	 * @throws Exception forwarded from {@link #createFilledParameters(boolean, double, double, double, double)} 
	 */
	public static NumericalPerformanceMeasureParameterSet createFilledParameters() throws Exception {
		return (NumericalPerformanceMeasureParameterSet) createFilledParameters( true, 0.999, 0.95, 0.95, 1 ); 
	}
	
	/**
	 * Creates a filled {@link AbstractPerformanceMeasureParameterSet} that can be used in
	 * {@link de.jstacs.classifiers.AbstractClassifier#evaluate(AbstractPerformanceMeasureParameterSet, boolean, de.jstacs.data.DataSet...)}.
	 * If <code>numerical = false</code>, the parameter set also contains curve measures (e.g. {@link PRCurve}, {@link ROCCurve}).
	 * 
	 * @param numerical if <code>true</code> the return type is {@link NumericalPerformanceMeasureParameterSet}, otherwise {@link AbstractPerformanceMeasureParameterSet} 
	 * @param spForSn the specificity for computing the sensitivity (cf. {@link SensitivityForFixedSpecificity})
	 * @param snForFPR the specificity for computing the sensitivity (cf. {@link FalsePositiveRateForFixedSensitivity})
	 * @param snForPPV the specificity for computing the sensitivity (cf. {@link PositivePredictiveValueForFixedSensitivity})
	 * @param beta the beta of the F-measure (cf. {@link MaximumFMeasure})
	 * 
	 * @return a filled {@link AbstractPerformanceMeasureParameterSet}
	 * 
	 * @throws Exception if a performance measure could not be created properly (e.g. wrong parameters: sensitivity &lt; 0, ...) 
	 */
	public static AbstractPerformanceMeasureParameterSet createFilledParameters( boolean numerical, double spForSn, double snForFPR, double snForPPV, double beta ) throws Exception {
		AbstractPerformanceMeasureParameterSet res;
		if( numerical ) {
			res = new NumericalPerformanceMeasureParameterSet( 2 );
		} else {
			res = new PerformanceMeasureParameterSet( 2 );
		}
		res.setMeasure( new ClassificationRate() );
		res.addMeasure( new SensitivityForFixedSpecificity( spForSn ) );
		res.addMeasure( new FalsePositiveRateForFixedSensitivity( snForFPR ) );
		res.addMeasure( new PositivePredictiveValueForFixedSensitivity( snForPPV ) );
		if( numerical ) {
			res.addMeasure( new AucROC() );
			res.addMeasure( new AucPR() );
		} else {
			res.addMeasure( new ROCCurve() );
			res.addMeasure( new PRCurve() );
		}
		res.addMeasure( new MaximumCorrelationCoefficient() );
		res.addMeasure( new MaximumFMeasure(beta) );
		return res;
	}
	
	/**
	 * Adds a performance measure to the set, i.e., enlarges the set by one and set the given measure.
	 * 
	 * @param measure the measure to be added
	 * 
	 * @throws CloneNotSupportedException if the template could not be cloned (forwarded from {@link ExpandableParameterSet#addParameterToSet()}
	 * @throws IllegalValueException forwarded from {@link #setMeasure(AbstractPerformanceMeasure)}
	 */
	//@Deprecated
	private void addMeasure( T measure) throws CloneNotSupportedException, IllegalValueException{
		this.addParameterToSet();
		setMeasure( measure );
	}
	
	/**
	 * Sets the given measure as content of the internally last {@link ParameterSetContainer}.
	 * 
	 * @param measure the measure to be computed
	 * 
	 * @throws IllegalValueException if the measure could not be set (forwarded from {@link SelectionParameter#setValue(Object)})
	 */
	protected void setMeasure( T measure ) throws IllegalValueException {
		ParameterSetContainer cont = (ParameterSetContainer)this.parameters.get( this.parameters.size()-1 );
		SelectionParameter sel = (SelectionParameter)cont.getValue().getParameterAt( 0 );
		sel.setValue( measure );
	}

	/**
	 * Removes the measure with index <code>index</code> from the set and returns this measure.
	 * 
	 * @param index the index of the measure to be removed
	 * 
	 * @return the removed measure.
	 */
	public AbstractPerformanceMeasure removeMeasure(int index){
		ParameterSetContainer cont = (ParameterSetContainer) this.parameters.remove( index );
		SelectionParameter cp = (SelectionParameter)cont.getValue().getParameterAt( 0 );
		return (AbstractPerformanceMeasure)cp.getValue();
	}
	
	/**
	 * Removes all measures of a specific class from the set.
	 * 
	 * @param clazz the specific class 
	 * 
	 * @return the removed measures as an array
	 * 
	 * @see Class#equals(Object)
	 */
	public AbstractPerformanceMeasure[] removeMeasures(Class<? extends AbstractPerformanceMeasure> clazz){
		LinkedList<AbstractPerformanceMeasure> list = new LinkedList<AbstractPerformanceMeasure>();
		for(int i=0;i<parameters.size();i++){
			if(((ParameterSet)this.parameters.get( i ).getValue()).getParameterAt( 0 ).getValue().getClass().equals( clazz )){
				list.add( removeMeasure( i ) );
				i--;
			}
		}
		return list.toArray( new AbstractPerformanceMeasure[0] );
	}
	
	/**
	 * Removes all measures with a specific name.
	 * 
	 * @param name the specific name 
	 * 
	 * @return the removed measures as an array
	 * 
	 * @see AbstractPerformanceMeasure#getName()
	 * @see String#equals(Object)
	 */
	public AbstractPerformanceMeasure[] removeMeasures(String name){
		LinkedList<AbstractPerformanceMeasure> list = new LinkedList<AbstractPerformanceMeasure>();
		for(int i=0;i<parameters.size();i++){
			AbstractPerformanceMeasure meas = (AbstractPerformanceMeasure)((ParameterSet)this.parameters.get( i ).getValue()).getParameterAt( 0 ).getValue();
			if(meas.getName().equals( name )){
				list.add( removeMeasure( i ) );
				i--;
			}
		}
		return list.toArray( new AbstractPerformanceMeasure[0] );
	}
	
	/**
	 * Returns an array of all contained performance measures.
	 * 
	 * @return an array of all contained performance measures.
	 */
	public AbstractPerformanceMeasure[] getAllMeasures(){
		AbstractPerformanceMeasure[] measures = new AbstractPerformanceMeasure[parameters.size()];
		for(int i=0;i<measures.length;i++){
			ParameterSetContainer cont = (ParameterSetContainer) this.parameters.get( i );
			SelectionParameter cp = (SelectionParameter)cont.getValue().getParameterAt( 0 );
			measures[i] = (AbstractPerformanceMeasure)cp.getValue();
		}
		return measures;
	}	
}
