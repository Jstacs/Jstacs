package de.jstacs.classifier.performanceMeasures;

import java.util.LinkedList;

import de.jstacs.classifier.AbstractClassifier;
import de.jstacs.classifier.assessment.ClassifierAssessment;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;


/**
 * This class implements a container of {@link AbstractPerformanceMeasure}s that can be used
 * in {@link de.jstacs.classifier.AbstractClassifier#evaluate(PerformanceMeasureParameterSet, boolean, de.jstacs.data.DataSet...)}.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class PerformanceMeasureParameterSet extends ExpandableParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link PerformanceMeasureParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link PerformanceMeasureParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public PerformanceMeasureParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	/**
	 * Constructs a new {@link PerformanceMeasureParameterSet} that can be used for binary classifiers.
	 *  
	 * @throws Exception if something went wrong
	 * 
	 * @see #PerformanceMeasureParameterSet(int)
	 */
	public PerformanceMeasureParameterSet() throws Exception {
		this( 2 );
	}
	
	/**
	 * Constructs a new {@link PerformanceMeasureParameterSet} that can be used for classifiers that
	 * handle the given number of classes.
	 * 
	 * @param numClasses the number of classes
	 *  
	 * @throws Exception if something went wrong
	 * 
	 * @see de.jstacs.classifier.AbstractClassifier#getNumberOfClasses()
	 * @see PerformanceMeasureParameterSet#PerformanceMeasureParameterSet(int, AbstractPerformanceMeasure...)
	 */
	public PerformanceMeasureParameterSet( int numClasses ) throws Exception {
		this( numClasses, new AbstractPerformanceMeasure[0] );
	}

	/**
	 * Constructs a new {@link PerformanceMeasureParameterSet} that can be used for classifiers that
	 * handle the given number of classes. The instance contains the given performance measures.
	 * 
	 * @param numClasses the number of classes
	 * @param measures the measures contained in the instance 
	 *  
	 * @throws Exception if something went wrong
	 * 
	 * @see de.jstacs.classifier.AbstractClassifier#getNumberOfClasses()
	 * @see PerformanceMeasureParameterSet#PerformanceMeasureParameterSet(int, SelectionParameter, AbstractPerformanceMeasure... )
	 * @see AbstractPerformanceMeasure#getCollectionOfAllMeasures(int, boolean)
	 */
	public PerformanceMeasureParameterSet( int numClasses, AbstractPerformanceMeasure... measures ) throws Exception {
		this( numClasses, AbstractPerformanceMeasure.getCollectionOfAllMeasures( numClasses, false ), measures );
	}

	private static ParameterSet[] getParameterSets(int numClasses, SelectionParameter selection, AbstractPerformanceMeasure... measures) throws Exception {
		ParameterSet template = new SimpleParameterSet( selection );
		ParameterSet[] pars = new ParameterSet[measures == null || measures.length== 0 ? 1 : measures.length];
		pars[0] = template;
		
		for(int i=0;i<measures.length;i++){
			if(measures[i].getAllowedNumberOfClasses() != 0 && measures[i].getAllowedNumberOfClasses() != numClasses){
				throw new Exception("Provided measure "+measures[i].getName()+" not allowed for "+numClasses+" classes.");
			}
			if( i != 0 ) {
				pars[i] = template.clone();
			}
			((SelectionParameter)pars[i].getParameterAt( 0 )).setValue( measures[i].getName() );
			ParameterSet coll = ((SelectionParameter)pars[i].getParameterAt( 0 )).getParametersInCollection();
			ParameterSetContainer cont = (ParameterSetContainer)coll.getParameterAt( ((SelectionParameter)pars[i].getParameterAt( 0 )).getSelected() );
			cont.setValue( measures[i] );
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
	protected PerformanceMeasureParameterSet( int numClasses, SelectionParameter selection, AbstractPerformanceMeasure... measures ) throws Exception {
		super( getParameterSets(numClasses, selection, measures), "Performance measures", "Performance measures for evaluating a classifier of "+numClasses+" classes" );
	}	
	
	/**
	 * Creates a filled {@link NumericalPerformanceMeasureParameterSet} that can be used in {@link AbstractClassifier#evaluate(PerformanceMeasureParameterSet, boolean, de.jstacs.data.DataSet...)} or in a {@link ClassifierAssessment}.
	 * 
	 * @return a filled {@link NumericalPerformanceMeasureParameterSet}
	 * 
	 * @throws Exception forwarded from {@link #createFilledParameters(boolean, double, double, double, double)} 
	 */
	public static NumericalPerformanceMeasureParameterSet createFilledParameters() throws Exception {
		return (NumericalPerformanceMeasureParameterSet) createFilledParameters( true, 0.999, 0.95, 0.95, 1 ); 
	}
	
	/**
	 * Creates a filled {@link PerformanceMeasureParameterSet} that can be used in {@link AbstractClassifier#evaluate(PerformanceMeasureParameterSet, boolean, de.jstacs.data.DataSet...)}.
	 * If <code>numerical = false</code>, the parameter set also contains curve measures (e.g. {@link PRCurve}, {@link ROCCurve}).
	 * 
	 * @param numerical if <code>true</code the return type is {@link NumericalPerformanceMeasureParameterSet}, otherwise {@link PerformanceMeasureParameterSet} 
	 * @param spForSn the specificity for computing the sensitivity (cf. {@link SensitivityForFixedSpecificity})
	 * @param snForFPR the specificity for computing the sensitivity (cf. {@link FalsePositiveRateForFixedSensitivity})
	 * @param snForPPV the specificity for computing the sensitivity (cf. {@link PositivePredictiveValueForFixedSensitivity})
	 * @param beta the beta of the F-measure (cf. {@link MaximumFMeasure})
	 * 
	 * @return a filled {@link PerformanceMeasureParameterSet}
	 * 
	 * @throws Exception if a performance measure could not be created properly (e.g. wrong parameters: sensitivity < 0, ...) 
	 */
	public static PerformanceMeasureParameterSet createFilledParameters( boolean numerical, double spForSn, double snForFPR, double snForPPV, double beta ) throws Exception {
		PerformanceMeasureParameterSet res;
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
	@Deprecated
	public void addMeasure(AbstractPerformanceMeasure measure) throws CloneNotSupportedException, IllegalValueException{
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
	protected void setMeasure( AbstractPerformanceMeasure measure ) throws IllegalValueException {
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
