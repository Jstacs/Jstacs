package de.jstacs.classifier.performanceMeasures;

import java.util.LinkedList;

import de.jstacs.NonParsableException;
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

	private static ParameterSet[] getParameterSets(int numClasses, SelectionParameter collection, AbstractPerformanceMeasure... measures) throws Exception {
		ParameterSet template = new SimpleParameterSet( collection );
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
	
	protected PerformanceMeasureParameterSet( int numClasses, SelectionParameter collection, AbstractPerformanceMeasure... measures ) throws Exception {
		super( getParameterSets(numClasses, collection, measures), "Performance measures", "Performance measures for evaluating a classifier of "+numClasses+" classes" );
	}	
	
	public static NumericalPerformanceMeasureParameterSet createFilledParameters() throws Exception {
		return (NumericalPerformanceMeasureParameterSet) createFilledParameters( true, 0.999, 0.95, 0.95, 1 ); 
	}
	
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
	
	@Deprecated//XXX
	public void addMeasure(AbstractPerformanceMeasure measure) throws CloneNotSupportedException, IllegalValueException{
		this.addParameterToSet();
		setMeasure( measure );
	}
	
	protected void setMeasure( AbstractPerformanceMeasure measure ) throws IllegalValueException {
		ParameterSetContainer cont = (ParameterSetContainer)this.parameters.get( this.parameters.size()-1 );
		SelectionParameter sel = (SelectionParameter)cont.getValue().getParameterAt( 0 );
		sel.setValue( measure );
	}

	public AbstractPerformanceMeasure removeMeasure(int index){
		ParameterSetContainer cont = (ParameterSetContainer) this.parameters.remove( index );
		SelectionParameter cp = (SelectionParameter)cont.getValue().getParameterAt( 0 );
		return (AbstractPerformanceMeasure)cp.getValue();
	}
	
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