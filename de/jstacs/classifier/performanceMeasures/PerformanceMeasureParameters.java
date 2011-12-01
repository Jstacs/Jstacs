package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;
import de.jstacs.parameters.CollectionParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;


public class PerformanceMeasureParameters extends ExpandableParameterSet {

	public PerformanceMeasureParameters( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}
	
	public PerformanceMeasureParameters( int numClasses ) throws Exception {
		this( numClasses, new AbstractPerformanceMeasure[0] );
	}

	public PerformanceMeasureParameters( int numClasses, AbstractPerformanceMeasure... measures ) throws Exception {
		this( numClasses, AbstractPerformanceMeasure.getCollectionOfAllMeasures( numClasses, false ), measures );
	}

	private static ParameterSet[] getParameterSets(int numClasses, CollectionParameter collection, AbstractPerformanceMeasure... measures) throws Exception {
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
			((CollectionParameter)pars[i].getParameterAt( 0 )).setValue( measures[i].getName() );
			ParameterSet coll = ((CollectionParameter)pars[i].getParameterAt( 0 )).getParametersInCollection();
			ParameterSetContainer cont = (ParameterSetContainer)coll.getParameterAt( ((CollectionParameter)pars[i].getParameterAt( 0 )).getSelected() );
			cont.setValue( measures[i] );
		}
		return pars;
	}
	
	protected PerformanceMeasureParameters( int numClasses, CollectionParameter collection, AbstractPerformanceMeasure... measures ) throws Exception {
		super( getParameterSets(numClasses, collection, measures), "Performance measures", "Performance measures for evaluating a classifier of "+numClasses+" classes" );
	}	
	
	public static PerformanceMeasureParameters createFilledParameters( boolean numerical, double spForSn, double snForFPR, double snForPPV ) throws Exception {
		PerformanceMeasureParameters res;
		if( numerical ) {
			res = new NumericalPerformanceMeasureParameters( 2 );
		} else {
			res = new PerformanceMeasureParameters( 2 );
		}
		res.loadParameters();
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
		return res;
	}
	
	@Deprecated//XXX
	public void addMeasure(AbstractPerformanceMeasure measure) throws CloneNotSupportedException, IllegalValueException{
		this.addParameterToSet();
		setMeasure( measure );
	}
	
	protected void setMeasure( AbstractPerformanceMeasure measure ) throws IllegalValueException {
		ParameterSetContainer cont = (ParameterSetContainer)this.parameters.get( this.parameters.size()-1 );
		CollectionParameter cp = (CollectionParameter)cont.getValue().getParameterAt( 0 );
		cp.setValue( measure.getName() );
		cont = (ParameterSetContainer)cp.getParametersInCollection().getParameterAt( cp.getSelected() );
		cont.setValue( measure );
	}

	public AbstractPerformanceMeasure removeMeasure(int index){
		ParameterSetContainer cont = (ParameterSetContainer) this.parameters.remove( index );
		CollectionParameter cp = (CollectionParameter)cont.getValue().getParameterAt( 0 );
		return (AbstractPerformanceMeasure)cp.getValue();
	}
	
	public int[] removeMeasures(Class<? extends AbstractPerformanceMeasure> clazz){
		return null;//TODO
	}
	
	public int[] removeMeasures(String name){
		return null;//TODO
	}
	
	public AbstractPerformanceMeasure[] getAllMeasures(){
		AbstractPerformanceMeasure[] measures = new AbstractPerformanceMeasure[parameters.size()];
		for(int i=0;i<measures.length;i++){
			ParameterSetContainer cont = (ParameterSetContainer) this.parameters.get( i );
			CollectionParameter cp = (CollectionParameter)cont.getValue().getParameterAt( 0 );
			measures[i] = (AbstractPerformanceMeasure)cp.getValue();
		}
		return measures;
	}
	
}
