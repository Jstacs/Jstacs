package de.jstacs.classifier.measures;

import de.jstacs.NonParsableException;
import de.jstacs.parameters.CollectionParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.SimpleParameterSet;


public class MeasureParameters extends ExpandableParameterSet {

	public MeasureParameters( int numClasses ) throws Exception {
		super( new SimpleParameterSet( AbstractMeasure.getCollectionOfAllMeasures( numClasses ) ), "Performance measures", "Performance measures for evaluating a classifier of "+numClasses+" classes" );
	}

	public MeasureParameters( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}
	
	public MeasureParameters( int numClasses, AbstractMeasure... measures ) throws Exception {
		super( getTemplate(numClasses, measures), "Performance measures", "Performance measures for evaluating a classifier of "+numClasses+" classes" );
	}

	private static ParameterSet[] getTemplate(int numClasses, AbstractMeasure... measures) throws Exception {
		ParameterSet[] pars = new ParameterSet[measures.length];
		ParameterSet template = new SimpleParameterSet( AbstractMeasure.getCollectionOfAllMeasures( numClasses ) );
		
		for(int i=0;i<pars.length;i++){
			if(measures[i].getAllowedNumberOfClasses() != 0 && measures[i].getAllowedNumberOfClasses() != numClasses){
				throw new Exception("Provided measure "+measures[i].getName()+" not allowed for "+numClasses+" classes.");
			}
			pars[i] = template.clone();
			((CollectionParameter)pars[i].getParameterAt( 0 )).setValue( measures[i].getName() );
			ParameterSet coll = ((CollectionParameter)pars[i].getParameterAt( 0 )).getParametersInCollection();
			ParameterSetContainer cont = (ParameterSetContainer)coll.getParameterAt( ((CollectionParameter)pars[i].getParameterAt( 0 )).getSelected() );
			cont.setValue( measures[i] );
		}
		return pars;
	}
	
	public static MeasureParameters createFilledParameters( boolean numeric, double spForSn, double snForFPR, double snForPPV ) throws Exception {
		MeasureParameters res = new MeasureParameters( 2 );
		res.loadParameters();
		res.setMeasure( new ClassificationRate() );
		res.addMeasure( new SensitivityForFixedSpecificity( spForSn ) );
		res.addMeasure( new FalsePositiveRateForFixedSensitivity( snForFPR ) );
		res.addMeasure( new PositivePredictiveValueForFixedSensitivity( snForPPV ) );
		res.addMeasure( new ROCCurve( !numeric ) );
		res.addMeasure( new PRCurve( !numeric ) );
		res.addMeasure( new MaximumCorrelationCoefficient() );
		return res;
	}
	
	public void addMeasure(AbstractMeasure measure) throws CloneNotSupportedException, IllegalValueException{
		this.addParameterToSet();
		setMeasure( measure );
	}
	
	private void setMeasure( AbstractMeasure measure ) throws IllegalValueException {
		ParameterSetContainer cont = (ParameterSetContainer)this.parameters.get( this.parameters.size()-1 );
		CollectionParameter cp = (CollectionParameter)cont.getValue().getParameterAt( 0 );
		cp.setValue( measure.getName() );
		cont = (ParameterSetContainer)cp.getParametersInCollection().getParameterAt( cp.getSelected() );
		cont.setValue( measure );
	}

	public AbstractMeasure removeMeasure(int index){
		ParameterSetContainer cont = (ParameterSetContainer) this.parameters.remove( index );
		CollectionParameter cp = (CollectionParameter)cont.getValue().getParameterAt( 0 );
		return (AbstractMeasure)cp.getValue();
	}
	
	public int[] removeMeasures(Class<? extends AbstractMeasure> clazz){
		return null;//TODO
	}
	
	public int[] removeMeasures(String name){
		return null;//TODO
	}
	
	public AbstractMeasure[] getAllMeasures(){
		AbstractMeasure[] measures = new AbstractMeasure[parameters.size()];
		for(int i=0;i<measures.length;i++){
			ParameterSetContainer cont = (ParameterSetContainer) this.parameters.get( i );
			CollectionParameter cp = (CollectionParameter)cont.getValue().getParameterAt( 0 );
			measures[i] = (AbstractMeasure)cp.getValue();
		}
		return measures;
	}
	
}
