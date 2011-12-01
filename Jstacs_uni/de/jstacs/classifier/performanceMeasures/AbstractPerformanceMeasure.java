package de.jstacs.classifier.performanceMeasures;

import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.NonParsableException;
import de.jstacs.parameters.CollectionParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.SubclassFinder;


public abstract class AbstractPerformanceMeasure extends ParameterSet {
	
	public AbstractPerformanceMeasure(){
	}
	
	public AbstractPerformanceMeasure(StringBuffer xml) throws NonParsableException{
		super(xml);
	}
	
	public abstract String getName();
	
	public abstract ResultSet compute(double[] classificationScoresFg, double[] classificationScoresBg);

	public abstract ResultSet compute(double[][][] classSpecificScores);
	
	public abstract int getAllowedNumberOfClasses();
	
	public static CollectionParameter getCollectionOfAllMeasures(int numClasses, boolean numerical) throws Exception {
		LinkedList<Class<? extends AbstractPerformanceMeasure>> list = SubclassFinder.findInstantiableSubclasses( AbstractPerformanceMeasure.class, "de.jstacs" );
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
		
		AbstractPerformanceMeasure[] ps = found.toArray( new AbstractPerformanceMeasure[0] );
		String[] keys = new String[ps.length];
		String[] comments = new String[ps.length];
		for(int i=0;i<ps.length;i++){
			keys[i] = ps[i].getName();
			comments[i] = "A performance measure that computes "+keys[i]+".";
		}
		
		CollectionParameter cp = new CollectionParameter( ps, keys, comments, "Performance Measures", "Performance measures that can be computed for "+(numClasses == 0 ? "any number of" : numClasses)+" classes.", true );
		return cp;
	}
	
}
