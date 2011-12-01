package de.jstacs.classifier.measures;

import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.NonParsableException;
import de.jstacs.parameters.CollectionParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.SubclassFinder;


public abstract class AbstractMeasure extends ParameterSet {
	
	public AbstractMeasure(){
	}
	
	public AbstractMeasure(StringBuffer xml) throws NonParsableException{
		super(xml);
	}
	
	public abstract String getName();
	
	public abstract ResultSet compute(double[] classificationScoresFg, double[] classificationScoresBg);

	public abstract ResultSet compute(double[][][] classSpecificScores);
	
	public abstract int getAllowedNumberOfClasses();
	
	public static CollectionParameter getCollectionOfAllMeasures(int numClasses) throws Exception {
		LinkedList<Class<? extends AbstractMeasure>> list = SubclassFinder.findInstantiableSubclasses( AbstractMeasure.class, "de.jstacs" );
		Iterator<Class<? extends AbstractMeasure>> it = list.iterator();
		LinkedList<AbstractMeasure> found = new LinkedList<AbstractMeasure>();
		while(it.hasNext()){
			Class<? extends AbstractMeasure> cl = it.next();
			try{
				AbstractMeasure mea =  cl.getConstructor().newInstance();
				if(mea.getAllowedNumberOfClasses() == 0 || mea.getAllowedNumberOfClasses() == numClasses){
					found.add( mea );
				}
			}catch(NoSuchMethodException e){
				
			}
		}
		
		AbstractMeasure[] ps = found.toArray( new AbstractMeasure[0] );
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
