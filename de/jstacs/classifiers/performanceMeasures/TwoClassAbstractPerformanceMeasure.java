package de.jstacs.classifiers.performanceMeasures;

import java.util.Arrays;

import de.jstacs.io.NonParsableException;
import de.jstacs.results.ResultSet;

/**
 * This class is the abstract super class of any performance measure that can only be computed for binary classifiers.
 *  
 * @author Jan Grau, Jens Keilwagen
 */
public abstract class TwoClassAbstractPerformanceMeasure extends AbstractPerformanceMeasure {

	/**
	 * Constructs a new {@link TwoClassAbstractPerformanceMeasure} with empty parameter values.
	 */
	protected TwoClassAbstractPerformanceMeasure() {}	

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link TwoClassAbstractPerformanceMeasure} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link TwoClassAbstractPerformanceMeasure} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	protected TwoClassAbstractPerformanceMeasure( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public ResultSet compute( double[][][] classSpecificScores ) {
		if(classSpecificScores.length != 2){
			throw new RuntimeException( "Only two classes possible for "+ getName() );
		}
		double[][] classificationScores = new double[2][];
		for(int i=0;i<classSpecificScores.length;i++){
			classificationScores[i] = new double[classSpecificScores[i].length];
			for(int j=0;j<classSpecificScores[i].length;j++){
				if(classSpecificScores[i][j].length != 2){
					throw new RuntimeException( "Only two classes possible for "+ getName() );
				}
				classificationScores[i][j] = classSpecificScores[i][j][0] - classSpecificScores[i][j][1];
			}
			Arrays.sort( classificationScores[i] );
		}
		return compute( classificationScores[0], classificationScores[1] );
	}

	@Override
	public final int getAllowedNumberOfClasses() {
		return 2;
	}
}