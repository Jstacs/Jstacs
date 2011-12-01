package de.jstacs.classifier.measures;

import de.jstacs.NonParsableException;
import de.jstacs.results.ResultSet;


public abstract class TwoClassAbstractMeasure extends AbstractMeasure {

	public TwoClassAbstractMeasure() {
		super();
	}	

	/**
	 * @param xml
	 * @throws NonParsableException
	 */
	public TwoClassAbstractMeasure( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public final ResultSet compute( double[][][] classSpecificScores ) {
		double[][] classificationScores = new double[classSpecificScores.length][];
		for(int i=0;i<classSpecificScores.length;i++){
			classificationScores[i] = new double[classSpecificScores[i].length];
			for(int j=0;j<classSpecificScores[i].length;j++){
				if(classSpecificScores[i][j].length != 2){
					throw new RuntimeException( "Only two classes possible for "+getClass().getSimpleName() );
				}
				classificationScores[i][j] = classSpecificScores[i][j][0] - classSpecificScores[i][j][1];
			}
		}
		return compute( classificationScores[0], classificationScores[1] );
	}

	@Override
	public final int getAllowedNumberOfClasses() {
		return 2;
	}

}
