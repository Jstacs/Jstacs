package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.utils.ToolBox;


public class ClassificationRate extends AbstractPerformanceMeasure implements NumericalPerformanceMeasure {

	public ClassificationRate() {
	}

	public ClassificationRate( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Classification rate";
	}

	@Override
	public NumericalResultSet compute( double[] classificationScoresFg, double[] classificationScoresBg ) {
		int i = 0, m = classificationScoresFg.length;
		while( i < m && classificationScoresFg[i] < 0 ) {
			i++;
		}

		int d = classificationScoresBg.length, j = d - 1;
		while( j >= 0 && classificationScoresBg[j] >= 0 ) {
			j--;
		}

		return new NumericalResultSet(new NumericalResult("Classification rate","Classification rate for two classes.",( ( m - i ) + j + 1 ) / (double)( m + d )));
	}

	@Override
	public NumericalResultSet compute( double[][][] classSpecificScores ) {
		int corr = 0, fals = 0;
		for(int i=0;i<classSpecificScores.length;i++){
			for(int j=0;j<classSpecificScores[i].length;j++){
				if(ToolBox.getMaxIndex( classSpecificScores[i][j] ) == i){
					corr++;
				}else{
					fals++;
				}
			}
		}
		return new NumericalResultSet(new NumericalResult(getName(),getName() + " for "+classSpecificScores[0][0].length+" classes.",(double)corr/(double)(corr+fals)));
	}

	@Override
	public int getAllowedNumberOfClasses() {
		return 0;
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList();
	}

}
