package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.utils.ToolBox;

/**
 * This class implements the classification rate.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ClassificationRate extends AbstractPerformanceMeasure implements NumericalPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link ClassificationRate}.
	 */
	public ClassificationRate() {}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ClassificationRate} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ClassificationRate} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public ClassificationRate( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Classification rate";
	}

	@Override
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] sortedScoresClass1 ) {
		int i = 0, m = sortedScoresClass0.length;
		while( i < m && sortedScoresClass0[i] < 0 ) {
			i++;
		}

		int d = sortedScoresClass1.length, j = d - 1;
		while( j >= 0 && sortedScoresClass1[j] >= 0 ) {
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
