package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.ToolBox;

/**
 * This class implements the performance measure confusion matrix.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ConfusionMatrix extends AbstractPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link ConfusionMatrix}.
	 */
	public ConfusionMatrix() {
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ConfusionMatrix} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ConfusionMatrix} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public ConfusionMatrix(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	@Override
	public ResultSet compute(double[] sortedScoresClass0, double[] sortedScoresClass1) {
		int i = 0, m = sortedScoresClass0.length;
		while( i < m && sortedScoresClass0[i] < 0 ) {
			i++;
		}

		int d = sortedScoresClass1.length, j = d - 1;
		while( j >= 0 && sortedScoresClass1[j] >= 0 ) {
			j--;
		}
		return new ResultSet(new ListResult( getName(), getName()+" for two classes.", null, 
				new NumericalResultSet( new NumericalResult[]{new NumericalResult( "TP", "true positives", sortedScoresClass0.length-i ),
				                            new NumericalResult( "FN", "false negatives", i )
				} ),
			    new NumericalResultSet( new NumericalResult[]{new NumericalResult( "FP", "false positives", j+1 ),
			                                new NumericalResult( "TN", "true negatives", d-(j+1) )
			                                
			    } )
		) );
	}

	@Override
	public ResultSet compute(double[][][] classSpecificScores) {
		int[][] res = new int[classSpecificScores.length][classSpecificScores.length];
		for(int i=0;i<classSpecificScores.length;i++){
			for(int j=0;j<classSpecificScores[i].length;j++){
				int predicted = ToolBox.getMaxIndex( classSpecificScores[i][j] );
				res[i][predicted]++;
			}
		}
		NumericalResultSet[] sets = new NumericalResultSet[classSpecificScores.length];
		for(int i=0;i<res.length;i++){
			NumericalResult[] temp = new NumericalResult[res[i].length];
			for(int j=0;j<res[i].length;j++){
				temp[j] = new NumericalResult( i+"/"+j, "correct class: "+i+", predicted class: "+j, res[i][j] );
			}
			sets[i] = new NumericalResultSet( temp );
		}
		return new ResultSet( new ListResult( getName(), getName()+" for "+classSpecificScores.length+" classes" , null, sets ) );
	}

	@Override
	public int getAllowedNumberOfClasses() {
		return 0;
	}

	@Override
	public String getName() {
		return "Confusion matrix";
	}
}