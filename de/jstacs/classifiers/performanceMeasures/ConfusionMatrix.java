/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 *
 * For more information on Jstacs, visit http://www.jstacs.de
 */
package de.jstacs.classifiers.performanceMeasures;

import de.jstacs.io.NonParsableException;
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
	public ResultSet compute(double[] sortedScoresClass0, double[] weightClass0, double[] sortedScoresClass1, double[] weightClass1) {
		double tp =0, fn=0, fp=0, tn = 0, w;
		
		w=1;
		for( int i = 0; i < sortedScoresClass0.length; i++ ) {
			if( weightClass0 != null ) {
				w = weightClass0[i];
			}
			if( sortedScoresClass0[i] >= 0 ) {
				tp += w;
			} else {
				fn += w;
			}
		}

		w=1;
		for( int i = 0; i < sortedScoresClass1.length; i++ ) {
			if( weightClass1 != null ) {
				w = weightClass1[i];
			}
			if( sortedScoresClass0[i] >= 0 ) {
				fp += w;
			} else {
				tn += w;
			}
		}
		
		return new ResultSet(new ListResult( getName(), getName()+" for two classes.", null, 
				new NumericalResultSet( new NumericalResult[]{
						new NumericalResult( "TP", "true positives", tp ),
						new NumericalResult( "FN", "false negatives", fn )
				} ),
			    new NumericalResultSet( new NumericalResult[]{
			    		new NumericalResult( "FP", "false positives", fp ),
			    		new NumericalResult( "TN", "true negatives", tn )
			                                
			    } )
		) );
	}

	@Override
	public ResultSet compute(double[][][] classSpecificScores, double[][] weights ) {
		double[][] res = new double[classSpecificScores.length][classSpecificScores.length];
		double w;
		for(int i=0;i<classSpecificScores.length;i++){
			w = 1;
			for(int j=0;j<classSpecificScores[i].length;j++){
				if( weights[i] != null ) {
					w = weights[i][j];
				}
				int predicted = ToolBox.getMaxIndex( classSpecificScores[i][j] );
				res[i][predicted] += w;
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