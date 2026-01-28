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
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.utils.ToolBox;

/**
 * This class implements the classification rate, i.e. \( \frac{TP + TN}{N} \).
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
	
	public NumericalResultSet compute(double[] sortedScoresClass0, double[] sortedScoresClass1) {
		return compute( sortedScoresClass0, null, sortedScoresClass1, null) ;
	}

	public NumericalResultSet compute(double[][][] classSpecificScores) {
		return compute( classSpecificScores, null );
	}

	@Override
	public NumericalResultSet compute( double[] sortedScoresClass0, double[] weightsClass0, double[] sortedScoresClass1, double[] weightsClass1 ) {
		double corr = 0, fals = 0, w;
		
		w = 1;
		for( int i = 0; i < sortedScoresClass0.length; i++ ) {
			if( weightsClass0 != null ) {
				w = weightsClass0[i];
			}
			if( sortedScoresClass0[i] >= 0 ) {
				corr += w;	
			} else {
				fals +=w;
			}
		}
		
		w = 1;
		for( int i = 0; i < sortedScoresClass1.length; i++ ) {
			if( weightsClass1 != null ) {
				w = weightsClass1[i];
			}
			if( sortedScoresClass1[i] < 0 ) {
				corr += w;	
			} else {
				fals +=w;
			}
		}
		return getResult( 2, corr, fals );
	}

	@Override
	public NumericalResultSet compute( double[][][] classSpecificScores, double[][] weights ) {
		double corr = 0, fals = 0, w;
		for(int i=0;i<classSpecificScores.length;i++){
			w = 1;
			for(int j=0;j<classSpecificScores[i].length;j++){
				if( weights != null && weights[i] != null ) {
					w = weights[i][j];
				}
				if(ToolBox.getMaxIndex( classSpecificScores[i][j] ) == i){
					corr += w;
				}else{
					fals += w;
				}
			}
		}
		return getResult( classSpecificScores.length, corr, fals );
	}
	
	private NumericalResultSet getResult( int classes, double corr, double fals ) {
		return new NumericalResultSet( new NumericalResult( getName(), getName() + " for "+classes+" classes.", corr/(corr+fals) ) );
	}

	@Override
	public int getAllowedNumberOfClasses() {
		return 0;
	}
}