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