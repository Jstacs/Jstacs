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

import java.util.ArrayList;

import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.io.NonParsableException;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.ToolBox;

/**
 * This class implements the Receiver Operating Characteristics curve and the area under the curve.
 * The Receiver Operating Characteristics curve is the plot of sensitivity ({@latex.inline $\\frac{TP}{TP+FN}$}) 
 * against the false positive rate ( {@latex.inline $\\frac{FP}{FP+TN}$}) for all possible classification thresholds.
 * 
 * If you are only interested in the area under this curve, you can use {@link AucROC} instead.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ROCCurve extends AbstractTwoClassPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link ROCCurve}.
	 */
	public ROCCurve() {
		super();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ROCCurve} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ROCCurve} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public ROCCurve( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/**
	 * The name of the performance measure return by {@link #getName()} 
	 */
	public static final String NAME = "Receiver Operating Characteristic curve";
	
	@Override
	public String getName() {
		return NAME;
	}

	@Override
	public ResultSet compute( double[] sortedScoresClass0, double[] weightsClass0, double[] sortedScoresClass1, double[] weightsClass1 ) {
		ArrayList<double[]> list = null;
		if( !(this instanceof NumericalPerformanceMeasure) ){
			list = new ArrayList<double[]>();
		}
		
		int i = 0, j = 0, d = sortedScoresClass1.length, m = sortedScoresClass0.length;
		double pos, neg, fn = 0, tn = 0;//TODO
		if( weightsClass0 == null ) {
			pos = m;
		} else {
			pos = ToolBox.sum( 0, m, weightsClass0 );
		}
		if( weightsClass1 == null ) {
			neg = d;
		} else {
			neg = ToolBox.sum( 0, d, weightsClass1 );
		}
		double erg = 0, help1, help2;
		double[] p = new double[]{ 1, 1 };

		if( list != null ) {
			list.add( p.clone() );
		}
		//which class defines the threshold
		boolean unique, fromMotif = false;
		if( sortedScoresClass0[i] == sortedScoresClass1[j] ) {
			unique = false;
		} else {
			unique = true;
			if( sortedScoresClass0[i] < sortedScoresClass1[j] ) {
				fromMotif = true;
			} else {
				fromMotif = false;
			}
		}
		while( i < m && j < d ) {
			// find next possible threshold
			// discard values that are not interesting
			if( unique ) {
				if( fromMotif ) {
					while( i < m && sortedScoresClass0[i] < sortedScoresClass1[j] ) {
						fn += getWeight( weightsClass0, i );
						i++;
					}
				} else {
					while( j < d && sortedScoresClass0[i] > sortedScoresClass1[j] ) {
						tn += getWeight( weightsClass1, j );
						j++;
					}
				}
			} else {
				while( i + 1 < m && sortedScoresClass0[i] == sortedScoresClass0[i + 1] ) {
					fn += getWeight( weightsClass0, i );
					i++;
				}
				while( j + 1 < d && sortedScoresClass1[j] == sortedScoresClass1[j + 1] ) {
					tn += getWeight( weightsClass1, j );
					j++;
				}
				fn += getWeight( weightsClass0, i );
				tn += getWeight( weightsClass1, j );
				i++;
				j++;
			}

			// erg += height * width
			help1 = (neg-tn) / neg;
			help2 = (pos-fn) / pos;
			erg += ( p[1] + help2 ) / 2d * ( p[0] - help1 );
			p[0] = help1;
			p[1] = help2;
			if( list != null ) {
				list.add( p.clone() );
			}

			if( i < m && j < d ) {
				//next
				if( sortedScoresClass0[i] == sortedScoresClass1[j] ) {
					unique = false;
				} else {
					unique = true;
					if( sortedScoresClass0[i] < sortedScoresClass1[j] ) {
						fromMotif = true;
					} else {
						fromMotif = false;
					}
				}
			}
		}
		if( list != null ) {
			list.add( new double[]{ 0, 0 } );
		}
		NumericalResult auc = new NumericalResult( "AUC-ROC", getName(), erg );
		if(list == null){
			return new NumericalResultSet( auc );
		}else{
			return new ResultSet( new Result[]{
               auc,
               new AbstractScoreBasedClassifier.DoubleTableResult(getName(), getName(), list)
			} );
		}
	}

}
