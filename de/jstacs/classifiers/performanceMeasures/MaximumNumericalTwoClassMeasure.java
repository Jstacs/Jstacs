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
 * This class prepares everything for an easy implementation of a maximum of any numerical performance measure.
 * 
 * @author Jens Keilwagen
 */
public abstract class MaximumNumericalTwoClassMeasure extends AbstractNumericalTwoClassPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link MaximumNumericalTwoClassMeasure}.
	 */
	protected MaximumNumericalTwoClassMeasure() {}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link MaximumNumericalTwoClassMeasure} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MaximumFMeasure} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	protected MaximumNumericalTwoClassMeasure(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.classifiers.performanceMeasures.AbstractPerformanceMeasure#getName()
	 */
	@Override
	public final String getName() {
		return "Maximum " + getMeasureName();
	}
	
	/**
	 * This method returns a short name of the measure without any parameters.
	 * @return a short name of the measure without any parameters
	 */
	protected abstract String getMeasureName();
	
	/**
	 * This method returns a specific name of the measure including any parameters.
	 * @return a specific name of the measure including any parameters
	 */
	protected abstract String getSpecificName();

	/**
	 * This measure compute the measure for a given confusion matrix
	 * 
	 * @param tp true positives
	 * @param fp false positives
	 * @param fn false negatives
	 * @param tn true negative
	 * 
	 * @return the value of the measure
	 */
	protected abstract double getMeasure( double tp, double fp, double fn, double tn );
	
	@Override
	public NumericalResultSet compute(double[] sortedScoresClass0, double[] weightsClass0, double[] sortedScoresClass1, double[] weightsClass1 ) {
		int i = 0, j = 0, d = sortedScoresClass1.length, m = sortedScoresClass0.length;
		double pos, neg, tn = 0, fn = 0, w;
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
		double tp = pos, fp = neg;
		double[] current = new double[2];
		current[0] = Math.min( sortedScoresClass0[i], sortedScoresClass1[j] );
		double[] erg = new double[]{ Double.NaN, Double.NaN };
		
		while( i < m || j < d ) {
			// compute CC
			current[1] = getMeasure( tp, fp, fn, tn );
			if( Double.isNaN(erg[1]) || current[1] > erg[1] ) {
				erg[0] = current[0];
				erg[1] = current[1];
			}
	
			//find next threshold
			while( j < d && sortedScoresClass1[j] == current[0] ) {
				w=getWeight( weightsClass1, j );
				tn += w;
				fp -= w;
				j++;
			}
			while( i < m && sortedScoresClass0[i] == current[0] ) {
				w = getWeight( weightsClass0, i );
				tp -= w;
				fn += w;
				i++;
			}
			if( i < m && j < d ) {
				if( sortedScoresClass0[i] <= sortedScoresClass1[j] ) {
					current[0] = sortedScoresClass0[i];
				} else // if( motifs_lr[ i ] > decoys_lr[ j ] )
				{
					current[0] = sortedScoresClass1[j];
				}
			} else if( i < m ) // j == d
			{
				current[0] = sortedScoresClass0[i];
			} else if( j < d ) // i == m
			{
				current[0] = sortedScoresClass1[j];
			}
			// else // i == m && j == d => the outer loop will break
		}
		if( Double.isNaN(erg[1]) ) {
			throw new IllegalArgumentException( "The measure could not be properly computed as it return always NaN. The arrays looked like: "
					+ "[" + sortedScoresClass0[0] + " .. " + sortedScoresClass0[sortedScoresClass0.length-1] + "] and "
					+ "[" + sortedScoresClass1[0] + " .. " + sortedScoresClass1[sortedScoresClass1.length-1] + "].");
		}
		return new NumericalResultSet(new NumericalResult[]{
				new NumericalResult( "Maximum "+ getSpecificName(), "The maximal value of the " + getSpecificName().toLowerCase(), erg[1] ),
				new NumericalResult( "Threshold for the maximum "  + getSpecificName().toLowerCase(), "", erg[0] )
		});
	}
}