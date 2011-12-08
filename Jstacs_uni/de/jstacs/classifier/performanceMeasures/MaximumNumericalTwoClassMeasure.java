package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;

/**
 * This class prepares everything for an easy implementation of a maximum of any numerical performance measure.
 * 
 * @author Jens Keilwagen
 */
public abstract class MaximumNumericalTwoClassMeasure extends TwoClassAbstractPerformanceMeasure implements NumericalPerformanceMeasure {

	public MaximumNumericalTwoClassMeasure() {}

	public MaximumNumericalTwoClassMeasure(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	protected abstract double getMeasure( double tp, double fp, double fn, double tn );
	
	@Override
	public NumericalResultSet compute(double[] sortedScoresClass0, double[] sortedScoresClass1) {
		int i = 0, j = 0, d = sortedScoresClass1.length, m = sortedScoresClass0.length;
		double[] current = new double[2];
		current[0] = Math.min( sortedScoresClass0[i], sortedScoresClass1[j] );
		double[] erg = new double[]{ Double.NaN, Double.NaN };
		
		while( i < m || j < d ) {
			// compute CC
			current[1] = getMeasure( m - i, d - j, i, j );
			if( Double.isNaN(erg[1]) || current[1] > erg[1] ) {
				erg[0] = current[0];
				erg[1] = current[1];
			}
	
			//find next threshold
			while( j < d && sortedScoresClass1[j] == current[0] ) {
				j++;
			}
			while( i < m && sortedScoresClass0[i] == current[0] ) {
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
				new NumericalResult( "Maximum "+ getName(), "The maximal value of the " + getName().toLowerCase(), erg[1] ),
				new NumericalResult( "Threshold for the maximum "  + getName().toLowerCase(), "", erg[0] )
		});
	}

	public NumericalResultSet compute(double[][][] classSpecificScores) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}
}