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
	 * @see de.jstacs.classifier.performanceMeasures.AbstractPerformanceMeasure#getName()
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
				new NumericalResult( "Maximum "+ getSpecificName(), "The maximal value of the " + getSpecificName().toLowerCase(), erg[1] ),
				new NumericalResult( "Threshold for the maximum "  + getSpecificName().toLowerCase(), "", erg[0] )
		});
	}

	public NumericalResultSet compute(double[][][] classSpecificScores) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}
}