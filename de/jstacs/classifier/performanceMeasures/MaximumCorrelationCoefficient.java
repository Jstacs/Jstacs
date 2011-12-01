package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;


public class MaximumCorrelationCoefficient extends TwoClassAbstractPerformanceMeasure implements NumericalPerformanceMeasure {

	public MaximumCorrelationCoefficient() {
	}

	public MaximumCorrelationCoefficient( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	public String getName() {
		return "Maximum correlation coefficient";
	}

	@Override
	public NumericalResultSet compute( double[] scoresClass0, double[] scoresClass1 ) {
		int i = 1, j = 1, d = scoresClass1.length, m = scoresClass0.length;
		double[] current = new double[2];
		double[] erg = new double[]{ Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };

		//find second smallest value => first threshold current[0]

		if( scoresClass0[0] == scoresClass1[0] ) {
			while( j < d && scoresClass1[j - 1] == scoresClass1[j] ) {
				j++;
			}
			while( i < m && scoresClass0[i - 1] == scoresClass0[i] ) {
				i++;
			}
			if( i != m && j != d ) {
				current[0] = Math.min( scoresClass1[j], scoresClass0[i] );
				// decoys_lr[j-1] = motifs_lr[i-1] < c <= Math.min( decoys_lr[j], motifs_lr[i] )
			} else if( i == m && j == d ) {
				throw new IllegalArgumentException( "All likelihood-ratios had the same values for both samples." );
			} else if( i == m ) // j != d
			{
				current[0] = scoresClass1[j];
				// motifs_lr[i-1] = decoys_lr[j-1] < c = decoys_lr[j], i == m
			} else
			// if( j == d && i != m )
			{
				current[0] = scoresClass0[i];
				// decoys_lr[j-1] = motifs_lr[i-1] < c = motifs_lr[i] , j == d
			}
		} else {
			if( scoresClass0[0] < scoresClass1[0] ) {
				while( i < m && scoresClass0[i - 1] == scoresClass0[i] ) {
					i++;
				}
				j = 0;
				if( i == m ) {
					current[0] = scoresClass1[j];
					// motifs_lr[i-1] < c == decoys_lr[j], i == m
				} else {
					if( scoresClass0[i] <= scoresClass1[j] ) {
						current[0] = scoresClass0[i];
						// motifs_lr[i-1] < c == motifs_lr[i] <= decoys_lr[j]
					} else {
						current[0] = scoresClass1[j];
						// motifs_lr[i-1] < c == decoys_lr[j] < motifs_lr[i]
					}
				}
			} else // ( motifs_lr[0] > decoys_lr[0] )
			{
				while( j < d && scoresClass1[j - 1] == scoresClass1[j] ) {
					j++;
				}
				i = 0;
				if( j == d ) {
					current[0] = scoresClass0[i];
					//decoys_lr[j-1] < c == motis_lr[i], j == d
				} else {
					if( scoresClass1[j] <= scoresClass0[i] ) {
						current[0] = scoresClass1[j];
						// decoys_lr[j-1] < c = decoys_lr[j] <= motifs_lr[i]
					} else {
						current[0] = scoresClass0[i];
						// motifs_lr[i-1] < c = motifs_lr[i] < decoys_lr[j]
					}
				}
			}
		}
		// forall 0 <= l < i <= k < m:  motifs_lr[l] < c <= motifs_lr[k]
		// forall 0 <= l < j <= k < d:  decoys_lr[l] < c <= decoys_lr[k]

		while( i < m || j < d ) {
			// compute CC
			current[1] = getCorrelationCoefficient( m - i, d - j, i, j );
			if( current[1] > erg[1] ) {
				erg[0] = current[0];
				erg[1] = current[1];
			}

			//find next threshold
			while( j < d && scoresClass1[j] == current[0] ) {
				j++;
			}
			while( i < m && scoresClass0[i] == current[0] ) {
				i++;
			}
			if( i < m && j < d ) {
				if( scoresClass0[i] <= scoresClass1[j] ) {
					current[0] = scoresClass0[i];
				} else // if( motifs_lr[ i ] > decoys_lr[ j ] )
				{
					current[0] = scoresClass1[j];
				}
			} else if( i < m ) // j == d
			{
				current[0] = scoresClass0[i];
			} else if( j < d ) // i == m
			{
				current[0] = scoresClass1[j];
			}
			// else // i == m && j == d => the outer loop will break
		}
		return new NumericalResultSet(new NumericalResult[]{
				new NumericalResult( getName(), "The value of the " + getName().toLowerCase(), erg[1] ),
				new NumericalResult( "Threshold", "Threshold for the "  + getName().toLowerCase(), erg[0] )
		});
	}
	
	public NumericalResultSet compute( double[][][] classSpecificScores ) {
		return (NumericalResultSet) super.compute( classSpecificScores );
	}
	
	// double since otherwise an overflow is easily possible
	private static double getCorrelationCoefficient( double tp, double fp, double fn, double tn ) {
		return ( tp * tn - fn * fp ) / Math.sqrt( ( tp + fn ) * ( tn + fp ) * ( tp + fp ) * ( tn + fn ) );
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList();
	}
}
