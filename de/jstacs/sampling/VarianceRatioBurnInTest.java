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

package de.jstacs.sampling;

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;

/**
 * In this class the Variance-Ratio method of Gelman and Rubin is implemented to
 * test the length of the burn-in phase of a multi-chain Gibbs Sampling (number
 * of chains >2). The number of initial iterations is calculated by comparing
 * the variances between the different chains with the variances within the
 * chains. The method returns the same length of the burn-in for all sampled
 * chains.
 * 
 * @author Berit Haldemann
 */
public class VarianceRatioBurnInTest extends AbstractBurnInTest {

	/**
	 * The threshold for testing the end of the burn-in phase with the
	 * Variance-Ratio burn-in test.
	 */
	private double threshold;

	/**
	 * Creates a new {@link VarianceRatioBurnInTest} instance.
	 * 
	 * @param parameters
	 *            set of parameters
	 *            
	 * @throws CloneNotSupportedException if the parameters can not be cloned
	 */
	public VarianceRatioBurnInTest( VarianceRatioBurnInTestParameterSet parameters ) throws CloneNotSupportedException {
		super( parameters );
		this.threshold = parameters.getThreshold();

	}

	/**
	 * The standard constructor for the {@link de.jstacs.Storable} interface.
	 * Creates a new {@link VarianceRatioBurnInTest} instance out of a
	 * {@link StringBuffer}.
	 * 
	 * @param rep
	 *            the {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} can not be parsed
	 */
	public VarianceRatioBurnInTest( StringBuffer rep ) throws NonParsableException {
		super( rep );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.mixture.gibbssampling.AbstractBurnInTest#getXMLTag()
	 */
	@Override
	protected String getXMLTag() {
		return getClass().getSimpleName();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.mixture.gibbssampling.AbstractBurnInTest#getFurtherInformation()
	 */
	@Override
	protected StringBuffer getFurtherInformation() {
		StringBuffer furtherinf = new StringBuffer( 2000 );
		XMLParser.appendObjectWithTags( furtherinf, threshold, "threshold" );
		return furtherinf;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.mixture.gibbssampling.AbstractBurnInTest#setFurtherInformation(java.lang.StringBuffer)
	 */
	@Override
	protected void setFurtherInformation( StringBuffer xml ) throws NonParsableException {

		threshold = XMLParser.extractObjectForTags( xml, "threshold", double.class );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.mixture.gibbssampling.BurnInTest#getInstanceName()
	 */
	public String getInstanceName() {
		return "Variance-Ratio burn-in test of Gelman and Rubin for " + values.length + " different chains with threshold " + threshold;
	}

	/**
	 * Computes and returns the length of the burn-in phase given by the
	 * Variance-Ratio burn-in test. To get an effective test the number of Gibbs
	 * iterations should be greater than 250. Fewer iterations may lead to
	 * inaccuracy in the calculation of the different required variances.
	 * 
	 * @see de.jstacs.sampling.AbstractBurnInTest#computeLengthOfBurnIn()
	 */
	@Override
	protected int computeLengthOfBurnIn() {

		if( values[0].length() < 250 ) {

			return values[0].length();

		} else {
			int m = values.length;
			DoubleList meanOfchains = new DoubleList( m );
			int minL = values[0].length();
			int maxL = values[0].length();
			for(int i=1;i<values.length;i++){
				if(values[i].length() < minL){
					minL = values[i].length();
				}
				if(values[i].length() > maxL){
					maxL = values[i].length();
				}
			}
			for( int it = 250; it < minL; it = it + 2 ) {

				// defining the window for the calculations
				int n = it / 2;
				int start = it - n;

				// calculation of the means of the different chains
				meanOfchains.clear();
				for( int i = 0; i < m; i++ ) {
					meanOfchains.add( values[i].mean( start, it ) );
				}

				// calculation of the means of the calculated means
				double meanOfChainMeans = meanOfchains.mean( 0, meanOfchains.length() );

				// calculation of the measure B
				double b = 0;
				for( int i = 0; i < m; i++ ) {
					b += ( meanOfchains.get( i ) - meanOfChainMeans ) * ( meanOfchains.get( i ) - meanOfChainMeans );
				}
				b = b / ( (double)m - 1d );

				// calculation of the measure W
				double w = 0;
				for( int i = 0; i < m; i++ ) {
					for( int j = start; j < it; j++ ) {
						w += ( values[i].get( j ) - meanOfchains.get( i ) ) * ( values[i].get( j ) - meanOfchains.get( i ) );
					}
				}
				w = w / ( (double)m * ( (double)n - 1d ) );

				// calculation of the estimator of sigma
				double sigma = ( 1d - ( 1d / (double)n ) ) * w;
				sigma += ( ( 1d + ( 1d / (double)m ) ) * b );

				// calculation of the potential scale reduction factor R
				double r = Math.sqrt( sigma / w );

				// test if the burn-in is ended
				if( r < threshold ) {
					return it;
				}
			}
			return maxL;
		}
	}
}
