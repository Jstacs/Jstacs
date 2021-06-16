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

package de.jstacs.clustering.distances;

/**
 * This class implements a p-norm. For p=1, this results in a Manhattan distance and in an Euclidean distance for p=2
 *  
 * @author Jens Keilwagen
 */
public class PNorm extends DistanceMetric<double[]> {

	double p, q;
	
	/**
	 * Creates a new {@link PNorm} instance.
	 * 
	 * @param p the value used for the norm
	 */
	public PNorm( double p) {
		if( p < 1 ) {
			throw new IllegalArgumentException("p has to be larger or equal to 1.");
		}
		this.p = p;
		this.q = 1d/p;
	}

	@Override
	public double getDistance(double[] o1, double[] o2) throws Exception {
		if( o1==null || o2 == null || o1.length!=o2.length ) {
			throw new IllegalArgumentException("The object have to have the same length.");
		}
		double diff, sum=0;
		for( int i = 0; i < o1.length; i++ ) {
			diff = Math.abs(o1[i]-o2[i]);
			sum += Math.pow(diff,p);
		}
		return Math.pow(sum,q);
	}
}