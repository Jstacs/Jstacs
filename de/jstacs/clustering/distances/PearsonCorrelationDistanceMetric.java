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

import de.jstacs.utils.ToolBox;

/**
 * Implements a distance metric based on the Pearson correlation of two double vectors.
 * 
 * @author Jan Grau
 *
 */
public class PearsonCorrelationDistanceMetric extends DistanceMetric<double[]> {

	private boolean abs;
	
	/**
	 * Creates a new distance based on the Pearson correlation.
	 * @param abs if <code>true</code>, the distance is computed as \( 1 - |cor(x,y)| \), otherwise as \( 1-cor(x,y) \).
	 */
	public PearsonCorrelationDistanceMetric(boolean abs){
		this.abs = abs;
	}
	
	@Override
	public double getDistance( double[] o1, double[] o2 ) throws Exception {
		double corr = ToolBox.pearsonCorrelation( o1, o2 );
		if(abs){
			corr = Math.abs( corr );
		}
		return 1.0 - corr;
	}

}
