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
	 * @param abs if <code>true</code>, the distance is computed as {@latex.inline $1 - |cor(x,y)|$}, otherwise as {@latex.inline $1-cor(x,y)$}.
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
