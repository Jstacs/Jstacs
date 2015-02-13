package de.jstacs.clustering.distances;

import de.jstacs.utils.ToolBox;


public class PearsonCorrelationDistanceMetric extends DistanceMetric<double[]> {

	private boolean abs;
	
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
