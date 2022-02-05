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

package de.jstacs.clustering.kmeans;

import java.util.Arrays;
import java.util.Random;

import de.jstacs.clustering.distances.PNorm;

/**
 * Simple K-means implementation for <code>double[]</code>.
 * 
 * @author Jens Keilwagen
 */
public class Kmeans {
	
	private PNorm ed;
	
	/**
	 * The constructor for a {@link Kmeans} instance.
	 */
	public Kmeans() {
		ed = new PNorm(2);
	}
	
	/**
	 * Computes the class assignment for the <code>data</code> given the initial <code>center</code>s.
	 * The <code>center</code>s will be updated iteratively.
	 * 
	 * @param data the data, element <code>i</code> is encodes in <code>data[i]</code>
	 * @param center the initial centers will be updated during the algorithm so that the finally it contains the centers of the clusters
	 *  
	 * @return the class assignment
	 * 
	 * @throws Exception forwarded from {@link PNorm}
	 * 
	 * @see PNorm
	 */
	public int[] cluster( double[][] data, double[][] center ) throws Exception {
		int[] res = new int[data.length];
		int[] stat = new int[center.length];
		Arrays.fill(res, -2);
		int changed=-1;
		while( true ) {
			//update class assignment
			changed=0;
			for( int i = 0; i<data.length; i++ ) {
				int old=res[i];
				res[i] = -1;
				double bestD=Double.POSITIVE_INFINITY;
				for( int j = 0; j<center.length; j++ ) {
					double current = ed.getDistance(data[i], center[j]);
					if( current < bestD ) {
						bestD=current;
						res[i]=j;
					}
				}
				if( old != res[i] ) changed++;
			}
			//check if something changed
			if( changed == 0 ) break;
			//update centers
			Arrays.fill(stat, 0);
			for( int i = 0; i<center.length; i++ ) {
				Arrays.fill( center[i], 0 );
			}
			for( int i = 0; i<data.length; i++ ) {
				stat[res[i]]++;
				for( int j = 0; j<data[i].length; j++ ) {
					center[res[i]][j] += data[i][j];
				}
			}
			for( int i = 0; i<center.length; i++ ) {
				for( int j = 0; j<center[i].length; j++ ) {
					center[i][j] /= stat[i];
				}
			}
		}
		return res;
	}
	
	/**
	 * K-means for given number <code>k</code> of clusters and given number <code>n</code> of restarts.
	 * 
	 * @param data the data
	 * @param k the number of clusters
	 * @param n the number of (re-)starts
	 * 
	 * @return the cluster assignment
	 * 
	 * @throws Exception forwarded from {@link PNorm}
	 * 
	 * @see PNorm
	 */
	public int[] cluster( double[][] data, int k, int n ) throws Exception {
		double[][] center = new double[k][data[0].length];
		int[] best = null;
		double sse=Double.POSITIVE_INFINITY;
		Random r = new Random();
		double[] maxmin = new double[data.length];
		for( int i = 0; i < n; i++ ) {
			//choose random clusters using Maxmin: https://doi.org/10.1016/j.patcog.2019.04.014
			Arrays.fill(maxmin, Double.POSITIVE_INFINITY);
			int first = r.nextInt(data.length);
			System.arraycopy(data[first],0,center[0],0,data[first].length);
			for( int j = 1; j < k; j++ ) {
				int max=-1;
				for( int l = 0; l < data.length; l++ ) {
					double current = ed.getDistance(center[j-1], data[l]);
					if( current < maxmin[l] ) {
						maxmin[l] = current;
					}
					if( max<0 || maxmin[l] > maxmin[max] ) {
						max=l;
					}
				}
				System.arraycopy(data[max],0,center[j],0,data[max].length);
			}
			
			//compute solution
			int[] assignment = cluster(data, center);
			
			//determine sse
			double current = getSSE(center, assignment, data);
			if( current < sse ) {
				sse = current;
				best = assignment;
			}
		}
		return best;
	}
	
	/**
	 * Computes the sum of squared errors (SSE).
	 * 
	 * @param center the cluster centers
	 * @param assignment the cluster assignment
	 * @param data the data
	 * 
	 * @return SSE
	 * 
	 * @throws Exception forwarded from {@link PNorm}
	 * 
	 * @see PNorm
	 */
	double getSSE( double[][] center, int[] assignment, double[][] data ) throws Exception {
		double sse = 0;
		for( int i = 0; i<data.length; i++ ) {
			double dist = ed.getDistance(data[i], center[assignment[i]] );
			sse += dist*dist;
		}
		return sse;
	}
	
	/**
	 * This method allows to rescale the data, which might be useful ahead of clustering.
	 * 
	 * @param x the data
	 * 
	 * @return the rescaled data
	 */
	public static double[][] rescale( double[][] x ) {
		//mean
		double[] mean = new double[x[0].length];
		for( int i = 0; i < x.length; i++ ) {
			for( int j = 0; j < x[0].length; j++ ) {
				mean[j] += x[i][j];
			}
		}
		for( int j = 0; j < mean.length; j++ ) {
			mean[j] /= x.length;
		}

		//sd
		double[] sd = new double[x[0].length];
		for( int i = 0; i < x.length; i++ ) {
			for( int j = 0; j < x[0].length; j++ ) {
				double v = x[i][j]-mean[j];
				sd[j] += v*v;
			}
		}
		for( int j = 0; j < mean.length; j++ ) {
			sd[j] = Math.sqrt(sd[j])/(x.length-1);
		}
		
		double[][] rescaled = new double[x.length][x[0].length];
		for( int i = 0; i < x.length; i++ ) {
			for( int j = 0; j < x[0].length; j++ ) {
				rescaled[i][j] = (x[i][j]-mean[j])/sd[j];
			}
		}
		
		return rescaled;
	}
}