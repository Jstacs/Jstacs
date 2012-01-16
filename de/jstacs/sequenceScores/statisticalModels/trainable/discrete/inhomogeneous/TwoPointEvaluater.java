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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous;

import java.awt.image.BufferedImage;

import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.graphs.tensor.SymmetricTensor;
import de.jstacs.data.DataSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.utils.REnvironment;

/**
 * This class is for visualizing two point dependency between sequence
 * positions.
 * 
 * @author Jens Keilwagen
 */
public class TwoPointEvaluater {

	private static final byte order = 1;

	/**
	 * This method computes the pairwise mutual information (in bits) between
	 * the sequence positions.
	 * 
	 * @param s
	 *            the {@link DataSet}
	 * @param weights
	 *            the weights for each sequence of the {@link DataSet} or
	 *            <code>null</code>
	 * 
	 * @return a matrix containing all pairwise mutual informations (in bits)
	 * 
	 * @throws IllegalArgumentException
	 *             if something went wrong (e.g. the
	 *             {@link de.jstacs.data.AlphabetContainer} is not discrete, the
	 *             length is not matching with the
	 *             {@link de.jstacs.data.AlphabetContainer}, ...)
	 */
	public static double[][] getMIInBits( DataSet s, double[] weights ) throws IllegalArgumentException {
		int l = s.getElementLength(), i = 0;
		int[] j = new int[1];
		StructureLearner sl = new StructureLearner( s.getAlphabetContainer(), l );
		SymmetricTensor tensor = null;
		try {
			tensor = sl.getTensor( s, weights, order, LearningType.ML_OR_MAP );
		} catch ( WrongAlphabetException w ) {
			//does not happen
			throw new IllegalArgumentException();
		}
		double sum = 0;
		if( weights != null ) {
			for( ; i < weights.length; i++ ) {
				sum += weights[i];
			}
		} else {
			sum = s.getNumberOfElements();
		}
		double[][] mi = new double[l][l];
		for( ; i < l; i++ ) {
			for( j[0] = i + 1; j[0] < l; j[0]++ ) {
				mi[i][j[0]] = mi[j[0]][i] = tensor.getValue( order, i, j ) / sum;
				//System.out.println( mi[i][j[0]] );
			}
		}
		return mi;
	}

	/**
	 * This method can be used to create an image of a mutual information
	 * matrix.
	 * 
	 * @param mi
	 *            the matrix of mutual information
	 * @param r
	 *            the R environment
	 * 
	 * @return an image of the matrix
	 * 
	 * @throws Exception
	 *             if an {@link Exception} is thrown from {@link REnvironment}
	 * 
	 * @see TwoPointEvaluater#getMIInBits(DataSet, double[])
	 */
	public static BufferedImage getImage( double[][] mi, REnvironment r ) throws Exception {
		r.createMatrix( "mi", mi );
		r.eval( "mini = min( mi ); maxi = max( mi );" );
		r.eval( "lim = c(mini, maxi);c = rev( gray((1:1000)/1000) ); l=1:dim(mi)[1];" );
		r.eval( "skala = seq( mini, maxi, length = 1000 ); m = matrix( skala, ncol = 1000 );" );
		BufferedImage b = r.plot( "layout( matrix( c(1,2), ncol=2 ), width=c(12,3) ); par( las=1 );" + "image( l, l, mi, xlab = \"position\",  ylab = \"position\", zlim = lim, col = c, main = \"mutual information\" );"
									+ "par( mar = c(5, 0, 4, 6) + 0.1);"
									+ "image( 1, skala, m, col = c, xlab = \"\", ylab = \"\", main = \"scale\", axes = FALSE );"
									+ "axis( 4, seq(mini,maxi,length=5), round( seq(mini,maxi,length=5), digits = 4 ) )",
				1080,
				864 );
		//r.eval( "layout( matrix( 1 ), height=1, width=1 );" );
		return b;
	}

	/**
	 * This method can be used to determine the maximal value of the matrix of
	 * mutual informations.
	 * 
	 * @param mi
	 *            the matrix of mutual informations
	 * 
	 * @return the maximal value
	 */
	public static double getMax( double[][] mi ) {
		double max = Double.NEGATIVE_INFINITY;
		int i = 0, j;
		for( ; i < mi.length; i++ ) {
			for( j = 0; j < mi[i].length; j++ ) {
				if( mi[i][j] > max ) {
					max = mi[i][j];
				}
			}
		}
		return max;
	}
}
