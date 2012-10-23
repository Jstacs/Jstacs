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
import java.util.Arrays;

import de.jstacs.algorithms.graphs.tensor.SymmetricTensor;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
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
	 * This method computes the pairwise mutual information between
	 * the sequence positions.
	 * 
	 * @param s
	 *            the {@link DataSet}
	 * @param weights
	 *            the weights for each sequence of the {@link DataSet} or
	 *            <code>null</code>
	 * 
	 * @return a matrix containing all pairwise mutual informations
	 * 
	 * @throws IllegalArgumentException
	 *             if something went wrong (e.g. the
	 *             {@link de.jstacs.data.AlphabetContainer} is not discrete, the
	 *             length is not matching with the
	 *             {@link de.jstacs.data.AlphabetContainer}, ...)
	 */
	public static double[][] getMI( DataSet s, double[] weights ) throws IllegalArgumentException {
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
		double[][] mi = getMI(s, weights);
		for( int i = 0; i < mi.length; i++ ) {
			for( int j = 0; j < mi[i].length; j++ ) {
				mi[i][j] = mi[i][j]  / Math.log(2);
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
	
	public static BufferedImage getImage(DataSet d, double[] weights, REnvironment r, double alpha, int... borders ) throws Exception {
		double[][] mi = getMI(d, weights);//TODO
		r.createMatrix( "mi", mi );
		r.createVector( "borders", borders );
		r.eval("L=dim(mi)[1];l=1:L;");
		double sum = 0;
		if( weights != null ) {
			for( int i = 0; i < weights.length; i++ ) {
				sum += weights[i];
			}
		} else {
			sum = d.getNumberOfElements();
		}
		r.eval("f=2*"+sum);
		int[] a = new int[d.getElementLength()];
		for( int i = 0; i < a.length; i++ ) {
			a[i] = (int) d.getAlphabetContainer().getAlphabetLengthAt(i);
		}
		r.createVector("a", a);
		r.voidEval( "p=matrix(NA,ncol=L,nrow=L);" );
		r.voidEval( "count=rep(0,L);" );
		r.voidEval( "correction=(L*(L-1)/2); alpha = " + alpha + ";" );
		r.voidEval( "for( i in l ) { for( j in i:L ) { "
					+ "mi[j,i]=NA;"
					+ "if( i != j ) { " 
						+ "p[j,i]=pchisq(f*mi[i,j],df=(a[i]-1)*(a[j]-1),lower.tail=F)*correction;"
						+ "if( p[j,i] > alpha ) {"
							+ "p[j,i]=NA;"
						+ "} else {"
							+ "count[i] = count[i] + 1;"
							+ "count[j] = count[j] + 1;"
						+ "}\n"
					+"}\n"
					+ "}\n }\n" );
		r.voidEval( "count[count==0]=NA;" );
		r.eval( "mini = min( mi, na.rm=T ); maxi = max( mi, na.rm=T ); lim = c(mini, maxi);c = rev( gray((1:1000)/1000) ); skala = seq( mini, maxi, length = length(c) ); m = matrix( skala, ncol = length(c) );" );
		r.eval( "miniP = min( p, na.rm=T ); maxiP = max( p, na.rm=T ); limP = c(miniP, maxiP);colP = rev( gray((800:1)/1000) ); skalaP = seq( miniP, maxiP, length = length(colP) ); mP = matrix( skalaP, ncol = length(colP) );" );
		//r.eval("seq=\"" + d.getElementAt(0) +"\"; s=rep(NA,L); for( i in l ) { s[i] = substr(seq,i,i); }\n" );
		BufferedImage b = r.plot(
				"layout( matrix( c(3,1,4,5,2,6), ncol=3, nrow=2, byrow=T ), width=c(2,11,2), height=c(4,1) ); par( las=1 );" 
									+ "image( l, l, mi, xlab = \"position\",  ylab = \"position\", zlim = lim, col = c, main = \"\" );"
									+ "image( l, l, p, add=T, col=colP );"
									+ "abline(h=borders+0.5,col=2,lwd=3);"
									+ "abline(v=borders+0.5,col=2,lwd=3);"
									
									+ "par( mar = c(5, 4, 0, 2) + 0.1);"
									+ "plot(l,count,type=\"h\",xlim=c(1,L),lwd=10,xlab=\"position\",ylab=\"# significant MI\");"
									+ "abline(v=borders+0.5,col=2,lwd=3);"
									//+ "axis( 3, l, s );"
									
									+ "par( mar = c(5, 6, 4, 0) + 0.1);"
									+ "image( 1, skala, m, col = c, xlab = \"\", ylab = \"\", main = \"MI\", axes = FALSE );"
									+ "axis( 2, seq(mini,maxi,length=5), round( seq(mini,maxi,length=5), digits = 4 ) );"
								
									+ "par( mar = c(5, 0, 4, 6) + 0.1);"
									+ "image( 1, skalaP, mP, col = colP, xlab = \"\", ylab = \"\", main = \"p-value\", axes = FALSE );"
									+ "axis( 4, seq(miniP,maxiP,length=5), signif( seq(miniP,maxiP,length=5), 2 ) );"
									+ "mtext(\"Bonferonie\", 1, 1);"
									+ "mtext(\"corrected\", 1, 2.5);"
									,

				1080,
				1080 );
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
