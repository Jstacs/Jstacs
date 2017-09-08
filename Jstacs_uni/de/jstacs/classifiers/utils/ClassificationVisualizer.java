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

package de.jstacs.classifiers.utils;

import java.util.ArrayList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.classifiers.AbstractClassifier;
import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ImageResult;
import de.jstacs.utils.REnvironment;

/**
 * This class enables you to visualize some classifier results.
 * 
 * @author Jens Keilwagen
 */
public class ClassificationVisualizer {

	private ClassificationVisualizer() {}

	private static String getPlotScoresCmd( AbstractScoreBasedClassifier cl, DataSet class0, DataSet class1, REnvironment e, int bins,
			double density, String plotOptions ) throws Exception {
		StringBuffer b = new StringBuffer( 100000 );
		if( cl.getNumberOfClasses() != 2 ) {
			throw new OperationNotSupportedException( "This method is only possible for 2-class-classifiers." );
		}

		e.createVector( "sample0", cl.getScores( class0 ) );
		e.createVector( "sample1", cl.getScores( class1 ) );

		b.append( "min0 = min( sample0 ); max0 = max( sample0 );\n" );
		b.append( "min1 = min( sample1 ); max1 = max( sample1 );\n" );
		b.append( "min = min( min0, min1 ); max = max( max0, max1 );\n" );

		// System.out.println( "[" + e.eval( "min0" ).asDouble() + ", " + e.eval( "max0" ).asDouble() + "]" );
		// System.out.println( "[" + e.eval( "min1" ).asDouble() + ", " + e.eval( "max1" ).asDouble() + "]" );
		// System.out.println( "[" + e.eval( "min" ).asDouble() + ", " + e.eval( "max" ).asDouble() + "]" );

		b.append( "s=(max-min)/" + bins + "; binBreaks = seq(min,max,by=s);\n" );
		b.append( "h0 = hist( sample0, breaks=c(min0,s+binBreaks[binBreaks>=min0 & binBreaks<=max0]), plot=F );\n" );
		b.append( "h1 = hist( sample1, breaks=c(min1,s+binBreaks[binBreaks>=min1 & binBreaks<=max1]), plot=F );\n\n" );

		if( plotOptions == null ) {
			plotOptions = "";
		} else {
			plotOptions = plotOptions.trim();
		}

		if( plotOptions.length() == 0 ) {
			// default
			plotOptions = ", xlim=c(min, max), ylim=c(0, max(h0$density,h1$density)), xlab=\"ratio of the scores\", ylab=\"density\", main=\"histogram for " + getClassifierName( cl )
							+ "\"";
		} else if( plotOptions.charAt( 0 ) != ',' ) {
			plotOptions = ", " + plotOptions;
		}

		// colors
		int c1 = 4, c2 = 2;
		b.append( "plot( h0, col=" + c1 + ", border=" + c1 + ", angle=45, density=" + density + ", freq=F " + plotOptions + ");\n" );
		b.append( "plot( h1, col=" + c2 + ", border=" + c2 + ", angle=-45, density=" + density + ", add=T, freq=F );\n" );

		//System.out.println(b);
		return b.toString();
	}

	/**
	 * This method returns an {@link ImageResult} containing a plot of the
	 * histograms of the scores. Further plotting options can be used to make
	 * such results comparable.
	 * 
	 * @param cl
	 *            the classifier
	 * @param class0
	 *            the sample for class 0
	 * @param class1
	 *            the sample for class 1
	 * @param e
	 *            the R environment
	 * @param bins
	 *            the number of bins in the plot
	 * @param density
	 *            the density of shading
	 * @param plotOptions
	 *            further plot options, e.g.
	 *            <code>xlim=c(-20,20), xlab=&quot;score&quot;</code>, if
	 *            <code>null</code> the default plotting options are chosen
	 * 
	 * @return the plot of the histograms
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see ImageResult#ImageResult(String, String,
	 *      java.awt.image.BufferedImage)
	 * @see REnvironment#plot(CharSequence, double, double)
	 */
	public static ImageResult plotScores( AbstractScoreBasedClassifier cl, DataSet class0, DataSet class1, REnvironment e, int bins,
			double density, String plotOptions ) throws Exception {
		return new ImageResult( "plot of the scores",
				"this plot shows the scores that are used to assign the classes",
				e.plot( getPlotScoresCmd( cl, class0, class1, e, bins, density, plotOptions ), 720, 360 ) );
	}

	/**
	 * This method creates a pdf containing a plot of the histograms of the
	 * scores. Further plotting options can be used to make such results
	 * comparable.
	 * 
	 * @param cl
	 *            the classifier
	 * @param class0
	 *            the sample for class 0
	 * @param class1
	 *            the sample for class 1
	 * @param e
	 *            the R environment
	 * @param bins
	 *            the number of bins in the plot
	 * @param density
	 *            the density of shading
	 * @param plotOptions
	 *            further plot options, e.g.
	 *            <code>xlim=c(-20,20), xlab=&quot;score&quot;</code>, if
	 *            <code>null</code> the default plotting options are chosen
	 * @param fName
	 *            the file name
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see REnvironment#plotToPDF(CharSequence, double, double, String, boolean)
	 */
	public static void plotScores( AbstractScoreBasedClassifier cl, DataSet class0, DataSet class1, REnvironment e, int bins, double density,
			String plotOptions, String fName ) throws Exception {
		e.plotToPDF( getPlotScoresCmd( cl, class0, class1, e, bins, density, plotOptions ), 10, 5, fName, true );
	}

	/**
	 * This method returns an {@link ImageResult} containing a scatter plot of
	 * the scores for the given classifiers <code>cl1</code> and
	 * <code>cl2</code>.
	 * 
	 * @param cl1
	 *            the first classifier
	 * @param cl2
	 *            the second classifier
	 * @param class0
	 *            the sample for class 0
	 * @param class1
	 *            the sample for class 1
	 * @param e
	 *            an R environment
	 * @param drawThreshold
	 *            indicates whether to draw the classification threshold or not
	 * 
	 * @return the scatter plot of the scores
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see ImageResult#ImageResult(String, String,
	 *      java.awt.image.BufferedImage)
	 * @see REnvironment#plot(CharSequence, double, double)
	 */
	public static ImageResult getScatterplot( AbstractScoreBasedClassifier cl1, AbstractScoreBasedClassifier cl2, DataSet class0,
			DataSet class1, REnvironment e, boolean drawThreshold ) throws Exception {
		if( cl1.getNumberOfClasses() != 2 || cl2.getNumberOfClasses() != 2 ) {
			throw new OperationNotSupportedException( "This method is only possible for 2-class-classifiers." );
		}

		e.createVector( "cl1_0", cl1.getScores( class0 ) );
		e.createVector( "cl2_0", cl2.getScores( class0 ) );
		e.createVector( "cl1_1", cl1.getScores( class1 ) );
		e.createVector( "cl2_1", cl2.getScores( class1 ) );

		e.voidEval( "xlim = c(min(cl1_0,cl1_1),max(cl1_0,cl1_1))" );
		e.voidEval( "ylim = c(min(cl2_0,cl2_1),max(cl2_0,cl2_1))" );

		String pltcmd = "plot( cl1_1, cl2_1, col=2, xlim = xlim, ylim = ylim, xlab=\"" + getClassifierName( cl1 )
						+ "\", ylab=\""
						+ getClassifierName( cl2 )
						+ "\", main=\""
						+ getTitle( class0, class1 )
						+ "\"); points( cl1_0, cl2_0, col=1, pch=2 );";
		if( drawThreshold ) {
			pltcmd += "lines( c(0,0), ylim, lty=2, col=3 ); lines( xlim, c(0,0), lty=2, col=3 )";
		}

		return new ImageResult( "scatterplot of the scores",
				"this plot shows the scores that are used to assign the classes scattered against each other",
				e.plot( pltcmd, 720, 720 ) );
	}
	
	/**
	 * Scatters the classification scores of two binary classifiers for given data.
	 * 
	 * @param cl1 the first classifier
	 * @param cl2 the second classifier
	 * @param e the {@link REnvironment}
	 * @param data the data
	 * 
	 * @return a scatter plot of the classification scores of two classifiers for given data
	 * 
	 * @throws Exception if the classifiers are not binary, the scores cannot be computed correctly or some problems occur in the communication with the {@link REnvironment}
	 */
	public static ImageResult getFancyScatterplot( AbstractScoreBasedClassifier cl1, AbstractScoreBasedClassifier cl2, REnvironment e, DataSet... data ) throws Exception {
		if( cl1.getNumberOfClasses() != 2 || cl2.getNumberOfClasses() != 2 ) {
			throw new OperationNotSupportedException( "This method is only possible for binary classifiers." );
		}
		
		e.voidEval( FANCY_SCATTTER );
		
		Sequence seq;
		double[] current;
		ArrayList<double[]> dlist = new ArrayList<double[]>();
		for( int d = 0; d < data.length; d++ ) {
			for( int n =  0; n < data[d].getNumberOfElements(); n++ ) {
				current = new double[3];
				seq = data[d].getElementAt( n );
				current[0] = cl1.getScore( seq, 0 ) - cl1.getScore( seq, 1 );
				current[1] = cl2.getScore( seq, 0 ) - cl2.getScore( seq, 1 );
				current[2] = d+1;
				dlist.add( current );
			}
		}
		
		e.createMatrix( "data", dlist.toArray( new double[0][0] ) );
		
		return new ImageResult( "scatterplot of the scores",
				"this plot shows the scores that are used to assign the classes scattered against each other",
				e.plot( "fancyScatter(data)", 720, 720 ) );
	}

	
	private static final String FANCY_SCATTTER =
			"fancyScatter <- function( data, breaks=20, density=10 ) {\n" +
				"m=matrix(c(2,1,4,3),ncol=2);\n" +
				"layout(m, widths = c(3,1), heights=c(1,3));\n" +
				"par(mar=c(3,3,1,1));\n" +
				"plot(data[,1],data[,2],col=data[,3],pch=16,xlab=\"\",ylab=\"\");\n" +
				
				"l = length(table(data[,3]));\n" +
				"degree = 180/l;\n" +
				"def.par <- par(no.readonly = TRUE); # save default, for resetting...\n" +
				"mar = par(\"mar\");\n\n" +
				
				"for( column in 1:2 ) {\n" +
					"myMar = mar;\n" +
					"myMar[column] = 0;\n" +
					"par(mar=myMar);\n" +
					"h = marginal(data,column=column,breaks=breaks);\n" +
					"bounds = c(min(data[,column]),max(data[,column]));\n\n" +
					
					"max=0;\n" +
					"for( i in 1:l ) {\n" +
						"max=max(max,h[[i]]$intensities);\n" +
					"}\n" +
					"for( i in 1:l ) {\n" +
						"if( column == 1 ) {\n" +
							"barplot(h[[i]]$intensities, add=i>1, axes=FALSE, ylim=c(0, max), space=0, col=i, angle=-90+i*degree, density=density);\n" +
						"} else {\n" +
							"barplot(h[[i]]$intensities, add=i>1, axes=FALSE, xlim=c(0, max), space=0, horiz=TRUE, col=i, angle=i*degree, density=density);\n" +
						"}\n" +
					"}\n" +
				"}\n" +
				"par( def.par );\n" +
			"}\n\n\n" 

			+
			
			"marginal <- function ( data, column=1, breaks=20 ) {\n" +
				"bounds = c(min(data[,column]),max(data[,column]));\n" +
				"t = table(data[,3]);\n" +
				"l = length(t);\n\n" +
				
				"h=list();\n" +
				"for( i in 1:l ) {\n" +
					"h[[i]] = hist( data[which(data[,3]==names(t)[i]),column], breaks=seq(bounds[1],bounds[2],length=breaks+1), plot=F );\n" +
				"}\n" +
				"return( h );\n" +
			"}\n";
	
	private static String getClassifierName( AbstractClassifier cl ) {
		CategoricalResult[] cat = cl.getClassifierAnnotation();
		String res = cat[0].getValue() + "(" + cat[1].getValue();
		for( int i = 2; i < cat.length; i++ ) {
			res += "; " + cat[i].getValue();
		}
		return res + ")";
	}

	private static String getTitle( DataSet class0, DataSet class1 ) {
		String res = "scatterplot for ";
		if( class0 != null && class1 != null ) {
			res += class0.getAnnotation() + " and " + class1.getAnnotation();
		} else {
			if( class0 == null ) {
				res += class1.getAnnotation();
			} else {
				res += class0.getAnnotation();
			}
		}
		return res;
	}
}
