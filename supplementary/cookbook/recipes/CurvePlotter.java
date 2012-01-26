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
package supplementary.cookbook.recipes;

import java.io.File;

import de.jstacs.classifiers.AbstractClassifier;
import de.jstacs.classifiers.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.classifiers.performanceMeasures.AbstractPerformanceMeasure;
import de.jstacs.classifiers.performanceMeasures.PRCurve;
import de.jstacs.classifiers.performanceMeasures.PerformanceMeasureParameterSet;
import de.jstacs.classifiers.performanceMeasures.ROCCurve;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.results.ImageResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.REnvironment;

/**
 * This class exemplarily shows how to compute and plot curve for PR and ROC.
 * You need a running Rserve instance.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class CurvePlotter {

	/**
	 * @param args 
	 * <ul>
	 * <li>args[0] contains the path to the classifier</li>
	 * <li>args[1] contains the path to the foreground data set</li>
	 * <li>args[2] contains the path to the background data set</li>
	 * </ul>
	 */
	public static void main(String[] args) throws Exception {
		//read XML-representation from disk
		StringBuffer buf2 = FileManager.readFile( new File(args[0]+"myClassifier.xml") );
		 
		//create new classifier from read StringBuffer containing XML-code
		AbstractClassifier trainedClassifier = (AbstractClassifier) XMLParser.extractObjectForTags(buf2, "classifier");	

		//create a DataSet for each class from the input data, using the DNA alphabet
		DataSet[] test = new DataSet[2];
		test[0] = new DNADataSet( args[1] );
		
		//the length of our input sequences
		int length = test[0].getElementLength();

		test[1] = new DataSet( new DNADataSet( args[2] ), length );
		
		 
		AbstractPerformanceMeasure[] m = { new PRCurve(), new ROCCurve() };
		PerformanceMeasureParameterSet mp = new PerformanceMeasureParameterSet( m );
		ResultSet rs = trainedClassifier.evaluate( mp, true, test );
		 
		REnvironment r = null;
		try {
			r = new REnvironment();//might be adjusted
			for( int i = 0; i < rs.getNumberOfResults(); i++ )  {
				Result res = rs.getResultAt(i);
				if( res instanceof DoubleTableResult ) {
					DoubleTableResult dtr = (DoubleTableResult) res;
					ImageResult ir = DoubleTableResult.plot( r, dtr );
					REnvironment.showImage( dtr.getName(), ir.getValue() );
				} else {
					System.out.println( res );
				}
			}
		} catch( Exception e ) {
			e.printStackTrace();
		} finally {
			if( r != null ) {
				r.close();
			}
		}
	}
}