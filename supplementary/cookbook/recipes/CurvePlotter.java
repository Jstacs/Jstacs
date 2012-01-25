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
			for( AbstractPerformanceMeasure s : m )  {
				DoubleTableResult dtr = (DoubleTableResult) rs.getResultAt( 1 );
				ImageResult ir = DoubleTableResult.plot( r, dtr );
				REnvironment.showImage( dtr.getName(), ir.getValue() );
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