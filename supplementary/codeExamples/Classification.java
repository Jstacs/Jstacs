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

package supplementary.codeExamples;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.classifier.ScoreBasedPerformanceMeasureDefinitions;
import de.jstacs.classifier.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.classifier.MeasureParameters.Measure;
import de.jstacs.classifier.ScoreBasedPerformanceMeasureDefinitions.ThresholdMeasurePair;
import de.jstacs.io.FileManager;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.REnvironment;

/**
 * This class contains a {@link #main(String[])}-method that allows to compute
 * some performance measures in a very simple way based on scores saved in files.
 * 
 * @author Jens Keilwagen
 * 
 * @see Classification#main(String[])
 */
public class Classification {

	private static double[] getSortedValues( String fName ) throws Exception {
		BufferedReader r = new BufferedReader( new FileReader( fName ) );
		String line;
		DoubleList dlist = new DoubleList();
		while( ( line = r.readLine() ) != null ) {
			dlist.add( Double.parseDouble( line ) );
		}
		r.close();
		double[] res = dlist.toArray();
		Arrays.sort( res );
		return res;
	}

	/**
	 * This {@link #main(String[])}-method computes some performance measures
	 * for two given score files. Note, the scores of class 0 should have the tendency to be higher
	 * than the scores of class 1. 
	 * 
	 * @param args
	 *            the name/path of the score files
	 *            <ul>
	 *            <li> <code>args[0]</code> file name of score file for class 0
	 *            <li> <code>args[1]</code> file name of score file for class 1
	 *            <li> optionally <code>args[2]</code> rserve server name
	 *            <li> optionally <code>args[3]</code> user name for rserve (only together with login)
	 *            <li> optionally <code>args[4]</code> login for rserve (only together with user name)
	 *            </ul>
	 * 
	 * @throws Exception
	 *             if the performance measures could not be computed
	 * 
	 * @see ScoreBasedPerformanceMeasureDefinitions#getSensitivityForSpecificity(double[], double[], double)
	 * @see ScoreBasedPerformanceMeasureDefinitions#getFPRForSensitivity(double[], double[], double)
	 * @see ScoreBasedPerformanceMeasureDefinitions#getPPVForSensitivity(double[], double[], double)
	 * @see ScoreBasedPerformanceMeasureDefinitions#getAUC_ROC(double[], double[], java.util.AbstractList)
	 * @see ScoreBasedPerformanceMeasureDefinitions#getAUC_PR(double[], double[], java.util.AbstractList)
	 */
	public static void main( String[] args ) throws Exception {
		double[] fg = getSortedValues( args[0] );
		double[] bg = getSortedValues( args[1] );
		
		boolean plot = args.length > 2;
		REnvironment r = null;
		LinkedList<double[]> list = null;
		DoubleTableResult dtr;
		try{
			if( plot ) {
				list = new LinkedList<double[]>();
				switch( args.length ) {
					case 3:
						//server
						r = new REnvironment( args[2], null, null );
						break;
					case 5:
						//server,user,login
						r = new REnvironment( args[2], args[3], args[4] );
						break;
					default:
						System.out.println( "You have to specify the server and optionally username and login." );
				}
			}
	
			ThresholdMeasurePair tmp = ScoreBasedPerformanceMeasureDefinitions.getSensitivityForSpecificity( fg, bg, 0.999 );
			System.out.println( "Sn (with Sp=99.9%) =\t" + tmp.getMeasure() + "\tthreshold =\t" + tmp.getThreshold() );
			tmp = ScoreBasedPerformanceMeasureDefinitions.getFPRForSensitivity( fg, bg, 0.95 );
			System.out.println( "FPR (with Sn=95%) =\t" + tmp.getMeasure() + "\tthreshold =\t" + tmp.getThreshold() );
			tmp = ScoreBasedPerformanceMeasureDefinitions.getPPVForSensitivity( fg, bg, 0.95 );
			System.out.println( "PPV with Sn=95%\t" + tmp.getMeasure() + "\tthreshold =\t" + tmp.getThreshold() );
			System.out.println( "AUC_ROC\t" + ScoreBasedPerformanceMeasureDefinitions.getAUC_ROC( fg, bg, list ) );
			if( plot ) {
				dtr = new DoubleTableResult("ROC","",list);
				list.clear();
				r.plotToPDF( DoubleTableResult.getPlotCommands( r, Measure.ReceiverOperatingCharacteristicCurve.getNameString(), dtr ).toString(), "./ROC-curve.pdf", true );
				FileManager.writeFile( new File("ROC.xml"), dtr.toXML() );
			}
			System.out.println( "AUC_PR\t" + ScoreBasedPerformanceMeasureDefinitions.getAUC_PR( fg, bg, list ) );
			if( plot ) {
				dtr = new DoubleTableResult("PR","",list);
				list.clear();
				r.plotToPDF( DoubleTableResult.getPlotCommands( r, Measure.PrecisionRecallCurve.getNameString(), dtr ).toString(), "./PR-curve.pdf", true );
				FileManager.writeFile( new File("PRC.xml"), dtr.toXML() );
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
