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

import de.jstacs.classifier.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.classifier.performanceMeasures.AbstractPerformanceMeasure;
import de.jstacs.classifier.performanceMeasures.PerformanceMeasureParameterSet;
import de.jstacs.io.FileManager;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.Measure;
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
	 * @see PerformanceMeasureParameterSet
	 */
	public static void main( String[] args ) throws Exception {
		double[] fg = getSortedValues( args[0] );
		double[] bg = getSortedValues( args[1] );
		
		boolean plot = args.length > 2;
		REnvironment r = null;
		DoubleTableResult dtr;
		try{
			if( plot ) {
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
			
			PerformanceMeasureParameterSet performance = PerformanceMeasureParameterSet.createFilledParameters( plot, 0.999, 0.95, 0.95, 1 );
			boolean isNumeric = true;
			AbstractPerformanceMeasure[] m = performance.getAllMeasures();
			ResultSet rs;
			Result res;
			for( AbstractPerformanceMeasure current : m ) {
				rs = current.compute( fg, bg );
				if( rs != null ) {
					System.out.println( rs );
					if( plot ) {
						for( int i = 0; i < rs.getNumberOfResults(); i++ ) {
							res = rs.getResultAt(i);
							if( res instanceof DoubleTableResult ) {
								dtr = (DoubleTableResult) res;
								String name = dtr.getName().replaceAll( " ", "-" );
								r.plotToPDF( DoubleTableResult.getPlotCommands( r, dtr.getName(), dtr ).toString(), "./"+name+".pdf", true );
								FileManager.writeFile( new File("./"+name+".xml"), dtr.toXML() );
							}
						}
					}
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
