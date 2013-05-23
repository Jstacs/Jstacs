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
import java.io.FileReader;

import de.jstacs.classifiers.performanceMeasures.AbstractPerformanceMeasure;
import de.jstacs.classifiers.performanceMeasures.AucPR;
import de.jstacs.classifiers.performanceMeasures.AucROC;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.ToolBox;

/**
 * This class contains a {@link #main(String[])}-method that allows to compute
 * some performance measures in a very simple way based on scores saved in files.
 * 
 * @author Jens Keilwagen
 */
public class AUC {

	private static void getValues( String fName, DoubleList score, DoubleList fg, DoubleList bg, double fgWeight ) throws Exception {
		BufferedReader r = new BufferedReader( new FileReader( fName ) );
		String line;
		while( ( line = r.readLine() ) != null ) {
			score.add( Double.parseDouble( line ) );
			fg.add( fgWeight );
			bg.add( 1-fgWeight );
		}
		r.close();
	}
	
	private static void getValues( String fName, DoubleList score, DoubleList fg, DoubleList bg ) throws Exception {
		BufferedReader r = new BufferedReader( new FileReader( fName ) );
		String line;
		String[] split;
		while( ( line = r.readLine() ) != null ) {
			split = line.split("\t");
			score.add( Double.parseDouble( split[0] ) );
			fg.add( Double.parseDouble( split[1] ) );
			bg.add( Double.parseDouble( split[2] ) );
		}
		r.close();
	}

	/**
	 * This {@link #main(String[])}-method computes some performance measures
	 * for two given score files. Note, the scores of class 0 should have the tendency to be higher
	 * than the scores of class 1. 
	 * 
	 * @throws Exception
	 *             if the performance measures could not be computed
	 * 
	 * @see PerformanceMeasure
	 * @see AucPR
	 * @see AucROC 
	 */
	public static void main( String[] args ) throws Exception {
		DoubleList score = new DoubleList();
		DoubleList fg = new DoubleList(), bg = new DoubleList();
		
		/*test case
		Random r = new Random();
		int n = 10;
		for( int i = 1; i < n; i++) {
			score.add( r.nextGaussian() );
			fg.add( i/(double)n );
			bg.add( 1d-i/(double)n );
		}*/
		
		//read
		switch( args.length ) {
			case 1:
				getValues(args[0], score, fg, bg);
				break;
			case 2:
				getValues(args[0], score, fg, bg, 1);
				getValues(args[1], score, fg, bg, 0);
				break;
			default:
				System.out.println("There are two options to start this program." );
				System.out.println();
				System.out.println("For unweighted data, please use:\njava -jar AUC.jar <fg> <bg>\nwhere <fg> and <bg> are files with one classification score per line.");
				System.out.println();
				System.out.println("For weighted data please use:\njava -jar AUC.jar <weighted>\nwhere <weighted> is a tab-delimited file with one classification score and the weights for fg and bg per line.");
				System.exit(0);
		}
		
		//sort
		double[] s = score.toArray();
		double[] fgWeights = fg.toArray();
		double[] bgWeights = bg.toArray();
		ToolBox.sortAlongWith(s, fgWeights, bgWeights);

		//compute
		AbstractPerformanceMeasure[] m = {
				new AucPR(),
				new AucROC()
		};
		
		ResultSet rs;
		for( AbstractPerformanceMeasure current : m ) {
			rs = current.compute( s, fgWeights, s, bgWeights );
			System.out.println( rs );
		}
	}

}
