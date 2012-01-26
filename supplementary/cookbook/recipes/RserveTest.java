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

import java.awt.image.BufferedImage;

import javax.swing.JFrame;

import org.rosuda.REngine.REXP;

import de.jstacs.utils.REnvironment;

/**
 * This class exemplarily shows how to use R from within Jstacs.
 * You need a running Rserve instance.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class RserveTest {

	/**
	 * @param args args[0] contains the path that is used to store plots
	 */
	public static void main(String[] args) throws Exception {
		REnvironment e = null;
		try {
			//create a connection to R with YOUR server name, login and password
			e = new REnvironment();//might be adjusted
		 
			System.out.println( "java: " + System.getProperty( "java.version" ) );
			System.out.println();
			System.out.println( e.getVersionInformation() );
		 
			// compute something in R
			REXP erg = e.eval( "sin(10)" );
			System.out.println( erg.asDouble() );
		 
			//create a histrgram in R in 3 steps
			//1) create the data
			e.voidEval( "a = 100;" );
			e.voidEval( "n = rnorm(a)" );
			//2) create the plot command
			String plotCmd = "hist(n,breaks=a/5)";
			//3a) plot as pdf
			e.plotToPDF( plotCmd, args[0]+"/test.pdf", true );
			//or
			//3b) create an image and show it
			BufferedImage i = e.plot( plotCmd, 640, 480 );
			REnvironment.showImage( "histogramm", i, JFrame.EXIT_ON_CLOSE );
		 
		} catch ( Exception ex ) {
			ex.printStackTrace();
		} finally {
			if( e != null ) {
				try {
					//close REnvironment correctly
					e.close();
				} catch ( Exception e1 ) {
					System.err.println( "could not close REnvironment." );
					e1.printStackTrace();
				}
			}
		}
	}
}