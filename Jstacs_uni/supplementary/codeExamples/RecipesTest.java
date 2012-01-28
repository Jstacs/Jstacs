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

import supplementary.cookbook.recipes.AlphabetCreation;
import supplementary.cookbook.recipes.CurvePlotter;
import supplementary.cookbook.recipes.DataLoader;
import supplementary.cookbook.recipes.GenDisMixClassifierTest;
import supplementary.cookbook.recipes.RserveTest;
import supplementary.cookbook.recipes.TrainSMBasedClassifierTest;

/**
 * Test class for <code>supplementary.cookbook.recipes</code>.
 * 
 * @author Jens Keilwagen
 */
public class RecipesTest {

	/**
	 * 
	 * @param args the arguments to be passed to the recipes;
	 * <li>args[0] contains the path to the foreground train data set</li>
	 * <li>args[1] contains the path to the background train data set</li>
	 * <li>args[2] contains the path to the foreground test data set</li>
	 * <li>args[3] contains the path to the background test data set</li>
	 * 
	 * @throws Exception if something went wrong
	 */
	public static void main(String[] args) throws Exception {
		
		AlphabetCreation.main(null);
		DataLoader.main( new String[]{"./supplementary/cookbook/recipes/"} );
		TrainSMBasedClassifierTest.main( new String[]{"./supplementary/cookbook/recipes/", args[0], args[1], args[2]} );
		GenDisMixClassifierTest.main( new String[]{args[0], args[1]} );
		RserveTest.main( new String[]{"./supplementary/cookbook/recipes/"} );
		CurvePlotter.main(  new String[]{"./supplementary/cookbook/recipes/", args[2], args[3]} );
	}
}