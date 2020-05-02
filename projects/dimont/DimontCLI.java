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

package projects.dimont;

import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ui.cli.CLI;
import projects.quickscan.QuickBindingSitePredictionTool;
import projects.slim.SlimDimontTool;

public class DimontCLI {

	public static void main(String[] args) throws Exception {
		JstacsTool[] tools = new JstacsTool[] {
				new ExtractSequencesTool(),
				new DimontTool(),
				new SlimDimontTool(),
				new DimontPredictorTool(),
				new QuickBindingSitePredictionTool()
		};
		boolean[] threads = new boolean[] {
				false,
				true,
				true,
				false,
				false
		};
		CLI cli = new CLI(threads,tools);
		
		cli.run(args);

	}

}
