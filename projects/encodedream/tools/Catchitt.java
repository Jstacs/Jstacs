package projects.encodedream.tools;
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
import de.jstacs.tools.ui.cli.CLI;

public class Catchitt {

	public static void main(String[] args) throws Exception {
		
		CLI cli = new CLI(new boolean[]{false,false,true,false,true,false}, new ChromatinAccessibility(), new MethylationLevels(), new MotifScores(), new DeriveLabelTool(), new IterativeTrainingTool(), new PredictionTool());
		//cli.wiki();
		cli.run(args);
		
	}

}
