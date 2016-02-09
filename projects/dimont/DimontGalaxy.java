package projects.dimont;

import de.jstacs.tools.ui.galaxy.Galaxy;

public class DimontGalaxy {

	public static void main(String[] args) throws Exception {
		
		DimontTool tool = new DimontTool();
		
		DimontPredictorTool pred = new DimontPredictorTool();
		
		ExtractSequencesTool est = new ExtractSequencesTool();
		
		Galaxy gal = new Galaxy("", true, est,tool,pred);

		gal.run(args);
		
	}

}
