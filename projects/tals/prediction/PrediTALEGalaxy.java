package projects.tals.prediction;

import de.jstacs.tools.ui.galaxy.Galaxy;

public class PrediTALEGalaxy {

	public static void main(String[] args) throws Exception {
		Galaxy gal = new Galaxy("", false, new QuickTBSPredictionTool());
		
		gal.run(args);

	}

}
