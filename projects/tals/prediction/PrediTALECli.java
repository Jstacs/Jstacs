package projects.tals.prediction;

import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ui.cli.CLI;
import projects.tals.rnaseq.DerTALE;

public class PrediTALECli {

public static void main( String[] args ) throws Exception {
		
		JstacsTool[] tools = new JstacsTool[]{
                new QuickTBSPredictionTool(),
                new DerTALE()};

		CLI cli = new CLI( tools );
		
		cli.run( args );

	}
	
}
