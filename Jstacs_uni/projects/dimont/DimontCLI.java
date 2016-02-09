package projects.dimont;

import de.jstacs.tools.ui.cli.CLI;

public class DimontCLI {

	public static void main(String[] args) throws Exception {
		
		DimontTool tool = new DimontTool();
		
		DimontPredictorTool pred = new DimontPredictorTool();
		
		ExtractSequencesTool est = new ExtractSequencesTool();
		
		CLI cli = new CLI(new boolean[]{false,true,false},est,tool,pred);
		
		cli.run(args);

	}

}
