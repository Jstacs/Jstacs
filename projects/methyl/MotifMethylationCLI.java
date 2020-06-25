package projects.methyl;

import de.jstacs.tools.ui.cli.CLI;
import projects.encodedream.tools.MotifScores;
import projects.quickscan.QuickBindingSitePredictionTool;

public class MotifMethylationCLI {

	public static void main(String[] args) throws Exception{
		CLI cli = new CLI(new boolean[]{false,true,false,false,true,false,false},
				new ExtractMethylatedSequencesTool(),
				new MethylSlimDimontTool(),
				new MotifScanningTool(), 
				new EvaluateScoringTool(),
				new MotifScores(), 
				new QuickBindingSitePredictionTool(), 
				new MethylationSensitivity());
		cli.run(args);
		//System.out.println(cli.wikiPage("MeDeMo-1.0.jar"));
	}

}
