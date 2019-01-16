package projects.xanthogenomes;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ui.cli.CLI;
import projects.tals.prediction.QuickTBSPredictionTool;
import projects.tals.rnaseq.DerTALE;
import projects.xanthogenomes.tools.ClassAssignmentTool;
import projects.xanthogenomes.tools.ClassBuilderTool;
import projects.xanthogenomes.tools.ClassPresenceTool;
import projects.xanthogenomes.tools.LoadAndViewClassesTool;
import projects.xanthogenomes.tools.PredictAndIntersectTargetsTool;
import projects.xanthogenomes.tools.RenameTool;
import projects.xanthogenomes.tools.TALEAnalysisTool;
import projects.xanthogenomes.tools.TALEComparisonTool;
import projects.xanthogenomes.tools.TALEPredictionTool;



public class AnnoTALEcli {

	public static void main( String[] args ) throws Exception {
		
		JstacsTool[] tools = new JstacsTool[]{
				new TALEPredictionTool(),
				new TALEAnalysisTool(), 
				new ClassBuilderTool(), 
				new LoadAndViewClassesTool(), 
				new ClassAssignmentTool(), 
				new RenameTool(), 
				new PredictAndIntersectTargetsTool(),
                new ClassPresenceTool(),
                new TALEComparisonTool(),
                new QuickTBSPredictionTool(),
                new DerTALE()};

		CLI cli = new CLI( tools );
		
		cli.run( args );

	}

}
