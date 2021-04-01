package projects.tals.epigenetic;

import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ui.cli.CLI;

public class EpiTALEcli {

	public static void main(String[] args) throws Exception {
		
		JstacsTool[] tools = new JstacsTool[] {
			new Bed2Bismark(),
			new BismarkMerge2Files(),
			new BismarkConvertToPromoterSearch(),
			new PileupTool(),
			new NormalizePileupOutput(),
			new PileupConvertToPromoterSearch(),
			new NarrowPeakConvertToPromoterSearch(),
			new QuickTBSPredictionToolEpigenetic()
		};
		
		CLI cli = new CLI(tools);
		
		cli.run(args);

	}

}
