package projects.tals.epigenetic;

import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ui.cli.CLI;

public class EpiTALEcli {

	public static void main(String[] args) throws Exception {
		
		JstacsTool[] tools = new JstacsTool[] {
			new Bed2Bismark(),
			new BismarkConvertToPromotorSearch(),
			new BismarkMerge2Files(),
			new PileupConvertToPromotorSearch(),
			new NormalizePileupOutput(),
			new NarrowPeakConvertToPromotorSearch(),
			new QuickTBSPredictionToolEpigenetic()	
		};
		
		CLI cli = new CLI(tools);
		
		cli.run(args);

	}

}
