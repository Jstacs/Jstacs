package projects.gemoseq;

import de.jstacs.tools.ui.galaxy.Galaxy;

public class GeMoSeqGalaxy {

	public static void main(String[] args) throws Exception {
		
		Galaxy galaxy = new Galaxy("", true, true, new TranscriptPrediction(), new MergeGeMoMaGeMoSeq());
		
		galaxy.run(args);

	}

}
