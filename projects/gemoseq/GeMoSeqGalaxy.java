package projects.gemoseq;

import de.jstacs.tools.ui.galaxy.Galaxy;
import projects.gemoma.GeMoMaAnnotationFilter;

public class GeMoSeqGalaxy {

	public static void main(String[] args) throws Exception {
		
		Galaxy galaxy = new Galaxy("", true, true, new TranscriptPrediction(), new MergeGeMoMaGeMoSeq(), new GeMoMaAnnotationFilter());
		
		galaxy.run(args);

	}

}
