package projects.tals.epigenetic;

import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.results.Result;

public class MethylationSequenceAnnotation extends SequenceAnnotation{

	private Methylationprofil MP;
	
	public static final String TYPE = "methylationprofil"; 
	
	public MethylationSequenceAnnotation(String identifier,Methylationprofil MP, Result... results) {
		super(TYPE, identifier, results);
		this.MP=MP;
		
	}
	
	public Methylationprofil getMethylationprofile() {
		return MP;
	}

}
