package projects.xanthogenomes;

import de.jstacs.data.sequences.Sequence;


public class TALEConsensus {

	public static final Sequence repeat=getRepeat(), lastRepeat=getLastRepeat(), start=getStart(), end=getEnd();
	
	
	private static Sequence getRepeat(){
		try{
			return Sequence.create( Tools.Translator.DEFAULT.getProteinAlphabet(), "LTPDQVVAIASXXGGKQALETVQRLLPVLCQDHG" );
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	private static Sequence getLastRepeat(){
		try{
			return Sequence.create( Tools.Translator.DEFAULT.getProteinAlphabet(), "LTPDQVVAIASXXGGKQALE" );
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	private static Sequence getStart(){
		try{
			return Sequence.create( Tools.Translator.DEFAULT.getProteinAlphabet(), "MDPIRSRTPSPARELLPGPQPDRVQPTADRGGAPPAGGPLDGLPARRTMSRTRLPSPPAPSPAFSAGSFSDLLRQFDPSLLDTSLLDSMPAVGTPHTAAAPAEWDEVQSGLRAADDPPPTVRVAVTAARPPRAKPAPRRRAAQPSDASPAAQVDLRTLGYSQQQQEKIKPKVRSTVAQHHEALVGHGFTHAHIVALSQHPAALGTVAVTYQDIIRALPEATHEDIVGVGKQWSGARALEALLTEAGELRGPPLQLDTGQLLKIAKRGGVTAVEAVHAWRNALTGAPLN" );
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	private static Sequence getEnd(){
		try{
			return Sequence.create( Tools.Translator.DEFAULT.getProteinAlphabet(), "SIVAQLSRPDPALAALTNDHLVALACLGGRPALDAVKKGLPHAPELIRRINRRIPERTSHRVRVADLAHVVRVLGFFQSHSHPAQAFDDAMTQFGMSRHGLVQLFRRVGVTEFEARCGTLPPASQRWDRILQASGMKRAKPSPTSAQTPDQASLHAFADSLERDLDAPSPMHEGDQTRASSRKRSRSDRAVTGPSAQQSFEVRVPEQRDALHLPLSWRVKRPRTRIGGGLPDPGTPIAADLAASSTVMWEQDAAPFAGAADDFPAFNEEELAWLMELLPQSGSVGGTI" );
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	
}
