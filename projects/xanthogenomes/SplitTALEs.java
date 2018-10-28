package projects.xanthogenomes;

import java.io.FileWriter;
import java.io.PrintWriter;

import de.jstacs.algorithms.alignment.Alignment;
import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.PairwiseStringAlignment;
import de.jstacs.algorithms.alignment.cost.AffineCosts;
import de.jstacs.algorithms.alignment.cost.Costs;
import de.jstacs.algorithms.alignment.cost.MatrixCosts;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.tools.Protocol;
import projects.xanthogenomes.TALE.Repeat;
import projects.xanthogenomes.Tools.Aligner;
import projects.xanthogenomes.Tools.ProteinAlphabetContainer;
import projects.xanthogenomes.Tools.Translator;


public class SplitTALEs {

	
	public static TALE[] split(String id, Sequence tale, Protocol protocol) throws Exception {
		Costs costs = new AffineCosts( 3, new MatrixCosts( Aligner.DEFAULT.getMatrix(), 2 ) ); 
		Alignment align = new Alignment( costs );
		
		AlphabetContainer protAlph = ProteinAlphabetContainer.SINGLETON;
		Sequence prot = null;
		if(tale.getAlphabetContainer().checkConsistency( protAlph )){
			prot = tale;
		}else{
			prot = Tools.Translator.DEFAULT.translate( tale, 0 );
		}
		
		int talelen = prot.getLength();
		
		
		Sequence repeat = TALEConsensus.repeat;
		Sequence lastRepeat = TALEConsensus.lastRepeat;
		Sequence start = TALEConsensus.start;
		Sequence end = TALEConsensus.end;
		
		
		int nreps = 0;
		
		StringBuffer first = new StringBuffer();
		first.append(start);
		/*for(int j=0;j<nreps;j++){
			first.append( repeat );
		}*/
		String last = lastRepeat.toString()+end.toString();
		
		Sequence template = Sequence.create(protAlph, first.toString()+last);
		PairwiseStringAlignment bestAl = align.getAlignment( AlignmentType.GLOBAL, template, prot );
		Double best = bestAl.getCost();
		
		while(true){
			
			template = Sequence.create(protAlph, first.toString()+repeat.toString()+last);
			
			PairwiseStringAlignment psa = align.getAlignment( AlignmentType.GLOBAL, template, prot );
			if(psa.getCost() > best){
				break;
			}else{
				best = psa.getCost();
				bestAl = psa;
				nreps++;
				first.append( repeat );
			}
			
		}
		
		//System.out.println(bestAl);
		
		String al1 = bestAl.getAlignedString( 0 );
		String al2 = bestAl.getAlignedString( 1 );
		int alignmentIndex = 0;
		int index = 0;
		while( index<start.getLength() ){
			if(al1.charAt( alignmentIndex ) != '-'){
				index++;
			}
			alignmentIndex++;
		}
		
		Sequence foundStart = prot.getSubSequence( 0, al2.substring( 0, alignmentIndex ).replaceAll( "-", "" ).length() );

		int off = alignmentIndex;
		int origOff = foundStart.getLength();
		
		Repeat[] foundRepeats = new Repeat[nreps+1];
		for(int j=0;j<nreps;j++){
			index = 0;
			while(index < repeat.getLength()){
				if(al1.charAt( alignmentIndex ) != '-'){
					index++;
				}
				alignmentIndex++;
			}
			
			int currL = al2.substring( off, alignmentIndex ).replaceAll( "-", "" ).length();
			foundRepeats[j] = new Repeat( prot.getSubSequence( origOff, currL) );
			origOff += currL;
			off = alignmentIndex;
		}
		
		index = 0;
		while(index<lastRepeat.getLength()){
			if(al1.charAt( alignmentIndex ) != '-'){
				index++;
			}
			alignmentIndex++;
		}
		
		int currL = al2.substring( off, alignmentIndex ).replaceAll( "-", "" ).length();
		if(currL == 0){
			Repeat[] reps2 = new Repeat[nreps];
			System.arraycopy( foundRepeats, 0, reps2, 0, reps2.length );
			foundRepeats = reps2;
		}else{
			foundRepeats[foundRepeats.length-1] = new Repeat( prot.getSubSequence( origOff, currL ) );
		}
		off = alignmentIndex;
		origOff += currL;
		
		index = 0;
		while(index<end.getLength()){
			if(al1.charAt( alignmentIndex ) != '-'){
				index++;
			}
			alignmentIndex++;
		}
		currL = al2.substring( off, alignmentIndex ).replaceAll( "-", "" ).length();
		Sequence foundEnd = prot.getSubSequence( origOff, currL );
		
		//TALE protTALE = new TALE( id, foundStart, foundRepeats, foundEnd, true );
		TALE dnaTALE = null;
		
		if(tale.getAlphabetContainer().checkConsistency( DNAAlphabetContainer.SINGLETON )){
			Sequence dnaStart = tale.getSubSequence( 0, foundStart.getLength()*3 );
			off = foundStart.getLength();
			Repeat[] dnaRepeat = new Repeat[foundRepeats.length];
			for(int j=0;j<dnaRepeat.length;j++){
				dnaRepeat[j] = new Repeat( tale.getSubSequence( off*3, foundRepeats[j].getRepeat().getLength()*3 ) );
				off += foundRepeats[j].getRepeat().getLength();
			}
			//Sequence dnaLastRepeat = dna.getSubSequence( off*3, (foundLastRepeat.getLength())*3 );
			//off += foundLastRepeat.getLength();
			Sequence dnaEnd = tale.getSubSequence( off*3, (foundEnd.getLength())*3 );



			dnaTALE = new TALE( id, dnaStart, dnaRepeat, dnaEnd, true );
		}
		
		TALE transl = null;
		if(dnaTALE != null){
			try{
				transl = dnaTALE.getTranslatedTALE( Translator.DEFAULT );
			}catch(Exception e){
				protocol.appendWarning("Translated version of TALE "+id+" could not be created. Reason\n"+e.getMessage()+"\n");
			}
		}else{
			transl = new TALE( id, foundStart, foundRepeats, foundEnd, true );
		}
		return new TALE[]{dnaTALE, transl };
	}
	
	public static void main( String[] args ) throws Exception {
		
		AlphabetContainer protAlph = Tools.Translator.DEFAULT.getProteinAlphabet();
		
		Sequence repeat = TALEConsensus.repeat;
		Sequence lastRepeat = TALEConsensus.lastRepeat;
		Sequence start = TALEConsensus.start;
		Sequence end = TALEConsensus.end;
		
		DNADataSet ds = new DNADataSet( args[0],'>', new SimpleSequenceAnnotationParser() );
		
		PrintWriter starts = new PrintWriter( new FileWriter( args[0]+"_starts.txt" ) );
		PrintWriter repeats = new PrintWriter( new FileWriter( args[0]+"_repeats.txt" ) );
		PrintWriter ends = new PrintWriter( new FileWriter( args[0]+"_ends.txt" ) );
		PrintWriter dnaStarts = new PrintWriter( new FileWriter( args[0]+"_dnastarts.txt" ) );
		PrintWriter dnaRepeats = new PrintWriter( new FileWriter( args[0]+"_dnarepeats.txt" ) );
		PrintWriter dnaEnds = new PrintWriter( new FileWriter( args[0]+"_dnaends.txt" ) );
		
		Costs costs = new AffineCosts( 3, new MatrixCosts( Aligner.DEFAULT.getMatrix(), 2 ) ); 
		Alignment align = new Alignment( costs );
		
		for(int i=0;i<ds.getNumberOfElements();i++){
			Sequence dna = ds.getElementAt( i );
			Sequence prot = Tools.Translator.DEFAULT.translate( dna, 0 );
			
			int talelen = prot.getLength();
			//talelen -= start.getLength()+end.getLength()+lastRepeat.getLength();
			//int nreps = (int)Math.floor(talelen/(double)repeat.getLength())-1;
			int nreps = 1;
			
			StringBuffer first = new StringBuffer();
			first.append(start);
			for(int j=0;j<nreps;j++){
				first.append( repeat );
			}
			String last = lastRepeat.toString()+end.toString();
			
			Double best = Double.POSITIVE_INFINITY;
			PairwiseStringAlignment bestAl = null;
			while(true){
				
				Sequence template = Sequence.create(protAlph, first.toString()+repeat.toString()+last);
				
				PairwiseStringAlignment psa = align.getAlignment( AlignmentType.GLOBAL, template, prot );
				if(psa.getCost() > best){
					break;
				}else{
					best = psa.getCost();
					bestAl = psa;
					nreps++;
					first.append( repeat );
				}
				
			}
			
			System.out.println(bestAl);
			
			String al1 = bestAl.getAlignedString( 0 );
			String al2 = bestAl.getAlignedString( 1 );
			int alignmentIndex = 0;
			int index = 0;
			while( index<start.getLength() ){
				if(al1.charAt( alignmentIndex ) != '-'){
					index++;
				}
				alignmentIndex++;
			}
			
			String foundStart = al2.substring( 0, alignmentIndex ).replaceAll( "-", "" );

			int off = alignmentIndex;
			
			String[] foundRepeats = new String[nreps];
			for(int j=0;j<nreps;j++){
				index = 0;
				while(index < repeat.getLength()){
					if(al1.charAt( alignmentIndex ) != '-'){
						index++;
					}
					alignmentIndex++;
				}
				
				foundRepeats[j] = al2.substring( off, alignmentIndex ).replaceAll( "-", "" );
				off = alignmentIndex;
			}
			
			index = 0;
			while(index<lastRepeat.getLength()){
				if(al1.charAt( alignmentIndex ) != '-'){
					index++;
				}
				alignmentIndex++;
			}
			
			String foundLastRepeat = al2.substring( off, alignmentIndex ).replaceAll( "-", "" );
			off = alignmentIndex;
			
			index = 0;
			while(index<end.getLength()){
				if(al1.charAt( alignmentIndex ) != '-'){
					index++;
				}
				alignmentIndex++;
			}
			String foundEnd = al2.substring( off, alignmentIndex ).replaceAll( "-", "" );
			
			
			String dnaStart = dna.toString().substring( 0, foundStart.length()*3 );
			off = foundStart.length();
			String[] dnaRepeat = new String[foundRepeats.length];
			for(int j=0;j<dnaRepeat.length;j++){
				dnaRepeat[j] = dna.toString().substring( off*3, (off+foundRepeats[j].length())*3 );
				off += foundRepeats[j].length();
			}
			String dnaLastRepeat = dna.toString().substring( off*3, (off+foundLastRepeat.length())*3 );
			off += foundLastRepeat.length();
			String dnaEnd = dna.toString().substring( off*3, (off+foundEnd.length())*3 );
			
			
			String id = (String)dna.getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultAt( 0 ).getValue();
			
			starts.println(">"+id+"\n"+foundStart);
			dnaStarts.println(">"+id+"\n"+dnaStart);
			repeats.println(">"+id);
			dnaRepeats.println(">"+id);
			for(int j=0;j<foundRepeats.length;j++){
				repeats.println(foundRepeats[j]);
				dnaRepeats.println(dnaRepeat[j]);
			}
			repeats.println(foundLastRepeat);
			dnaRepeats.println(dnaLastRepeat);
			ends.println(">"+id+"\n"+foundEnd);
			dnaEnds.println(">"+id+"\n"+dnaEnd);
		}

		starts.close();
		dnaStarts.close();
		repeats.close();
		dnaRepeats.close();
		ends.close();
		dnaEnds.close();
		
	}
	

}
