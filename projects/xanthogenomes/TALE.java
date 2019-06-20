package projects.xanthogenomes;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.naming.OperationNotSupportedException;

import de.jstacs.Singleton;
import de.jstacs.Storable;
import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.PairwiseStringAlignment;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import projects.xanthogenomes.Tools.Aligner;
import projects.xanthogenomes.Tools.Translator;


public class TALE implements Storable {
	
	public enum Type{
		NORMAL,
		LONG,
		SHORT,
		THIRTYFIVE,
		UNKNOWN
	}
	
	public static class Repeat implements Storable{
		
		private Sequence repeat;
		private String rvd;
		private Sequence[] maskedRepeats;
		private int rvdPosition;
		
		
		public Repeat(StringBuffer xml) throws NonParsableException, IllegalArgumentException, WrongAlphabetException{
			xml = XMLParser.extractForTag( xml, "Repeat" );
			repeat = XMLParser.extractSequencesWithTags( xml, "repeat" )[0];
			rvd = (String)XMLParser.extractObjectForTags( xml, "rvd" );
			maskedRepeats = XMLParser.extractSequencesWithTags( xml, "maskedRepeats" );
			rvdPosition = (Integer)XMLParser.extractObjectForTags( xml, "rvdPosition" );
		}
		
		public StringBuffer toXML(){
			StringBuffer xml = new StringBuffer();
			XMLParser.appendSequencesWithTags( xml, "repeat", repeat );
			XMLParser.appendObjectWithTags( xml, rvd, "rvd" );
			XMLParser.appendSequencesWithTags( xml, "maskedRepeats", maskedRepeats );
			XMLParser.appendObjectWithTags( xml, rvdPosition, "rvdPosition" );
			XMLParser.addTags( xml, "Repeat" );
			return xml;
		}
		
		public int getRvdPosition() {
			return rvdPosition;
		}


		
		public int getRvdLength() {
			if(rvd != null){
				return rvd.charAt(1)=='*' ? 1 : 2;
			}else{
				return 0;
			}
		}

		private int rvdLength;
		
		public Repeat(Sequence repeat){
			this.repeat = repeat;
		}

		public Type getType(){
			if(this.repeat != null && this.repeat.getAlphabetContainer().checkConsistency(Translator.DEFAULT.getProteinAlphabet())){
				int expected = rvd.endsWith("*") ? 33 : 34;
				int diff = repeat.getLength()-expected;
				if(diff == 0){
					return Type.NORMAL;
				}else if(diff < 0){
					return Type.SHORT;
				}else if(diff == 1){
					return Type.THIRTYFIVE;
				}else{
					return Type.LONG;
				}
			}else{
				return Type.UNKNOWN;
			}
		}
		
		public Sequence getRepeat() {
			return repeat;
		}

		
		public String getRvd() {
			return rvd;
		}

		
		public Sequence[] getMaskedRepeats() {
			return maskedRepeats;
		}
		
		
		
		private int[] getRVDIndexAndLength(Sequence consensus) throws OperationNotSupportedException{
			PairwiseStringAlignment al = Aligner.DEFAULT.align( repeat.reverse(), consensus.reverse(), AlignmentType.GLOBAL );
			//System.out.println(al);
			String ca = new StringBuffer(al.getAlignedString( 1 )).reverse().toString();
			int k=0,j=0;
			while( k<11 ){
				if(ca.charAt( j ) != '-'){
					k++;
				}
				j++;
			}
			int l=j;
			while( k<12 ){
				if(ca.charAt( l ) != '-'){
					k++;
				}
				l++;
			}
			
			String ra = new StringBuffer(al.getAlignedString( 0 )).reverse().toString();
			//System.out.println(ca+"\n"+ra+"\n");
			
			char twelve = ra.charAt( j );
			char thirteen = ra.charAt( l );
			
			int idx = ra.substring( 0, j+1 ).replaceAll( "-", "" ).length()-1;
			
			int length = 2;
			if(twelve == '-' && thirteen == '-'){
				throw new RuntimeException("Both RVD positions are gaps.\n"+repeat+"\n"+consensus+"\n");
			}else if(twelve == '-'){
				idx++;
				length = 1;
			}else if(thirteen == '-'){
				length = 1;
			}
			
			return new int[]{idx,length};
		}
		
		private Sequence[] getBlankedRepeats(int[] idl) throws WrongAlphabetException, WrongSequenceTypeException{
			int index = idl[0];
			int length = idl[1];
			
			if(length == 1){

				int[] cont1 = new int[index];
				int[] cont2 = new int[repeat.getLength()-index-1];
				int j=0;
				for(int i=0;i<repeat.getLength();i++){
					if( i<index ){
						cont1[j] = repeat.discreteVal( i );
						j++;
					}else if(i>index){
						cont2[j] = repeat.discreteVal( i );
						j++;
					}else{
						j=0;
					}
				}
				
				return new Sequence[]{new IntSequence( repeat.getAlphabetContainer(), cont1 ),new IntSequence( repeat.getAlphabetContainer(), cont2 )};
			}else{
				int[] cont1 = new int[index];
				int[] cont2 = new int[repeat.getLength()-index-2];
				int j=0;
				for(int i=0;i<repeat.getLength();i++){
					if( i<index){
						cont1[j] = repeat.discreteVal( i );
						j++;
					}else if( i>index+length-1 ){
						cont2[j] = repeat.discreteVal( i );
						j++;
					}else{
						j=0;
					}
				}
				
				return new Sequence[]{new IntSequence( repeat.getAlphabetContainer(), cont1 ),new IntSequence( repeat.getAlphabetContainer(), cont2 )};
			}
		}
		
		public Sequence[] getBlankedRepeats(Sequence consensus) throws WrongAlphabetException, WrongSequenceTypeException, OperationNotSupportedException{
			//Sequence repeat = repeats[i].getRepeat();

			
			int[] idl = getRVDIndexAndLength( consensus );
			
			return getBlankedRepeats( idl );
			
		}
		
		//idl: 0: index, 1: length
		private String getRVD(int[] idl){
			int index = idl[0];
			int length = idl[1];
			
			StringBuffer reps = new StringBuffer();
			
			if(length == 1){
				reps.append( repeat.getSubSequence( index, 1 )+"*" );
			}else{
				reps.append( repeat.getSubSequence( index, 2 ).toString() );
			}
			//System.out.println(reps);
			return reps.toString();
		}
		
		public String getRVD(Sequence consensus) throws OperationNotSupportedException{
			//Sequence repeat = repeats[i].getRepeat();
			int[] idl = getRVDIndexAndLength( consensus );
			
			return getRVD( idl );
		}
		
		
		
	}
	
	
	
	private String strain;
	private String accession;
	private Integer startPos;
	private Integer endPos;
	private Boolean strand;
	
	void addAnnotation(String strain, String accession, Integer startPos, Integer endPos, Boolean strand){
		this.strain = strain;
		if(dnaOriginal != null){
			dnaOriginal.strain = strain;
		}
		this.accession = accession;
		if(dnaOriginal != null){
			dnaOriginal.accession = accession;
		}
		this.startPos = startPos;
		if(dnaOriginal != null){
			dnaOriginal.startPos = startPos;
		}
		this.endPos = endPos;
		if(dnaOriginal != null){
			dnaOriginal.endPos = endPos;
		}
		this.strand = strand;
		if(dnaOriginal != null){
			dnaOriginal.strand = strand;
		}
	}
	
	public void setStrain(String strain){
		this.strain = strain;
		if(dnaOriginal != null){
			dnaOriginal.strain = strain;
		}
	}
	
	private String annotationToColumns(){
		StringBuffer buf = new StringBuffer();
		if(strain != null){
			buf.append(strain);
		}
		buf.append("\t");
		if(accession != null){
			buf.append(accession);
		}
		buf.append("\t");
		if(startPos != null){
			buf.append(startPos);
		}
		buf.append("\t");
		if(endPos != null){
			buf.append(endPos);
		}
		buf.append("\t");
		if(strand != null){
			buf.append((strand ? "+1" : "-1"));
		}
		return buf.toString();
	}
	
	public String annotationToString(){
		StringBuffer buf = new StringBuffer();
		/*if(strain != null){
			buf.append(strain);
		}
		if(accession != null){
			buf.append(" ("+accession+")");
		}*/
		if(startPos != null && endPos != null){
			buf.append("["+(accession != null ? accession+": " : "")+startPos+"-"+endPos+":"+(strand != null ? (strand ? "+1" : "-1") : "")+"]");
		}
		return buf.toString();
	}
	
	private String id;
	private Sequence start;


	private Repeat[] repeats;
	private Sequence end;
	
	private Sequence rvdSequence;
	
	private boolean isNew;
	
	public Sequence getRvdSequence() {
		return rvdSequence;
	}

	private TALE dnaOriginal;
	
	
	public TALE getDnaOriginal() {
		return dnaOriginal;
	}


	
	
	
	private void setDnaOriginal( TALE dnaOriginal ) {
		this.dnaOriginal = dnaOriginal;
	}

	
	public String toString(){
		//return ">"+id+"\n"+getRvdSequence();
		return id;
	}
	
	public TALE( String id, Sequence rvds, boolean parseId, boolean isNew ) throws IllegalArgumentException, WrongAlphabetException{
		this.id = parseId ? parse(id) : id;
		this.start = this.end = null;
		this.isNew = isNew;
		if(!RVDAlphabetContainer.SINGLETON.checkConsistency(rvds.getAlphabetContainer())){
			throw new WrongAlphabetException();
		}
		this.rvdSequence = rvds;
		this.repeats = new Repeat[rvdSequence.getLength()];
		for(int i=0;i<rvdSequence.getLength();i++){
			this.repeats[i] = new Repeat((Sequence)null);
			this.repeats[i].rvdPosition = -1;
			this.repeats[i].rvd = RVDAlphabetContainer.SINGLETON.getSymbol(i,rvdSequence.discreteVal(i) );
		}
		
	}
	
	public TALE( String id, Sequence start, Repeat[] repeats, Sequence end) throws OperationNotSupportedException, WrongAlphabetException, WrongSequenceTypeException{
		this(true, id, start, repeats, end);
	}
	
	public TALE( boolean parseId, String id, Sequence start, Repeat[] repeats, Sequence end) throws OperationNotSupportedException, WrongAlphabetException, WrongSequenceTypeException{
		this( parseId, id, start, repeats, end, false );
	}
	
	public TALE( String id, Sequence start, Repeat[] repeats, Sequence end, boolean isNew ) throws OperationNotSupportedException, WrongAlphabetException, WrongSequenceTypeException{
		this(true, id, start, repeats, end, isNew);
	}
	
	public TALE( boolean parseId, String id, Sequence start, Repeat[] repeats, Sequence end, boolean isNew ) throws WrongAlphabetException, WrongSequenceTypeException, OperationNotSupportedException{
		this.id = parseId ? parse(id) : id;
		this.start = start;
		this.repeats = repeats;
		this.end = end;
		this.isNew = isNew;
		
		if(start.getAlphabetContainer().checkConsistency( Translator.DEFAULT.getProteinAlphabet() )){
			StringBuffer sb = new StringBuffer();
			for(int i=0;i<this.repeats.length;i++){

				int[] idl = repeats[i].getRVDIndexAndLength( i<repeats.length-1 ? TALEConsensus.repeat : TALEConsensus.lastRepeat );

				this.repeats[i].rvd = repeats[i].getRVD( idl );
				this.repeats[i].maskedRepeats = repeats[i].getBlankedRepeats( idl );
				this.repeats[i].rvdPosition = idl[0];
				this.repeats[i].rvdLength = idl[1];
				/*System.out.println(this.repeats[i].rvd);
			System.out.println(this.repeats[i].maskedRepeat);
			System.out.println("+++++++++++++");*/
				sb.append(this.repeats[i].rvd);
				if(i<this.repeats.length-1){
					sb.append("-");
				}

			}
			this.rvdSequence = Sequence.create( RVDAlphabetContainer.SINGLETON, sb.toString(),"-" );
		}
	}
	
	private String parse(String id2) {
		String oldId = id2;
		
		//System.out.println("oldId: >>"+oldId+"<<");
		
		Pattern p = Pattern.compile("\\s*\\[([0-9]+)\\-([0-9]+)\\:(\\-?\\+?[0-9])\\]");
		Matcher m = p.matcher(oldId);
		
		if(m.find()){
			//System.out.println("matches");
			int start = Integer.parseInt( m.group(1) );
			int end = Integer.parseInt(m.group(2));
			boolean strand = Integer.parseInt(m.group(3)) >= 0;
			this.addAnnotation(null, null, start, end, strand);
			//System.out.println("added annotation");
			oldId = m.replaceFirst("");
		}
		return oldId;
	}

	public TALE(StringBuffer xml) throws NonParsableException{
		xml = XMLParser.extractForTag( xml, "TALE" );
		dnaOriginal = (TALE)XMLParser.extractObjectForTags( xml, "dnaOriginal" );
		id = (String)XMLParser.extractObjectForTags( xml, "id" );
		isNew = (Boolean)XMLParser.extractObjectForTags( xml, "isNew" );
		repeats = (Repeat[])XMLParser.extractObjectForTags( xml, "repeats" );
		try{
			end = XMLParser.extractSequencesWithTags( xml, "end" )[0];
			rvdSequence = XMLParser.extractSequencesWithTags( xml, "rvdSequence" )[0];
			if(rvdSequence != null){
				if(! Singleton.class.isAssignableFrom( rvdSequence.getAlphabetContainer().getClass() ) ){
					this.rvdSequence = Sequence.create( RVDAlphabetContainer.SINGLETON, rvdSequence.toString("-", 0, rvdSequence.getLength()),"-" );
				}
			}
			start = XMLParser.extractSequencesWithTags( xml, "start" )[0];
		}catch(WrongAlphabetException e){
			NonParsableException ex = new NonParsableException();
			ex.setStackTrace( e.getStackTrace() );
			throw ex;
		}
		try{
			strain = (String) XMLParser.extractObjectForTags(xml, "strain");
			accession = (String) XMLParser.extractObjectForTags(xml, "accession");
			startPos = (Integer) XMLParser.extractObjectForTags(xml, "startPos");
			endPos = (Integer) XMLParser.extractObjectForTags(xml, "endPos");
			strand = (Boolean) XMLParser.extractObjectForTags(xml, "strand");
			if(dnaOriginal != null){
				if(dnaOriginal.strain == null){dnaOriginal.strain = strain;}
				if(dnaOriginal.accession == null){dnaOriginal.accession = accession;}
				if(dnaOriginal.startPos == null){dnaOriginal.startPos = startPos;}
				if(dnaOriginal.endPos == null){dnaOriginal.endPos = endPos;}
				if(dnaOriginal.strand == null){dnaOriginal.strand = strand;}
			}
		}catch(NonParsableException e){

		}
	}
	
	public StringBuffer toXML(){
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, dnaOriginal, "dnaOriginal" );
		XMLParser.appendObjectWithTags( xml, id, "id" );
		XMLParser.appendObjectWithTags( xml, isNew, "isNew" );
		XMLParser.appendObjectWithTags( xml, repeats, "repeats" );
		XMLParser.appendSequencesWithTags( xml, "end", end );
		XMLParser.appendSequencesWithTags( xml, "rvdSequence", rvdSequence );
		XMLParser.appendSequencesWithTags( xml, "start", start );
		XMLParser.appendObjectWithTags(xml, strain, "strain");
		XMLParser.appendObjectWithTags(xml, accession, "accession");
		XMLParser.appendObjectWithTags( xml, startPos, "startPos");
		XMLParser.appendObjectWithTags(xml, endPos, "endPos");
		XMLParser.appendObjectWithTags(xml, strand, "strand");
		XMLParser.addTags( xml, "TALE" );
		return xml;
	}
	
	public TALE getTranslatedTALE(Translator t) throws IllegalArgumentException, WrongAlphabetException, WrongSequenceTypeException, EmptyDataSetException, IOException, DoubleSymbolException, OperationNotSupportedException{
		
		Repeat[] trans = new Repeat[repeats.length];
		for(int i=0;i<repeats.length;i++){
			Sequence temp = repeats[i].getRepeat();
			temp = Translator.DEFAULT.translate( temp, 0 );
			trans[i] = new Repeat( temp );
		}
		
		TALE tale = new TALE(id,t.translate( start, 0 ),trans,t.translate( end, 0 ));
		
		tale.setDnaOriginal( this );
		
		tale.accession = accession;
		tale.strain = strain;
		tale.startPos = startPos;
		tale.endPos = endPos;
		tale.strand = strand;
		
		return tale;
	}
	
	public String getId() {
		return id;
	}

	public Sequence getStart() {
		return start;
	}


	public int getNumberOfRepeats(){
		return repeats.length;
	}
	
	public Repeat getRepeat(int i){
		return repeats[i];
	}
	
	public Repeat[] getRepeats() {
		return repeats;
	}

	
	
	public Sequence getRVDSequence(Sequence consensus, Sequence lastConsensus) throws IllegalArgumentException, WrongAlphabetException, OperationNotSupportedException{
		
		StringBuffer reps = new StringBuffer();
		for(int i=0;i<repeats.length;i++){
			reps.append(repeats[i].getRVD( i<repeats.length-1 ? consensus : lastConsensus ));
			if(i<repeats.length){
				reps.append("-");
			}
		}
		
		return Sequence.create( RVDAlphabetContainer.SINGLETON, reps.toString(),"-" );
		
	}
	
	public Sequence getEnd() {
		return end;
	}
	
	public static TALE[] translateTALEs(TALE[] tales, Translator t) throws IllegalArgumentException, IOException, DoubleSymbolException, WrongAlphabetException, WrongSequenceTypeException, EmptyDataSetException, OperationNotSupportedException{
		
		TALE[] trans = new TALE[tales.length];
		for(int i=0;i<tales.length;i++){
			trans[i] = tales[i].getTranslatedTALE( t );
		}
		
		return trans;
	}
	
	public static TALE[] readTALEs(String startFile, String repeatFile, String endFile ) throws IllegalArgumentException, IOException, WrongAlphabetException, OperationNotSupportedException, WrongSequenceTypeException{
		
		HashMap<String, Sequence> starts = new HashMap<String, Sequence>();
		BufferedReader read = new BufferedReader(  new FileReader( startFile ) );
		String str = null, key = null;
		LinkedList<String> keyList = new LinkedList<String>();
		while( (str = read.readLine()) != null ){
			if(str.startsWith(">")){
				key = str.substring(1).trim();
				keyList.add( key );
			}else{
				starts.put( key, Sequence.create(DNAAlphabetContainer.SINGLETON, str.trim() ) );
			}
		}
		read.close();
		
		HashMap<String, Sequence> ends = new HashMap<String, Sequence>();
		read = new BufferedReader(  new FileReader( endFile ) );
		str = null; key = null;
		while( (str = read.readLine()) != null ){
			if(str.startsWith(">")){
				key = str.substring(1).trim();
			}else{
				ends.put( key, Sequence.create(DNAAlphabetContainer.SINGLETON, str.trim() ) );
			}
		}
		read.close();
		
		HashMap<String,Repeat[]> repeats = new HashMap<String,Repeat[]>();
		
		read = new BufferedReader(  new FileReader( repeatFile ) );
		str = null; key = null;
		LinkedList<Repeat> rep = new LinkedList<Repeat>();
		while( (str = read.readLine()) != null ){
			if(str.startsWith(">")){
				if(key != null){
					repeats.put( key, rep.toArray(new Repeat[0]) );
					rep.clear();
				}
				key = str.substring(1).trim();
			}else{
				rep.add( new Repeat( Sequence.create(DNAAlphabetContainer.SINGLETON, str.trim() ) ) );
			}
		}
		if(key != null){
			repeats.put( key, rep.toArray(new Repeat[0]) );
			rep.clear();
		}
		read.close();
		
		if(starts.size() != ends.size() || starts.size() != repeats.size()){
			System.out.println(starts.size()+" "+ends.size()+" "+repeats.size());
		}
		
		Iterator<String> keys = keyList.iterator(); 
		
		TALE[] tales = new TALE[starts.size()];
		
		int i=0;
		while(keys.hasNext()){
			key = keys.next();
			tales[i] = new TALE( key, starts.get( key ), repeats.get( key ), ends.get( key ) );
			i++;
		}
		return tales;
	}



	public String getRVD( int i, Sequence repeatCons ) {
		return repeats[i].getRvd();
	}


	public Sequence[] getBlankedRepeats( int i, Sequence repeatCons ) throws OperationNotSupportedException, WrongAlphabetException, WrongSequenceTypeException {
		return repeats[i].getBlankedRepeats( repeatCons );
	}



	public void setId( String newId ) {
		this.id = newId;
		if(this.dnaOriginal != null){
			this.dnaOriginal.setId( newId );;
		}
	}





	public void drawLeaf( Graphics2D graphics, int xoff, int yoff, int height, int margin ) {
		graphics = (Graphics2D)graphics.create();
		
		Rectangle2D rect = graphics.getFontMetrics().getStringBounds( id, graphics );
		Sequence rvds = getRvdSequence();
		Rectangle2D rect2 = graphics.getFontMetrics().getStringBounds( rvds.toString( "-", 0, rvds.getLength() ), graphics );
		
		double h = (rect.getHeight()+rect2.getHeight())*1.2;
		double rat = height/h;
		
		int size = graphics.getFont().getSize();
		size *= rat;
		
		graphics.setFont( new Font(graphics.getFont().getFontName(),graphics.getFont().getStyle(),size) );
		
		rect = graphics.getFontMetrics().getStringBounds( id, graphics );
		rect2 = graphics.getFontMetrics().getStringBounds( rvds.toString( "-", 0, rvds.getLength() ), graphics );
		
		int w = this.getWidth( graphics, height );
		if(this.isNew){
			graphics.setColor( Color.BLUE );
		}
		graphics.drawRect( xoff+margin/4, yoff, w+margin/2, height );
		
		graphics.drawString( id, xoff+margin/2, yoff-margin/4+(int)(height/2.0-rect.getCenterY()) );
		graphics.drawString( rvds.toString("-", 0, rvds.getLength()), xoff+margin/2, yoff-margin/4+(int)(height-rect2.getCenterY()) );
		
	}




	public int getWidth( Graphics2D graphics, int height ) {
		
		Rectangle2D rect = graphics.getFontMetrics().getStringBounds( id, graphics );
		Sequence rvds = getRvdSequence();
		Rectangle2D rect2 = graphics.getFontMetrics().getStringBounds( rvds.toString( "-", 0, rvds.getLength() ), graphics );
		
		double h = (rect.getHeight()+rect2.getHeight())*1.2;
		double rat = height/h;
		
		double w = Math.max( rect.getWidth(), rect2.getWidth() );
		w *= rat;
		
		return (int)w;
	}





	public void setIsNew( boolean b ) {
		this.isNew = b;
	}





	public boolean isNew() {
		return this.isNew;
	}





	public String getCodon(int repeat, int position) {
		Repeat r = repeats[repeat];
		if(r.getMaskedRepeats() != null){
			int before = r.getMaskedRepeats()[0].getLength()*3;
			TALE dna = getDnaOriginal();
			if(dna == null){
				return null;
			}
			String rvd = r.getRvd();
			if(position == 0){
				return dna.getRepeat(repeat).getRepeat().getSubSequence(before, 3).toString();
			}else if(position == 1){
				if(rvd.endsWith("*")){
					return "---";
				}else{
					return dna.getRepeat(repeat).getRepeat().getSubSequence(before+3, 3).toString();
				}
			}else{
				return null;
			}
		}else{
			return null;
		}
	}





	public boolean containsAberrantRepeat() {
		for(int i=0;i<repeats.length-1;i++){
			if(repeats[i].getType() != Type.UNKNOWN && repeats[i].getType() != Type.NORMAL){
				return true;
			}
		}
		return false;
	}

	public ResultSet annotationToResultSet() {
		return  new ResultSet(new Result[]{
				new CategoricalResult("ID", "", getId()),
				new CategoricalResult("Strain", "", strain == null ? "" : strain),
				new CategoricalResult("Accession", "", accession == null ? "" : accession),
				new NumericalResult("Start", "", startPos == null ? -1 : startPos ),
				new NumericalResult("End", "", endPos == null ? -1 : endPos),
				new CategoricalResult("Strand", "", strand == null ? "" : (strand ? "+1" : "-1"))
		});
	}

	public void setAccession(String accession2) {
		this.accession = accession2;
		if(dnaOriginal != null){
			dnaOriginal.accession = accession2;
		}
	}

	public String getStrain() {
		return strain;
	}

	public String getAccession() {
		return accession;
	}
	
}
