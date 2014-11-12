package de.jstacs.data;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.data.GFFParser.GFFEntry.GFFType;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;


public class GFFParser {

	public static void main(String[] args) throws Exception {
		
		BufferedReader read = new BufferedReader(new FileReader( "/Users/dev/Downloads/Galaxy43-[UCSC_Main_on_Human__knownGene_(genome)].gtf.txt" ));
		
		GFFList list = GFFParser.parseGTF( read, null, null );
		
		/*BufferedReader read = new BufferedReader(new FileReader( "/Users/dev/Desktop/TAL-Chips/GA/ath/TAIR10_GFF3_genes.gff" ));
		
		
		
		DataSet set = new DataSet( new AlphabetContainer( new DiscreteAlphabet( true, "A","C","G","T","N", "Y","W","M","K","S","R","D" ) ), new SparseStringExtractor( "/Users/dev/Desktop/TAL-Chips/GA/ath/TAIR10.fa", '>', new SimpleSequenceAnnotationParser() ) );
		
		GFFList list = new GFFList( read );
		
		String[] keys = new String[set.getNumberOfElements()];
		for(int i=0;i<set.getNumberOfElements();i++){
			keys[i] = (String)set.getElementAt( i ).getAnnotation()[0].getResultAt( 0 ).getValue();
			keys[i] = keys[i].substring( 0, keys[i].indexOf( ' ' ) );
			System.out.println(keys[i]);
		}
		
		list.attachData( set, keys );
		*/
		list = list.getSubList( "chr4" );
		
		list = list.getSubList( GFFType.exon );
		
		ArrayList<GFFEntry> ens = list.getEntriesOverlapping( 196560 );
		
		System.out.println(ens);
		
		//System.out.println(list.getSequenceFor( ens.get( 0 ), 2, 2 ) );
		
	}
	
	public static GFFList parse(String filename) throws IOException{
		int idx = filename.lastIndexOf( '.' );
		if(idx >= 0){
			return parse(filename,filename.substring( idx+1 ));
		}else{
			return parse(filename,null);
		}
	}
	
	public static GFFList parse(String filename, String extension) throws IOException {
		if(extension == null || extension.equalsIgnoreCase( "gff" ) || extension.equalsIgnoreCase( "gff3" ) ){
			return parseGFF( filename );
		}else if(extension.equalsIgnoreCase( "gtf" )){
			return parseGTF( new BufferedReader( new FileReader( filename ) ), null, null );
		}else{
			return parseUCSCKnownGenesBED( new BufferedReader( new FileReader( filename ) ), null, null );
		}
	}
	
	public static GFFList parseGFF(String filename) throws IOException{
		return parseGFF(filename, null, null);
	}
	
	public static GFFList parseGFF(String filename, DataSet data, String[] keys) throws IOException{
		BufferedReader reader = new BufferedReader( new FileReader( filename ) );
		return parseGFF(reader, data, keys);
	}
	
	public static GFFList parseGFF(BufferedReader reader) throws IOException{
		return parseGFF(reader,null,null);
	}
	
	public static GFFList parseGFF(BufferedReader reader, DataSet data, String[] keys) throws IOException{
		GFFList list = new GFFList( reader );		
		
		if(data != null){
			list.attachData( data, keys );		
		}
		return list;
		
	}
	
	public static GFFList parseGTF(BufferedReader reader, DataSet data, String[] keys) throws IOException {
		LinkedList<GFFEntry> list = new LinkedList<GFFParser.GFFEntry>();
		String str = null;
		StringBuffer sb = new StringBuffer();
		while( (str = reader.readLine()) != null ){
			int idx = str.lastIndexOf( '\t' );
			String first = str.substring( 0, idx );
			String second = str.substring( idx+1 );
			String[] parts = second.split( ";" );
			sb.delete( 0, sb.length() );
			for(int i=0;i<parts.length;i++){
				parts[i] = parts[i].trim();
				if(parts[i].length() > 0){
					idx = parts[i].indexOf( ' ' );
					sb.append( parts[i].substring( 0, idx ) );
					sb.append( "=" );
					sb.append( parts[i].substring( idx+1 ) );
					if(i < parts.length-1){
						sb.append( ";" );
					}
				}
				
			}
			list.add( new GFFEntry( first+"\t"+sb.toString() ) );
		}
		GFFList gff = new GFFList( list );
		if(data != null){
			gff.attachData( data, keys );		
		}
		return gff;
	}
	
	public static GFFList parseUCSCKnownGenesBED(BufferedReader reader, DataSet data, String[] keys) throws IOException {
		
		ArrayList<GFFEntry> entries = new ArrayList<GFFParser.GFFEntry>();
		LinkedList<GFFEntry> subs = new LinkedList<GFFParser.GFFEntry>();
		String str = null;
		while( (str = reader.readLine()) != null ){
			String[] parts = str.split( "\t" );
			String seqid = parts[0];
			int start = Integer.parseInt( parts[1] )+1;
			int end = Integer.parseInt( parts[2] )+1;
			String id = parts[3];
			double score = Double.parseDouble( parts[4] );
			parts[5] = parts[5].trim();
			Strand strand;
			if(parts[5].length() == 0 ||parts[5].equals( "." )){
				strand = Strand.UNKNOWN;
			}else if(parts[5].equals( "+" )){
				strand = Strand.FORWARD;
			}else if(parts[5].equals( "-" )){
				strand = Strand.REVERSE;
			}else{
				strand = Strand.UNKNOWN;
			}
			int cdsStart = Integer.parseInt( parts[6] );
			int cdsEnd = Integer.parseInt( parts[7] );
			int numExons = Integer.parseInt( parts[9] );
			GFFEntry gene = new GFFEntry( seqid, "known_genes", GFFType.gene.name(), start, end, score, strand, -1, "" );
			gene.id = id;
			entries.add( gene );
			String[] lengthParts = parts[10].split( "," );
			String[] startParts = parts[11].split( "," );
			subs.clear();
			for(int i=0;i<numExons;i++){
				int exOff = Integer.parseInt( startParts[i] );
				int exLen = Integer.parseInt( lengthParts[i] );
				int exStart = start + exOff;
				int exEnd = exStart + exLen;
				GFFEntry exon = new GFFEntry( seqid, "known_genes", GFFType.exon.name(), exStart, exEnd, score, strand, -1, "" );
				subs.add( exon );
				if(cdsStart != cdsEnd){
					
					if(exStart>=cdsStart && exEnd <= cdsEnd){
						GFFEntry cds = new GFFEntry( seqid, "known_genes", GFFType.CDS.name(), exStart, exEnd, score, strand, -1, "" );
						subs.add( cds );
					}else{
						if(exStart < cdsStart){
							if(exEnd < cdsStart){
								GFFEntry utr5 = new GFFEntry( seqid, "known_genes", GFFType.five_prime_UTR.name(), exStart, exEnd, score, strand, -1, "" );
								subs.add( utr5 );
							}else{
								GFFEntry utr5 = new GFFEntry( seqid, "known_genes", GFFType.five_prime_UTR.name(), exStart, cdsStart-1, score, strand, -1, "" );
								GFFEntry cds = new GFFEntry( seqid, "known_genes", GFFType.CDS.name(), cdsStart, exEnd, score, strand, -1, "" );
								subs.add( utr5 );
								subs.add( cds );
							}
						}else if(exEnd > cdsEnd){
							if(exStart > cdsEnd){
								GFFEntry utr3 = new GFFEntry( seqid, "known_genes", GFFType.three_prime_UTR.name(), exStart, exEnd, score, strand, -1, "" );
								subs.add( utr3 );
							}else{
								GFFEntry utr3 = new GFFEntry( seqid, "known_genes", GFFType.three_prime_UTR.name(), cdsEnd+1, exEnd, score, strand, -1, "" );
								GFFEntry cds = new GFFEntry( seqid, "known_genes", GFFType.CDS.name(), exStart, cdsEnd, score, strand, -1, "" );
								subs.add( utr3 );
								subs.add( cds );
							}
						}
					}
				}
			}
			Iterator<GFFEntry> it = subs.iterator();
			while(it.hasNext()){
				GFFEntry sub = it.next();
				sub.parent = new LinkedList<GFFParser.GFFEntry>();
				sub.parent.add( gene );
				if(gene.children == null){
					gene.children = new LinkedList<GFFParser.GFFEntry>();
				}
				gene.children.add( sub );
				entries.add( sub );
			}
		}
		GFFList list = new GFFList( entries );
		if(data != null){
			list.attachData( data, keys );
		}
		return list;
	}
	
	private static class GFFEntryComparator implements Comparator<GFFEntry>{

		private static final GFFEntryComparator COMP = new GFFEntryComparator();
		
		@Override
		public int compare( GFFEntry o1, GFFEntry o2 ) {
			int i = o1.seqid.compareTo( o2.seqid );
			if(i == 0){
				i = o1.start > o2.start ? 1 : (o1.start < o2.start ? -1 : 0);
			}
			if(i == 0){
				i = o1.end > o2.end ? 1 : (o1.end < o2.end ? -1 : 0);
			}
			return i;
		}

		
	}
	
	public static class GFFList{
		
		private ArrayList<GFFEntry> entries;
		private HashMap<String, GFFEntry> byID;
		private HashMap<String, GFFList> bySeq;
		private HashMap<GFFParser.GFFEntry.GFFType, GFFList> byType;
		private HashMap<String, Sequence> data;
		
		
		public GFFList(Collection<GFFEntry> entries){
			this.entries = new ArrayList<GFFParser.GFFEntry>();
			this.entries.addAll( entries );
			parse(true);
		}
		
		private GFFList(ArrayList<GFFEntry> entries){
			this.entries = entries;
			parse(false);
		}
		
		public void attachData(DataSet data, String[] keys){
			if(keys.length != data.getNumberOfElements()){
				throw new IllegalArgumentException();
			}
			Iterator<String> it = bySeq.keySet().iterator();
			this.data = new HashMap<String, Sequence>();
			outerloop:
			while(it.hasNext()){
				String str = it.next();
				for(int i=0;i<keys.length;i++){
					if(str.equals( keys[i] )){
						this.data.put( keys[i], data.getElementAt( i ) );
						if(bySeq.get( keys[i] ).data == null){
							bySeq.get( keys[i] ).attachData( data, keys );
						}
						continue outerloop;
					}
				}
				throw new RuntimeException();
			}
			
			Iterator<GFFType> it2 = byType.keySet().iterator();
			while(it2.hasNext()){
				GFFType type = it2.next();
				if(byType.get( type ).data == null){
					byType.get( type ).attachData( data, keys );
				}
			}
		}
		
		public Sequence getSequenceFor(GFFEntry entry, int offleft, int offright){
			Sequence seq = data.get( entry.getSeqid() ).getSubSequence( entry.getStart()-offleft-1, entry.getEnd()-entry.getStart()+1+offleft+offright );
			return seq;
		}
		
		private void parse(boolean parseParents) {
			
			HashMap<String, ArrayList<GFFEntry>> bySeq = new HashMap<String, ArrayList<GFFEntry>>();
			HashMap<GFFParser.GFFEntry.GFFType,ArrayList<GFFEntry>> byType = new HashMap<GFFParser.GFFEntry.GFFType, ArrayList<GFFEntry>>();
			
			byID = new HashMap<String, GFFParser.GFFEntry>();
			
			for(int i=0;i<entries.size();i++){
				GFFEntry en = entries.get( i );
				if(en.getID() != null){
					byID.put( en.getID(), en );
				}
				if(bySeq.get( en.getSeqid() ) == null){
					bySeq.put( en.getSeqid(), new ArrayList<GFFParser.GFFEntry>() );
				}
				ArrayList<GFFEntry> list = bySeq.get( en.getSeqid() );
				list.add( en );
				if(byType.get( en.getParsedType() ) == null){
					byType.put( en.getParsedType(), new ArrayList<GFFEntry>() );
				}
				list = byType.get( en.getParsedType() );
				list.add( en );
			}
			if(parseParents){
				for(int i=0;i<entries.size();i++){
					GFFEntry en = entries.get( i );
					en.parseParent(byID);
				}
			}
			sort(bySeq);
			sort(byType);
			
			this.bySeq = new HashMap<String, GFFParser.GFFList>();
			this.byType = new HashMap<GFFParser.GFFEntry.GFFType, GFFParser.GFFList>();
			
			if(bySeq.keySet().size() > 1){
				descend( bySeq, this.bySeq );
			}else{
				this.bySeq.put( bySeq.keySet().iterator().next(), this );
			}
			
			if(byType.keySet().size() > 1){
				descend( byType, this.byType );
			}else{
				this.byType.put( byType.keySet().iterator().next(), this );
			}
			
		}

		private void descend(HashMap map, HashMap map2){
			Iterator<? extends Object> it = map.keySet().iterator();
			while(it.hasNext()){
				Object key = it.next();
				ArrayList<GFFEntry> en = (ArrayList<GFFParser.GFFEntry>)map.get( key );
				map2.put( key, new GFFList( en ) );
			}
		}
		
		private void sort( HashMap<? extends Object, ArrayList<GFFEntry>> bySeq2 ) {
			Iterator<? extends Object> it = bySeq2.keySet().iterator();
			while(it.hasNext()){
				ArrayList<GFFEntry> en = bySeq2.get( it.next() );
				GFFEntry[] ens = en.toArray( new GFFEntry[0] );
				Arrays.sort( ens, GFFEntryComparator.COMP);
				en.clear();
				for(int i=0;i<ens.length;i++){
					en.add( ens[i] );
				}
			}
		}

		public GFFList(BufferedReader reader) throws IOException{
			entries = new ArrayList<GFFParser.GFFEntry>();
			String line;
			while( (line = reader.readLine()) != null ){
				entries.add( new GFFEntry( line ) );
			}
			parse(true);
		}
		
		public GFFList getSubList(String seqId){
			return bySeq.get( seqId );
		}
		
		public GFFList getSubList(de.jstacs.data.GFFParser.GFFEntry.GFFType type){
			return byType.get( type );
		}
		
		public ArrayList<GFFEntry> getEntriesBetween(int start, int end){
			if(bySeq.keySet().size() > 1){
				throw new RuntimeException( "get sublist for seqid first" );
			}
			
			
			ArrayList<GFFEntry> res = new ArrayList<GFFEntry>();
			
			int idx = binarySearch(entries, start);
			for(int i=idx;i<entries.size();i++){
				GFFEntry en = entries.get( i );
				if(en.start > end){
					break;
				}else{
					if(en.start >= start && en.end <= end){
						res.add( en );
					}
				}
			}
			return res;
		}
		
		
		
		public ArrayList<GFFEntry> getEntriesOverlapping(int start, int end){
			if(bySeq.keySet().size() > 1){
				throw new RuntimeException( "get sublist for seqid first" );
			}
			
			int idx2 = binarySearch(entries,end);
			ArrayList<GFFEntry> res = new ArrayList<GFFEntry>();
			
			for(int i=0;i<idx2;i++){
				GFFEntry en = entries.get( i );
				if(en.start <= start && en.end >= start || en.start <= end && en.end >= end || en.start >= start && en.end <= end){
					res.add( en );
				}
			}
			
			return res;
		}
		
		public ArrayList<GFFEntry> getEntriesOverlapping(int idx){
			if(bySeq.keySet().size() > 1){
				throw new RuntimeException( "get sublist for seqid first" );
			}
			
			int idx2 = binarySearch(entries,idx);
			ArrayList<GFFEntry> res = new ArrayList<GFFEntry>();
			
			for(int i=0;i<idx2;i++){
				GFFEntry en = entries.get( i );
				if(en.start <= idx && en.end >= idx){
					res.add( en );
				}
			}
			
			return res;
		}
		
		private int binarySearch(ArrayList<GFFEntry> en, int idx){
			return binarySearch(en,idx,0,en.size());
		}
		
		private int binarySearch(ArrayList<GFFEntry> en, int idx, int left, int right){
			if(left == right){
				return left;
			}
			int pivot = en.get( (left+right)/2 ).start;
			if(pivot >= idx){
				return binarySearch(en,idx,left,(left+right)/2);
			}else{
				return binarySearch(en,idx,(left+right)/2+1,right);
			}
		}

		
		
		public ListResult toListResult(String name, String comment, boolean gff3) {
			
			ArrayList<ResultSet> set = new ArrayList<ResultSet>(); 
			
			Iterator<GFFEntry> it = entries.iterator();
			while(it.hasNext()){
				GFFEntry en = it.next();
				String strand = en.getStrand() == Strand.UNKNOWN ? "." : (en.getStrand() == Strand.FORWARD ? "+" : "-");
				String phase = en.getPhase() == -1 ? "." : en.getPhase()+"";
				set.add( new ResultSet(new Result[]{
				                                    new CategoricalResult( "SequenceID", "", en.getSeqid() ),
				                                    new CategoricalResult( "Source", "", en.getSource() ),
				                                    new CategoricalResult( "Type", "", en.getType() ),
				                                    new NumericalResult( "Start", "", en.getStart() ),
				                                    new NumericalResult( "End", "", en.getEnd() ),
				                                    new NumericalResult( "Score", "", en.getScore() ),
				                                    new CategoricalResult( "Strand", "", strand ),
				                                    new CategoricalResult( "Phase", "", phase ),
				                                    new CategoricalResult( "Attributes", "", en.getAttributeString(gff3) )
				                                    
				}) );
			}
			
			return new ListResult( name, comment, null, set.toArray( new ResultSet[0] ) );
		}
		
	}
	
	public static class GFFEntry{
		
		public enum GFFType{
			chromosome,
			gene,
			mRNA,
			protein,
			exon,
			five_prime_UTR,
			CDS,
			three_prime_UTR,
			other;
			
			public static GFFType fromString(String text) {
			    if (text != null) {
			      for (GFFType b : GFFType.values()) {
			        if (text.equalsIgnoreCase(b.name())) {
			          return b;
			        }
			      }
			    }
			    return null;
			  }
		}
		
		private String seqid;
		private String source;
		private String type;
		private GFFType parsedType;
		private int start;
		private int end;
		private String id;
		
		
		public String getID() {
			return id;
		}


		private String getAttributeString( boolean gff3 ) {
			
			if(gff3){
				StringBuffer sb = new StringBuffer();
				Iterator<String> keys = attributes.keySet().iterator();
				while(keys.hasNext()){
					String key = keys.next();
					sb.append( key+"="+attributes.get( key )+";" );
				}

				return sb.toString();
			}else{
				if(attributes.containsKey( "group" )){
					return attributes.get( "group" );
				}else{
					return "";
				}
			}
		}


		private void parseParent( HashMap<String, GFFEntry> byID ) {
			
			String parent = attributes.get( "parent" );
			if(parent != null){
				this.parent = new LinkedList<GFFParser.GFFEntry>();
				String[] parents = parent.split( "," );
				for(int i=0;i<parents.length;i++){
					GFFEntry par = byID.get( parents[i] );
					this.parent.add( par );
					if(par.children == null){
						par.children = new LinkedList<GFFParser.GFFEntry>();
					}
					par.children.add( this );
				}
			}
			
		}


		private double score;
		private Strand strand;
		private int phase;
		private HashMap<String, String> attributes;
		
		private LinkedList<GFFEntry> parent;
		private LinkedList<GFFEntry> children;
		
		/**
		 * @param seqid
		 * @param source
		 * @param type
		 * @param start
		 * @param end
		 * @param score
		 * @param strand
		 * @param phase
		 * @param attributes
		 */
		public GFFEntry( String seqid, String source, String type, int start, int end, double score, Strand strand, int phase,
							String attributes ) {
			this.seqid = seqid;
			this.source = source;
			this.type = type;
			parsedType = GFFType.fromString( type );
			this.start = start;
			this.end = end;
			this.score = score;
			this.strand = strand;
			this.phase = phase;
			parseAttributes( attributes );
		}
		
		private void setID() {
			id = attributes.get( "ID" );
		}

		public GFFEntry(String line){
			String[] parts = line.split( "\t" );
			seqid = parts[0];
			source = parts[1];
			type = parts[2];
			parsedType = GFFType.fromString( type );
			parts[3] = parts[3].trim();
			if(parts[3].length() > 0 && !parts[3].equals( "." )){
				start = Integer.parseInt( parts[3] );
			}else{
				start = -1;
			}
			parts[4] = parts[4].trim();
			if(parts[4].length() > 0 && !parts[4].equals( "." )){
				end = Integer.parseInt( parts[4] );
			}else{
				end = -1;
			}
			parts[5] = parts[5].trim();
			if(parts[5].length() > 0 && !parts[5].equals( "." )){
				score = Double.parseDouble( parts[5] );
			}else{
				score = Double.NaN;
			}
			parts[6] = parts[6].trim();
			if(parts[6].length() == 0 ||parts[6].equals( "." )){
				strand = Strand.UNKNOWN;
			}else if(parts[6].equals( "+" )){
				strand = Strand.FORWARD;
			}else if(parts[6].equals( "-" )){
				strand = Strand.REVERSE;
			}else{
				strand = Strand.UNKNOWN;
			}
			
			parts[7] = parts[7].trim();
			if(parts[7].length() > 0 && !parts[7].equals( "." )){
				phase = Integer.parseInt( parts[7] );
			}else{
				phase = -1;
			}		
			
			parseAttributes( parts[8] );
			
			
		}
		
		private void parseAttributes(String attrs){
			attributes = new HashMap<String, String>();
			if(attrs == null || attrs.length() == 0){
				return;
			}
			if(attrs.indexOf( "=" ) < 0){//GFF2
				attributes.put( "group", attrs );
				return;
			}else{//GFF3
				String[] parts = attrs.split( ";" );
				for(int i=0;i<parts.length;i++){
					String[] temp = parts[i].split( "=" );
					attributes.put( temp[0].trim(), temp[1].trim() );
				}
				setID();
			}
		}
		
		private boolean hasAttribute(String key){
			return attributes != null && attributes.get( key ) != null;
		}
		
		private String getAttribute(String key){
			if(attributes == null){
				return null;
			}else{
				return attributes.get( key );
			}
		}
		
		public String getSeqid() {
			return seqid;
		}

		
		public String getSource() {
			return source;
		}

		
		public String getType() {
			return type;
		}

		
		public GFFType getParsedType() {
			return parsedType;
		}

		
		public int getStart() {
			return start;
		}

		
		public int getEnd() {
			return end;
		}

		
		public double getScore() {
			return score;
		}

		
		public Strand getStrand() {
			return strand;
		}

		
		public int getPhase() {
			return phase;
		}

		
		public HashMap<String, String> getAttributes() {
			return attributes;
		}

		
		public LinkedList<GFFEntry> getParents() {
			return parent;
		}

		
		public LinkedList<GFFEntry> getChildren() {
			return children;
		}
		
		
		public String toString(){
			StringBuffer str = new StringBuffer( type );
			if(id != null){
				str.append( " " );
				str.append( id );
			}
			str.append( " ("+start+", "+end+")" );
			if(parent != null && parent.size() > 0){
				str.append(  " of ["+parent.get( 0 ) );
				for(int i=1;i<parent.size();i++){
					str.append( " and "+parent.get( i ) );
				}
				str.append( "]" );
			}
			return str.toString();
		}
		
	}
	
	
	
}
