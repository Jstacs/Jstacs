/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 *
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package projects.talen;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import projects.talen.GFFParser.GFFEntry.GFFType;

/**
 * Class for a rudimentary GFF parser.
 * @author Jan Grau
 *
 */
public class GFFParser {
	
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
	
	public static GFFList getCDSandAddParents(GFFList list, String transcriptID, String geneID, String transcriptPrefix, String genePrefix, boolean stripOtherAttributes){
		
		HashMap<String, GFFEntry> transcripts = new HashMap<String, GFFParser.GFFEntry>();
		HashMap<String, GFFEntry> genes = new HashMap<String, GFFParser.GFFEntry>();
		
		
		GFFList subCDS = list.getSubList( GFFType.CDS, GFFType.stop_codon );
		ArrayList<GFFEntry> cds = subCDS.entries;
		
		
		for(int i=0;i<cds.size();i++){
			GFFEntry curr = cds.get( i );
			
			String transcript = transcriptPrefix+curr.attributes.get( transcriptID );
			String gene = genePrefix+curr.attributes.get( geneID );
			
			if(transcripts.containsKey( transcript )){
				GFFEntry te = transcripts.get( transcript );
				if(te.start > curr.start){
					te.start = curr.start;
				}
				if(te.end < curr.end){
					te.end = curr.end;
				}
			}else{
				GFFEntry te = new GFFEntry( curr.seqid, curr.source, "mRNA",curr.start, curr.end, Double.NaN, curr.strand, -1, "ID="+transcript+";Parent="+gene );
				transcripts.put( transcript, te );
			}
			
			if(genes.containsKey( gene )){
				GFFEntry ge = genes.get( gene );
				if(ge.start > curr.start){
					ge.start = curr.start;
				}
				if(ge.end < curr.end){
					ge.end = curr.end;
				}
			}else{
				GFFEntry ge = new GFFEntry( curr.seqid, curr.source, "gene",curr.start, curr.end, Double.NaN, curr.strand, -1, "ID="+gene );
				genes.put( gene, ge );
			}
			
			GFFEntry te = transcripts.get( transcript );
			GFFEntry ge = genes.get( gene );
			
			curr.attributes.remove( transcriptID );
			curr.attributes.remove( geneID );
			
			if(stripOtherAttributes){
				curr.attributes.clear();
			}
			
			curr.attributes.put( "Parent", transcript );
			
			if(curr.parent == null){
				curr.parent = new LinkedList<GFFParser.GFFEntry>();
			}
			curr.parent.add( te );
			if(te.children == null){
				te.children = new LinkedList<GFFParser.GFFEntry>();
			}
			te.children.add( curr );
			if(te.parent == null){
				te.parent = new LinkedList<GFFParser.GFFEntry>();
			}
			te.parent.add( ge );
			if(ge.children == null){
				ge.children = new LinkedList<GFFParser.GFFEntry>();
			}
			if(!ge.children.contains( te )){
				ge.children.add( te );
			}
		}
		
		
		transcripts.clear();
		cds.clear();
		
		GFFEntry[] geneAr = genes.values().toArray( new GFFEntry[0] );
		genes.clear();
		
		LinkedList<GFFEntry> all = new LinkedList<GFFParser.GFFEntry>();
		
		Arrays.sort( geneAr, GFFEntryComparator.COMP );
				
		for(int i=0;i<geneAr.length;i++){
			GFFEntry gene = geneAr[i];
			all.add( gene );
			Iterator<GFFEntry> transIt = gene.children.iterator();
			while(transIt.hasNext()){
				GFFEntry transcript = transIt.next();
				
				all.add( transcript );
				
				GFFEntry[] children = transcript.children.toArray( new GFFEntry[0] );
				Arrays.sort( children, new GFFEntryComparator( transcript.strand == Strand.REVERSE ) );
				
				int k=children.length-1;
				while(k>0 && children[k].parsedType == GFFType.stop_codon){
					if(children[k-1].parsedType == GFFType.CDS){
						if(transcript.strand == Strand.FORWARD && children[k].start == children[k-1].end+1){
							children[k-1].end = children[k].end;
							children[k] = null;
						}else if(transcript.strand == Strand.REVERSE && children[k].end == children[k-1].start-1){
							children[k-1].start = children[k].start;
							children[k] = null;
						}
					}
					if(children[k] != null){
						children[k].parsedType = GFFType.CDS;
						children[k].type = "CDS";
					}
					k--;
				}
				
				for(int j=0;j<children.length;j++){
					if(children[j] != null){
						all.add( children[j] );
					}
				}
				
			}
			/*GFFEntry[] transcriptAr = gene.children.toArray( new GFFEntry[0] );
			Arrays.sort( transcriptAr, GFFEntryComparator.COMP );
			for(int j=0;j<transcriptAr.length;j++){
				GFFEntry transcript = transcriptAr[j];
				all.add( transcript );
				GFFEntry[] cdsAr = transcript.children.toArray( new GFFEntry[0] );
				Arrays.sort( cdsAr, GFFEntryComparator.COMP );
				for(int k=0;k<cdsAr.length;k++){
					all.add( cdsAr[k] );
				}
			}*/
		}
		
		return new GFFList( all );
		
	}
	
	public static GFFList parseEnsembl(BufferedReader reader) throws NumberFormatException, IOException{
		String str = reader.readLine();
		String[] header = str.split( "\t" );
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		for(int i=0;i<header.length;i++){
			map.put( header[i], i );
		}
		
		ArrayList<GFFEntry> cds = new ArrayList<GFFParser.GFFEntry>();
		HashMap<String, GFFEntry> transcripts = new HashMap<String, GFFParser.GFFEntry>();
		HashMap<String, GFFEntry> genes = new HashMap<String, GFFParser.GFFEntry>();
		
		while( (str = reader.readLine()) != null ){
			String[] parts = str.split( "\t" );
			
			int startIdx = map.get( "Genomic coding start" );
			int endIdx = map.get( "Genomic coding end" );
			
			
			if(startIdx < parts.length && endIdx < parts.length && 
					parts[ startIdx ].trim().length() > 0 && parts[ endIdx ].trim().length() > 0){
				int start = Integer.parseInt( parts[ startIdx ] );
				int end = Integer.parseInt( parts[ endIdx ] );
				String gene = parts[ map.get( "Ensembl Gene ID" ) ];
				String transcript = parts[ map.get( "Ensembl Transcript ID" ) ];
				int geneStart = Integer.parseInt( parts[ map.get( "Gene Start (bp)" ) ] );
				int geneEnd = Integer.parseInt( parts[ map.get( "Gene End (bp)" ) ] );
				int transcriptStart = Integer.parseInt( parts[ map.get( "Transcript Start (bp)" ) ] );
				int transcriptEnd = Integer.parseInt( parts[ map.get( "Transcript End (bp)" ) ] );
				Strand strand = "-1".equals( parts[ map.get( "Strand" ) ] ) ? Strand.REVERSE : Strand.FORWARD;
				String chromosome = parts[ map.get( "Chromosome Name" ) ];
				int phase = Integer.parseInt( parts[ map.get( "phase" ) ] );
				
				if(!genes.containsKey( gene )){
					GFFEntry geneEn = new GFFEntry( chromosome, "ensembl", GFFType.gene.name(), geneStart, geneEnd, Double.NaN, strand, 0, "" );
					geneEn.setID(gene);
					genes.put( gene, geneEn );
				}
				GFFEntry geneEn = genes.get( gene );
				if(geneEn.start != geneStart || geneEn.end != geneEnd || !geneEn.seqid.equals( chromosome ) || geneEn.strand != strand){
					throw new RuntimeException();
				}
				if(!transcripts.containsKey( transcript )){
					GFFEntry transcriptEn = new GFFEntry( chromosome, "ensembl", GFFType.mRNA.name(), transcriptStart, transcriptEnd, Double.NaN, strand, 0, "" );
					transcriptEn.addParent( geneEn );
					geneEn.addChild( transcriptEn );
					transcriptEn.setID( transcript );
					transcripts.put( transcript, transcriptEn );
					
				}
				GFFEntry transcriptEn = transcripts.get( transcript );
				if(transcriptEn.start != transcriptStart || transcriptEn.end != transcriptEnd || !transcriptEn.seqid.equals( chromosome ) || transcriptEn.strand != strand){
					throw new RuntimeException();
				}
				
				GFFEntry cdsEn = new GFFEntry( chromosome, "ensembl", GFFType.CDS.name(), start, end, Double.NaN, strand, phase, "" );
				cdsEn.addParent( transcriptEn );
				transcriptEn.addChild( cdsEn );
				
				cds.add( cdsEn );
				
			}
		}
		
		
		transcripts.clear();
		cds.clear();
		
		GFFEntry[] geneAr = genes.values().toArray( new GFFEntry[0] );
		genes.clear();
		
		LinkedList<GFFEntry> all = new LinkedList<GFFParser.GFFEntry>();
		
		Arrays.sort( geneAr, GFFEntryComparator.COMP );
				
		for(int i=0;i<geneAr.length;i++){
			GFFEntry gene = geneAr[i];
			all.add( gene );
			Iterator<GFFEntry> transIt = gene.children.iterator();
			while(transIt.hasNext()){
				GFFEntry transcript = transIt.next();
				
				all.add( transcript );
				
				GFFEntry[] children = transcript.children.toArray( new GFFEntry[0] );
				Arrays.sort( children, new GFFEntryComparator( transcript.strand == Strand.REVERSE ) );
				
				int k=children.length-1;
				while(k>0 && children[k].parsedType == GFFType.stop_codon){
					if(children[k-1].parsedType == GFFType.CDS){
						if(transcript.strand == Strand.FORWARD && children[k].start == children[k-1].end+1){
							children[k-1].end = children[k].end;
							children[k] = null;
						}else if(transcript.strand == Strand.REVERSE && children[k].end == children[k-1].start-1){
							children[k-1].start = children[k].start;
							children[k] = null;
						}
					}
					if(children[k] != null){
						children[k].parsedType = GFFType.CDS;
						children[k].type = "CDS";
					}
					k--;
				}
				
				for(int j=0;j<children.length;j++){
					if(children[j] != null){
						all.add( children[j] );
					}
				}
				
			}
			/*GFFEntry[] transcriptAr = gene.children.toArray( new GFFEntry[0] );
			Arrays.sort( transcriptAr, GFFEntryComparator.COMP );
			for(int j=0;j<transcriptAr.length;j++){
				GFFEntry transcript = transcriptAr[j];
				all.add( transcript );
				GFFEntry[] cdsAr = transcript.children.toArray( new GFFEntry[0] );
				Arrays.sort( cdsAr, GFFEntryComparator.COMP );
				for(int k=0;k<cdsAr.length;k++){
					all.add( cdsAr[k] );
				}
			}*/
		}
		
		return new GFFList( all );
		
		
	}
	
	public static GFFList parseGTF(BufferedReader reader, DataSet data, String[] keys) throws IOException {
		LinkedList<GFFEntry> list = new LinkedList<GFFParser.GFFEntry>();
		String str = null;
		StringBuffer sb = new StringBuffer();
		while( (str = reader.readLine()) != null ){
			if(!str.startsWith("#!")){
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
						String temp = parts[i].substring( idx+1 ).trim();
						if((temp.startsWith( "\"" ) || temp.startsWith( "\'" ))
								&& (temp.endsWith( "\"" ) || temp.endsWith("\'"))){
							temp = temp.substring( 1, temp.length()-1 );
						}
						sb.append( temp );
						if(i < parts.length-1){
							sb.append( ";" );
						}
					}

				}
				list.add( new GFFEntry( first+"\t"+sb.toString() ) );
			}
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
		
		private boolean reverse;
		
		public GFFEntryComparator(){
			reverse = false;
		}
		
		public GFFEntryComparator(boolean reverse){
			this.reverse = reverse;
		}
		
		@Override
		public int compare( GFFEntry o1, GFFEntry o2 ) {
			int i = o1.seqid.compareTo( o2.seqid );
			if(i == 0){
				i = o1.start > o2.start ? 1 : (o1.start < o2.start ? -1 : 0);
			}
			if(i == 0){
				i = o1.end < o2.end ? 1 : (o1.end > o2.end ? -1 : 0);
			}
			if(i == 0){
				if(o1.parsedType == null){
					if(o2.parsedType == null){
						i = 0;
					}else{
						i = -1;
					}
				}else{
					if(o2.parsedType == null){
						i = 1;
					}else{
						i = o1.parsedType.compareType( o2.parsedType );
					}
				}
			}
			return reverse ? -i : i;
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
		
		public GFFList getSubList(GFFType... type){
			if(type.length == 0){ return null; }
			if(type.length == 1){ return byType.get( type[0] ); }
			ArrayList<GFFEntry> all = byType.get( type[0] ).entries;
			for(int i=1;i<type.length;i++){
				all.addAll( byType.get( type[i] ).entries );
			}
			return new GFFList( all );
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
				String score = Double.isNaN( en.getScore() ) ? "." : en.getScore()+"";
				set.add( new ResultSet(new Result[]{
				                                    new CategoricalResult( "SequenceID", "", en.getSeqid() ),
				                                    new CategoricalResult( "Source", "", en.getSource() ),
				                                    new CategoricalResult( "Type", "", en.getType() ),
				                                    new NumericalResult( "Start", "", en.getStart() ),
				                                    new NumericalResult( "End", "", en.getEnd() ),
				                                    new CategoricalResult( "Score", "", score ),
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
			stop_codon,
			other;
			
			public int compareType(GFFType t2){
				return this.ordinal() < t2.ordinal() ? -1 : (this.ordinal() > t2.ordinal() ? 1 : 0);
			}
			
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
		
		private void setID(String id){
			this.id = id;
			if(this.attributes == null){
				this.attributes = new HashMap<String, String>();
			}
			this.attributes.put( "ID", id );
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
		
		public void addChild(GFFEntry child){
			if(children == null){
				children = new LinkedList<GFFParser.GFFEntry>();
			}
			children.add( child );
		}
		
		public void addParent(GFFEntry parent){
			if(this.parent == null){
				this.parent = new LinkedList<GFFParser.GFFEntry>();
			}
			this.parent.add( parent );
			if(this.attributes == null){
				this.attributes = new HashMap<String, String>();
			}
			this.attributes.put( "Parent", parent.id );
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
