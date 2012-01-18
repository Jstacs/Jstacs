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

package projects.dispom;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.StringTokenizer;

import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.MotifAnnotation;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.AbstractStringExtractor;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.StringExtractor;
import de.jstacs.io.SymbolExtractor;
import de.jstacs.results.NumericalResult;
import de.jstacs.utils.ComparableElement;

/**
 * This class implements some parser for the results of de-novo motif
 * discoverers. The results are parsed and packed in instances of
 * {@link de.jstacs.data.sequences.annotation.SequenceAnnotation}.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class ParserToolBox {

	public static DataSet annotateWithMoAnResults( DataSet data, boolean addAnnot, File... files ) throws Exception{
		
		AbstractStringExtractor best = null;
		double bestScore = Double.NEGATIVE_INFINITY;
		
		String str = null;
		for(File f : files){
			AbstractStringExtractor curr = new SparseStringExtractor( f.getAbsolutePath(), '#' );
			while(curr.hasMoreElements() && !(str = curr.nextElement()).startsWith( "Score:" ));
			if(str != null && str.startsWith( "Score:" )){
				double score = Double.parseDouble( str.split( ":\\s+" )[1] );
				if(score > bestScore){
					bestScore = score;
					best = curr;
				}
			}
		}
		if(best != null){
			return annotateWithMoAnResults( data, addAnnot, best );
		}else{
			return null;
		}
	}
	
	public static DataSet annotateWithMoAnResults( DataSet data, boolean addAnnot, AbstractStringExtractor moAn ) throws Exception {

		Sequence[] moAnSeqs = new Sequence[data.getNumberOfElements()];
		
		for(int i=0;i<moAnSeqs.length;i++){
			moAnSeqs[i] = data.getElementAt( i );
		}
		String str = null;
		int num = 1;
		while(moAn.hasMoreElements()){
			while( moAn.hasMoreElements() && !(str = moAn.nextElement()).startsWith( "Hits" ) );
			
			moAn.nextElement();
			
			while(moAn.hasMoreElements() && !(str = moAn.nextElement()).startsWith( "---" )){
				str = str.trim();
				String[] split = str.split( ":?\\s+" );
				int idx = Integer.parseInt( split[0] );
				int pos = Integer.parseInt( split[1] );
				String subseq = split[2];
				double score = Double.parseDouble( split[3] );
				//TODO
				Strand strand = Strand.UNKNOWN;
				moAnSeqs[idx] = moAnSeqs[idx].annotate( addAnnot, new MotifAnnotation( "motif"+num, pos, subseq.length(), strand, new NumericalResult("score","score of element",score) ) );
			}
			num++;
		}
		
		return new DataSet("Sample annotated by MoAn",moAnSeqs);
		
	}
	
	public static DataSet annotateBestCoBindResults(DataSet data, String cobindFileName) throws Exception {
		StringExtractor cobind = new StringExtractor( new File(cobindFileName), 100,'&' );
		
		double best = Double.NEGATIVE_INFINITY;
		DataSet bestS = null;
		while(cobind.hasMoreElements()){
			ComparableElement<DataSet, Double> el = annotateSingleCoBindResults( data, cobind );
			if(el.getWeight() > best){
				best = el.getWeight();
				bestS = el.getElement();
			}
		}
		
		return bestS;		
	}
	
	public static ComparableElement<DataSet, Double> annotateSingleCoBindResults(DataSet data, StringExtractor cobind) throws EmptyDataSetException, WrongAlphabetException{
		
		Sequence[] cobindSeqs = new Sequence[data.getNumberOfElements()];
		
		for(int i=0;i<cobindSeqs.length;i++){
			cobindSeqs[i] = data.getElementAt( i );
		}
		
		String str = null;
		
	//	while(cobind.hasMoreElements() && !(str = cobind.nextElement()).startsWith( "START OF RUN:" ));
		
		int num = 1;
		double annScore = Double.NEGATIVE_INFINITY;
		outerloop:
		while(true){
			while(cobind.hasMoreElements() && !(str = cobind.nextElement()).startsWith( "# BEST_SITES" )){
				//System.out.println(str);
				if(str.startsWith( "# >> Maximum conbined energy for both units (MaxUU):" )){
					annScore = Double.parseDouble( str.substring( "# >> Maximum conbined energy for both units (MaxUU):".length() ).trim() );
					cobind.nextElement();
					cobind.nextElement();
					break outerloop;
				}
			}
			//System.out.println(str);
			if(cobind.hasMoreElements()){
				cobind.nextElement();
			}else{
				break outerloop;
			}
			while(cobind.hasMoreElements()){
				str = cobind.nextElement();
				if(str.startsWith( "# ave_log_S" )){
					num++;
					break;
				}
				String[] sp = str.split( "[\\s:]+" );
				int idx = Integer.parseInt( sp[2] )-1;
				double score = Double.parseDouble( sp[3] );
				int pos = 0;
				Strand strand = Strand.FORWARD;
				sp[4] = sp[4].trim();
				if(sp[4].endsWith( "'" )){
					pos = Integer.parseInt( sp[4].substring( 0, sp[4].length()-1 ) )-1;
					strand = Strand.REVERSE;
				}else{
					pos = Integer.parseInt( sp[4] );
				}
				cobindSeqs[idx] = cobindSeqs[idx].annotate( true, new MotifAnnotation( "motif"+num, pos, sp[5].length(), strand, new NumericalResult("score","score of element",score) ) );	
			}
		}
		return new ComparableElement<DataSet, Double>( new DataSet("Sample annotated by Co-Bind",cobindSeqs), annScore );
	}
	
	public static DataSet annotateWithCisModuleResults( DataSet data, boolean addAnnot, String cismodFileName ) throws Exception {
		StringExtractor cismod = new StringExtractor( new File(cismodFileName), 100 );
		
		Sequence[] cismodSeqs = new Sequence[data.getNumberOfElements()];
		
		for(int i=0;i<cismodSeqs.length;i++){
			cismodSeqs[i] = data.getElementAt( i );
		}
		
		String str = null;
		
		while(cismod.hasMoreElements() && !(str = cismod.nextElement()).startsWith( "Motif" ));
		
		int num = 1;
		
		while(true) {
			
			int w = 0;
			while(cismod.hasMoreElements() && !(str = cismod.nextElement()).startsWith( ">" )){
				if(str.startsWith( (w+1)+"" )){
					w++;
				}
			}
			if(!cismod.hasMoreElements()){
				break;
			}
			
			do{
				String[] sp = str.split( "\\s+" );
				//TODO
				int idx = Integer.parseInt( sp[0].substring( 4 ) );
				Strand strand = Strand.FORWARD;
				if(sp[1].trim().equals( "r" )){
					strand = Strand.REVERSE;
				}
				int pos = Integer.parseInt( sp[2] )-1;
				double score = Double.parseDouble( sp[3] );
				
				cismodSeqs[idx] = cismodSeqs[idx].annotate( true, new MotifAnnotation( "motif"+num, pos, w, strand, new NumericalResult("score","score of element",score) ) );
				
				cismod.nextElement();
				
			}while( cismod.hasMoreElements() && (str = cismod.nextElement()).startsWith( ">" ) );
			
			num++;
			if(!cismod.hasMoreElements()){
				break;
			}
		}
		
		return new DataSet( "Sample annotated by CisModule", cismodSeqs );
		
	}
	
	public static DataSet annotateWithGibbsSamplerResults( DataSet data, boolean addAnnot, String gibbsFileName ) throws FileNotFoundException, IOException, IllegalArgumentException, EmptyDataSetException, WrongAlphabetException{
		
		StringExtractor gibbs = new StringExtractor(new File(gibbsFileName),100);
		
		Sequence[] gibbsSeqs = new Sequence[data.getNumberOfElements()];
		
		int idx = 0;
		String line = gibbs.getElement( idx );
		while( !line.matches( "\\s+[0-9]+,\\s+[0-9]+\\s+[0-9]+\\s+[acgt]*\\s*[ACGT]+\\s*[acgt]*\\s+[0-9]+.*" ) ){
			idx++;
			if(idx < gibbs.getNumberOfElements()){
				line = gibbs.getElement( idx );
			}else{
				return data;
			}
		}
		while( line.matches( "\\s+[0-9]+,\\s+[0-9]+\\s+[0-9]+\\s+[acgt]*\\s*[ACGT]+\\s*[acgt]*\\s+[0-9]+.*" ) ){
			String[] temp = line.split( "[\\s,]+" );
			int seq = Integer.parseInt( temp[1] );
			int left = Integer.parseInt( temp[3] );
			int i=4;
			while(temp[i].matches( "[acgtACGT]+" )){
				i++;
			}
			int right = Integer.parseInt( temp[i] );
			Boolean fwd = null;
			if(temp.length >= i+3){
				fwd = temp[i+2].equals( "F" );
			}
			double val = Double.parseDouble( temp[i+1] );
			
			MotifAnnotation ann = new MotifAnnotation("motif",Math.min( left,right )-1,Math.max( left,right)-Math.min(left,right)+1,fwd==null ? Strand.UNKNOWN : (fwd ? Strand.FORWARD : Strand.REVERSE), new NumericalResult("score","probability of element",val));
			if(gibbsSeqs[seq-1] == null){
				gibbsSeqs[seq-1] = data.getElementAt( seq-1 ).annotate( addAnnot, ann );
			}else{
				gibbsSeqs[seq-1] = gibbsSeqs[seq-1].annotate( true, ann );
			}
			idx++;
			line = gibbs.getElement( idx );
		}
		for(int i=0;i<gibbsSeqs.length;i++){
			if(gibbsSeqs[i] == null){
				gibbsSeqs[i] = data.getElementAt( i );
			}
		}
		return new DataSet("gibbsSampler",gibbsSeqs);
	}
	
	
	public static DataSet annotateWithMemeResults( DataSet data, boolean addAnnot, String memeFileName ) throws IllegalArgumentException, EmptyDataSetException, FileNotFoundException, IOException, WrongAlphabetException{
		StringExtractor meme = new StringExtractor(new File(memeFileName),100);
		
		Sequence[] memeSeqs = new Sequence[data.getNumberOfElements()];
		int start = -1;
		int idx = -1;
		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		LinkedList<MotifAnnotation> annots = new LinkedList<MotifAnnotation>();
		for(int i=0;i<meme.getNumberOfElements();i++){
			String line = meme.getElement( i );
			if(line.matches( "MOTIF\\s+\\d\\s+width.*" )){
				int w = line.indexOf( "width" );
				Integer num = Integer.parseInt( line.substring( line.indexOf( "MOTIF" )+"MOTIF".length(), w ).trim() );
				int length2 = Integer.parseInt( line.substring( line.indexOf( "=" )+1, line.indexOf( "sites" ) ).trim() );
				map.put( num, length2 );
			}
			if(meme.getElement( i ).indexOf( "Combined block diagrams" ) > -1){
				start = 3;
				idx = 0;
			}else if(start > 0){
				start--;
			}else if(start == 0){
				if(meme.getElement( i ).indexOf( "---" ) > -1){
					start = -1;
				}else{
					String s = meme.getElement( i );
					String[] temp = s.split( "\\s+" );
					s = temp[2];
					
					temp = s.split( "_" );
					if(temp.length > 1){
						annots.clear();
						int off = 0;
						for(int j=0;j<temp.length;j++){
							if(temp[j].matches( "\\[[+-]?\\d\\(.+\\)\\]" )){
								Integer num = null;
								Strand strand = Strand.FORWARD;
								if(temp[j].indexOf( "[+" ) > -1){
									num = Integer.parseInt( temp[j].substring( temp[j].indexOf( "[" )+2, temp[j].indexOf( "(" ) ) );
								}else if(temp[j].indexOf( "[-" ) > -1){
									num = Integer.parseInt( temp[j].substring( temp[j].indexOf( "[" )+2, temp[j].indexOf( "(" ) ) );
									strand = Strand.REVERSE;
								}else{
									num = Integer.parseInt( temp[j].substring( temp[j].indexOf( "[" )+1, temp[j].indexOf( "(" ) ) );
								}
								double val = Double.parseDouble( temp[j].substring( temp[j].indexOf( "(" )+1,temp[j].indexOf( ")" ) - 1 ) );
								annots.add( new MotifAnnotation("Motif "+num,off,map.get( num ),strand, new NumericalResult("score","1-(p-value)",1-val)) );
							}else{
								off += Integer.parseInt( temp[j] );
							}
						}
						memeSeqs[idx] = data.getElementAt( idx ).annotate( false, /*new CisRegulatoryModuleAnnotation("CRM1",*/annots.toArray( new MotifAnnotation[0] )/*)*/ );

					}else{
						memeSeqs[idx] = data.getElementAt( idx );
					}
					idx++;
				}
			}
		}
		
		DataSet memeSample = new DataSet("meme",memeSeqs);
		
		return memeSample;
	}
	
	public static DataSet annotateWithDemeResults( DataSet data, boolean addAnnot, String demeFileName ) throws IllegalArgumentException, EmptyDataSetException, FileNotFoundException, IOException, WrongAlphabetException{
		StringExtractor deme = new StringExtractor(new File(demeFileName),100,'$');
		
		Sequence[] demeSeqs = new Sequence[data.getNumberOfElements()];
		
		
		int start = -1;
		int idx = -1;
		for(int i=0;i<deme.getNumberOfElements();i++){
			String line = deme.getElement( i );
			if(line.indexOf( "# Training sequences" ) > -1){

				start = 1;
				idx = 0;
			}else if(start > 0){
				start--;
			}else if(start == 0){
				if(line.length() == 0){
					start = -1;
				}else{
					String[] temp = line.split( "\\s+" );
					if(temp[3].trim().equals( "1" )){
						int pos = Integer.parseInt( temp[4] )-1;
						Strand strand = Strand.FORWARD;
						if(temp[5].trim().equals( "-" )){
							strand = Strand.REVERSE;
						}
						int length2 = temp[6].trim().length();
						double val = Double.parseDouble( temp[8].trim() );
						//demeSeqs[idx] = data.getElementAt( idx ).annotate( addAnnot, new CisRegulatoryModuleAnnotation("CRM1",new MotifAnnotation[]{new MotifAnnotation("motif1",pos,length2,strand)}) );
						demeSeqs[idx] = data.getElementAt( idx ).annotate( addAnnot, new MotifAnnotation("motif1",pos,length2,strand, new NumericalResult("score","P(C|X,theta)",val)) );
						idx++;
					}
				}
			}
		}
		
		DataSet demeSample = new DataSet("deme",demeSeqs);
		
		return demeSample;
	}
	
	
	public static DataSet annotateWithImprobizerResults( DataSet data, boolean addAnnot, String improbFileName, int improbLength) throws IllegalArgumentException, EmptyDataSetException, FileNotFoundException, IOException, WrongAlphabetException{
		StringExtractor imp = new StringExtractor(new File(improbFileName),100);
		
		Sequence[] impSeqs = new Sequence[data.getNumberOfElements()];
		LinkedList<MotifAnnotation> annots = new LinkedList<MotifAnnotation>();
		int lastseq = -1;
		for(int i=0;i<imp.getNumberOfElements();i++){

			String line = imp.getElement( i );
			String[] parts = line.split( "\\s+" );
			int seq = Integer.parseInt( parts[2].substring( 3 ) );
			double val = Double.parseDouble( parts[1] );
			boolean fwd = parts[0].equals( "1" );
			int pos = Integer.parseInt( parts[3] );
			MotifAnnotation ann = new MotifAnnotation("Motif",pos,improbLength, fwd ? Strand.FORWARD : Strand.REVERSE, new NumericalResult("score","score",val));

			if(seq != lastseq){
				
				LinkedList<MotifAnnotation> toRemove = new LinkedList<MotifAnnotation>();
				for(int j=0;j<annots.size();j++){
					MotifAnnotation a = annots.get( j );
					for(int k=j+1;k<annots.size();k++){
						MotifAnnotation b = annots.get( k );
						if(a.getPosition() == b.getPosition() && a.getLength() == b.getLength()){
							toRemove.add( b );
						}
					}
				}
				annots.removeAll( toRemove );
				if(lastseq>=0){
					impSeqs[lastseq] = data.getElementAt( lastseq ).annotate( false, annots.toArray( new MotifAnnotation[0] ) );
				}
				annots.clear();
				if(seq != lastseq){
					annots.add( ann );
				}
			}else{
				annots.add( ann );
			}
			if(i==imp.getNumberOfElements()-1){
				LinkedList<MotifAnnotation> toRemove = new LinkedList<MotifAnnotation>();
				for(int j=0;j<annots.size();j++){
					MotifAnnotation a = annots.get( j );
					for(int k=j+1;k<annots.size();k++){
						MotifAnnotation b = annots.get( k );
						if(a.getPosition() == b.getPosition() && a.getLength() == b.getLength()){
							toRemove.add( b );
						}
					}
				}
				annots.removeAll( toRemove );
				impSeqs[seq] = data.getElementAt( seq ).annotate( false, annots.toArray( new MotifAnnotation[0] ) );
			}
			lastseq = seq;
		}
		
		for(int i=0;i<impSeqs.length;i++){
			if(impSeqs[i] == null){
				impSeqs[i] = data.getElementAt( i );
			}
		}
		
		DataSet impSample = new DataSet("improbizer",impSeqs);
		
		return impSample;
	}
	
	public static int getImprobizerMotifLength( String fName ) throws IOException {
		BufferedReader r = new BufferedReader( new FileReader( new File( fName ) ) );
		String line = r.readLine().trim();
		r.close();
		int idx = line.lastIndexOf( ' ' );
		return line.substring( idx+1 ).length();
	}
	
	
	//TODO skip first lines
	public static DataSet annotateWithWeederResults( DataSet data, boolean addAnnot, String weederFileName) throws FileNotFoundException, IOException, IllegalArgumentException, EmptyDataSetException, WrongAlphabetException{
		StringExtractor weed = new StringExtractor(new File(weederFileName),100);
		
		int i=0;
		while( !weed.getElement( i ).startsWith( "Best occurrences:" ) ){
			i++;
		}
		i+=2;
		
		Sequence[] weedSeqs = new Sequence[data.getNumberOfElements()];
		int start = -1;
		int idx = -1;
		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		LinkedList<MotifAnnotation> annots = new LinkedList<MotifAnnotation>();
		int lastseq = -1;
		for(;i<weed.getNumberOfElements();i++){

			String line = weed.getElement( i );
			String[] parts = line.split( "\\s+" );
			int off = parts[0].length() == 0 ? 1 : 0;
			int seq = Integer.parseInt( parts[off] );
			boolean fwd = parts[off+1].equals( "+" );
			String str = parts[off+2];
			int length = str.length();
			if(str.startsWith( "[" )){
				length--;
			}
			if(str.endsWith( "]" )){
				length--;
			}
			int pos = Integer.parseInt( parts[off+3] );
			
			double val = Double.parseDouble( parts[off+4].substring( 1, parts[off+4].length()-1 ) );
			MotifAnnotation ann = new MotifAnnotation("Motif",pos-1,length, fwd ? Strand.FORWARD : Strand.REVERSE, new NumericalResult("score","score",val));
			if(seq != lastseq){
				
				LinkedList<MotifAnnotation> toRemove = new LinkedList<MotifAnnotation>();
				for(int j=0;j<annots.size();j++){
					MotifAnnotation a = annots.get( j );
					for(int k=j+1;k<annots.size();k++){
						MotifAnnotation b = annots.get( k );
						if(a.getPosition() == b.getPosition() && a.getLength() == b.getLength()){
							toRemove.add( b );
						}
					}
				}
				annots.removeAll( toRemove );
				if(lastseq>0){
					weedSeqs[lastseq-1] = data.getElementAt( lastseq-1 ).annotate( addAnnot, annots.toArray( new MotifAnnotation[0] ) );
				}
				annots.clear();
				if(seq != lastseq){
					annots.add( ann );
				}
			}else{
				annots.add( ann );
			}
			if(i == weed.getNumberOfElements()-1 || weed.getElement( i+1 ).length()==0){
				LinkedList<MotifAnnotation> toRemove = new LinkedList<MotifAnnotation>();
				for(int j=0;j<annots.size();j++){
					MotifAnnotation a = annots.get( j );
					for(int k=j+1;k<annots.size();k++){
						MotifAnnotation b = annots.get( k );
						if(a.getPosition() == b.getPosition() && a.getLength() == b.getLength()){
							toRemove.add( b );
						}
					}
				}
				annots.removeAll( toRemove );
				weedSeqs[seq-1] = data.getElementAt( seq-1 ).annotate( addAnnot, annots.toArray( new MotifAnnotation[0] ) );
				break;
			}
			lastseq = seq;
		}
		
		for(i=0;i<weedSeqs.length;i++){
			if(weedSeqs[i] == null){
				weedSeqs[i] = data.getElementAt( i );
			}
		}
		
		DataSet weederSample = new DataSet("weeder",weedSeqs);
		
		return weederSample;
	}
	
	private static final String AGLAM_RESULTS = "! Best possible alignment:";

	private static final String AGLAM_ANCHOR = "Anchors are given as:";

	public static DataSet annotateWithAGlamResults( DataSet data, boolean addAnnot, String ouputFileName, String motifName, boolean evalues ) throws Exception {
		Sequence[] seqs = data.getAllElements();
		boolean[] add = new boolean[seqs.length];
		Arrays.fill( add, addAnnot );
		BufferedReader r = new BufferedReader( new FileReader( ouputFileName ) );
		String line, annot, anchor = null;
		StringTokenizer st;
		while( ( line = r.readLine() ) != null && !line.startsWith( AGLAM_RESULTS ) ) {
			if( line.startsWith( AGLAM_ANCHOR ) ) {
				anchor = line.substring( AGLAM_ANCHOR.length() );
			}
		}
		int[] anchors = new int[data.getNumberOfElements()];
		if( anchor == null ) {
			Arrays.fill( anchors, 0 );
		} else {
			st = new StringTokenizer( anchor );
			for( int i = 0; i < anchors.length; i++ ) {
				anchors[i] = Integer.parseInt( st.nextToken() );
			}
		}

		char seqComment = '>';
		int idx = 0, v1, v2, h;
		MotifAnnotation ma;
		Strand s;
		while( ( line = r.readLine() ) != null ) {
			if( line.length() > 0 && line.charAt( 0 ) == seqComment ) {
				annot = r.readLine().trim();
				st = new StringTokenizer( annot, " " );
				if( st.countTokens() > 1 ) {
					//System.out.println( line + "\n" + annot );

					v1 = anchors[idx] + Integer.parseInt( st.nextToken() );
					st.nextToken();
					v2 = anchors[idx] + Integer.parseInt( st.nextToken() );

					s = v1 < v2 ? Strand.FORWARD : Strand.REVERSE;
					if( s == Strand.REVERSE ) {
						h = v1;
						v1 = v2 - 1;
						v2 = h;
					} else {
						v1--;
					}
					
					NumericalResult res = null;

					st.nextToken();
					if(evalues){
						st.nextToken();
					}
					String temp = st.nextToken();
					temp = temp.substring( 1, temp.length()-1 );
					double val = Double.parseDouble( temp );
					if(evalues){
						res = new NumericalResult("score","negative E-value",-val);
					}else{
						res = new NumericalResult("score","alignment score",val);
					}

					//System.out.println( seqs[idx].toString().substring( v1, v2 ) );

					ma = new MotifAnnotation( motifName, v1, v2 - v1, s, res );//TODO more information?
					seqs[idx] = seqs[idx].annotate( add[idx], ma );
					add[idx] = true;
				}
				idx++;
			}
		}
		r.close();
		return new DataSet( "data with annotation from A-GLAM", seqs );
	}

	private static final String DME_MOTIF = "AC";

	private static final String DME_BINDINGSITE = "BS";

	public static DataSet annotateWithDMEResults( DataSet data, boolean addAnnot, String ouputFileName, int index ) throws Exception {
		Sequence[] seqs = data.getAllElements();
		boolean[] add = new boolean[seqs.length];
		Arrays.fill( add, addAnnot );
		BufferedReader r = new BufferedReader( new FileReader( ouputFileName ) );
		String line, motifName, seq;
		int i = 0;
		do {
			while( ( line = r.readLine() ) != null && !line.startsWith( DME_MOTIF ) );
			i++;
		} while( i < index && line != null );

		if( line == null ) {
			return null;
		} else {
			motifName = "motif " + index + " (" + line.substring( 2 ).trim() + ")";
			while( ( line = r.readLine() ) != null && !line.startsWith( DME_BINDINGSITE ) );

			char sep = ';';
			int idx = 0, v1, v2, h;
			MotifAnnotation ma;
			Strand s;
			do {
				//System.out.println( line );

				idx = line.indexOf( sep ) + 1;
				h = line.indexOf( sep, idx );
				seq = line.substring( idx, h ).trim();
				idx = Integer.parseInt( seq.replace( "seq", "" ) );//TODO
				v1 = h + 1;
				h = line.indexOf( sep, v1 );
				v1 = Integer.parseInt( line.substring( v1, h ).trim() );
				v2 = h + 1;
				h = line.indexOf( sep, v2 );
				v2 = Integer.parseInt( line.substring( v2, h ).trim() );
				h = line.indexOf( sep, h + 1 ) + 1;
				s = line.substring( h ).trim().startsWith( "p" ) ? Strand.FORWARD : Strand.REVERSE;
				h = line.indexOf( sep, h+1 );
				double val = Double.parseDouble( line.substring( h+1 ) );
				//System.out.println( seqs[idx].toString().substring( v1, v1+v2 ) + " " + s );

				ma = new MotifAnnotation( motifName, v1, v2, s, new NumericalResult("score","score",val) );//TODO more information?
				seqs[idx] = seqs[idx].annotate( add[idx], ma );
				add[idx] = true;
			} while( ( line = r.readLine() ) != null && line.startsWith( DME_BINDINGSITE ) );
			r.close();
			return new DataSet( "data with annotation from DME", seqs );
		}
	}
	
	public static DataSet annotateWithAmadeusResults( DataSet data, boolean addAnnot, String targetFile, int motifIndex, int length ) throws Exception {
		Sequence[] seqs = data.getAllElements();
		boolean[] add = new boolean[seqs.length];
		Arrays.fill( add, addAnnot );
		BufferedReader r = new BufferedReader( new FileReader( targetFile ) );
		String line, motif, set, pos, p, m = "" + motifIndex, motifName = "motif " + motifIndex;
		char s;
		Strand strand;
		SymbolExtractor se = new SymbolExtractor( "\t" );
		MyList<String> list = new MyList<String>();
		int j, i, l, seqIdx;
		while( ( line = r.readLine() ) != null ) {
			se.setStringToBeParsed( line );
			seqIdx = Integer.parseInt( se.nextElement() );//TODO
			motif = se.nextElement();
			if( motif.equals( m ) ) {
				se.nextElement();
				se.nextElement();
				se.nextElement();
				set = se.nextElement();
				if( set.equalsIgnoreCase( "Target" ) ) {
					//System.out.print( line );
					list.clear();
					while( se.hasMoreElements() ) {
						list.add( se.nextElement() );
					}
					list.sort();
					l = list.size()-1;
					for( i = 0; i <= l; i++ ) {
						pos = list.get( i );
						j = pos.indexOf( '(' );
						s = pos.charAt( j + 1 );
						pos = pos.substring( 0, j );
						strand = null;
						//System.out.println( pos + " " + s );
						while( i < l && list.get(1+i).startsWith(pos) ) {
							p = list.get( 1+i++ );
							//System.out.println( s + " vs " + p.charAt( p.indexOf( '(' ) + 1 ) + "\t" + (s != p.charAt( p.indexOf( '(' ) + 1 )) );
							if( s != p.charAt( p.indexOf( '(' ) + 1 ) ) {
								strand = Strand.UNKNOWN;
							}
						}
						if( strand == null ) {
							strand = s == '+' ? Strand.FORWARD : Strand.REVERSE;
						}
						//Sequence xxx = seqs[seqIdx].getSubSequence( seqs[seqIdx].getLength()-1+Integer.parseInt(pos), length );
						//System.out.print( "\t-> " + xxx + "\t" + xxx.reverseComplement() +"\t" +  strand );
						seqs[seqIdx] = seqs[seqIdx].annotate( add[seqIdx], new MotifAnnotation( motifName, seqs[seqIdx].getLength()-1+Integer.parseInt(pos), length, strand ) );//TODO more information?
						add[seqIdx] = true;
					}
					//System.out.println();
				}
			}
		}

		return new DataSet( "data with annotation from Amadeus", seqs );
	}
}


class MyList<T extends Comparable> {

	private Comparable[] comp;
	private int size;
	
	public MyList(){
		comp = new Comparable[10];
		size = 0;
	}
	
	public void add( T item ) {
		if( comp.length == size ){
			//extend
			Comparable[] help = new Comparable[comp.length*2];
			System.arraycopy(comp, 0, help, 0, size);
			comp = help;
		}
		comp[size++] = item;
	}
	
	public T get( int index ) throws IndexOutOfBoundsException {
		if( index < size ){
			return (T) comp[index];
		} else {
			throw new IndexOutOfBoundsException();
		}
	}
	
	public int size() {
		return size;
	}
	
	public void sort(){
		Arrays.sort( comp, 0, size );
	}
	
	public void clear() {
		size = 0;
	}
}