package de.jstacs.utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;

/**
 * Class for reading large DNA sequences (e.g. genomes) from disk sequentially.
 * @author Jan Grau
 *
 */
public class LargeSequenceReader {

	/**
	 * Returns the next chunk of input sequences. Input sequences are split at non-ACGT characters and parts between those are returned independently.
	 * For longer sequences that are split in this manner, the offsets of the individual chunks within the complete sequence are also returned.
	 * @param read the source of the sequence data
	 * @param lastHeader the header of previous chunk, may be empty initially and is modified internally
	 * @param minimumLength the minimum length of sequences or chunks considered, shorter sequences are skipped
	 * @return a {@link Pair} with the offsets of the individual chunks and the sequences
	 * @throws IOException if the sequence could not be read from <code>read</code>
	 */
	public static Pair<IntList,ArrayList<Sequence>> readNextSequences(BufferedReader read, StringBuffer lastHeader,  int minimumLength ) throws IOException {
		//System.out.println("started reading");
		String str = null;
		
		StringBuffer line = new StringBuffer();
		
		IntList starts = new IntList();
		ArrayList<Sequence> seqs = new ArrayList<>();
		
		Pattern acgt = Pattern.compile( "[ACGT]+", Pattern.CASE_INSENSITIVE );
		
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		
		int size = 0;
		
		while( (str = read.readLine()) != null || line.length() > 0 ){
			if(str != null){
				str = str.trim();
			}
			if(str == null || str.startsWith( ">" )){//next sequence
				String header = lastHeader.toString();
				if(str != null){
					lastHeader.delete( 0, lastHeader.length() );
					lastHeader.append( str.substring( 1 ).trim() );
				}
				if(line.length() > 0){//we have a sequence
					int idx = header.indexOf(" ");
					if( idx > 0 ) {//TODO new
						header = header.substring(0, idx);
					}
					SequenceAnnotation annotation = new SequenceAnnotation( "id", header );
					
					String seqStr = line.toString();
					line.delete( 0, line.length() );
					Matcher match = acgt.matcher( seqStr );
					while(match.find()){
						int start = match.start();
						int end = match.end();
						int l = end-start;
						if( l >= minimumLength ) {
							try{
								Sequence seq = Sequence.create( con, seqStr.substring( start, end ) );
								seq = seq.annotate( false, annotation );
								seqs.add( seq );
							}catch(WrongAlphabetException ex){
								throw new RuntimeException(ex);
							}
							size += l;
							starts.add( start );
							//	ends.add( seqStr.length()-end );
						}
					}
					if(size > 1E7 || str == null){
						return new Pair<IntList, ArrayList<Sequence>>(starts, seqs);
					}
				}
			}else{
				line.append( str );
			}	
		}
		return null;
	}

	
}
