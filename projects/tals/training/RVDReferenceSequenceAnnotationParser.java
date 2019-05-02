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
package projects.tals.training;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import projects.tals.RVDSequence;


/**
 * This class implements an {@link SequenceAnnotationParser} that parses a {@link ReferenceSequenceAnnotation} from the comment lines of a sequence.
 * 
 * @author Jan Grau
 */
public class RVDReferenceSequenceAnnotationParser extends SplitSequenceAnnotationParser {

	private String key;
	private AlphabetContainer alphabet12, alphabet13;
	
	/**
	 * Creates a new {@link RVDReferenceSequenceAnnotationParser} with user-specified delimiters.
	 * 
	 * @param key the key for the {@link ReferenceSequenceAnnotation}
	 * @param keyValueDelimiter the delimiter between key and corresponding value
	 * @param annotationDelimiter the delimiter between different {@link SequenceAnnotation}s
	 * @param delim the delimiter between symbols a the reference sequence
	 * 
	 * @throws IllegalArgumentException if the delimiters are identical
	 */
	public RVDReferenceSequenceAnnotationParser( String key, String keyValueDelimiter, String annotationDelimiter, AlphabetContainer alphabet12, AlphabetContainer alphabet13 ) throws IllegalArgumentException {
		super( keyValueDelimiter, annotationDelimiter );
		this.key = key;
		this.alphabet12 = alphabet12;
		this.alphabet13 = alphabet13;
	}

	
	protected ReferenceSequenceAnnotation getSequenceAnnotation(String seqString) throws IllegalArgumentException, WrongAlphabetException, WrongSequenceTypeException, DoubleSymbolException{
		return new ReferenceSequenceAnnotation( key, new RVDSequence(alphabet12,alphabet13,seqString) );
	}	
	
	protected void add( String type, String identifier ) {
		if( type.equalsIgnoreCase( this.key ) ) {
			try{
				annot.add( getSequenceAnnotation(identifier) );
			}catch(Exception e){
				RuntimeException re = new RuntimeException( e.getMessage() );
				re.setStackTrace( e.getStackTrace() );
				throw re;
			}
		} else {
			super.add( type, identifier );
		}
	}
	
	@Override
	public String parseAnnotationToComment( char commentChar, SequenceAnnotation... annotations ) {
		StringBuffer res = new StringBuffer();
		res.append(  commentChar );
		if( annotations != null && annotations.length > 0 ) {
			if(annotations[0] instanceof ReferenceSequenceAnnotation){
				Sequence rs = ((ReferenceSequenceAnnotation)annotations[0]).getReferenceSequence();
				res.append( annotations[0].getIdentifier() + keyValueDelimiter + rs.toString( rs.getAlphabetContainer().getDelim(), 0, rs.getLength() ) );
			}else{
				res.append( annotations[0].getType() + keyValueDelimiter + annotations[0].getIdentifier() );
			}
			for( int i = 1; i < annotations.length; i++ ) {
				if(annotations[i] instanceof ReferenceSequenceAnnotation){
					Sequence rs = ((ReferenceSequenceAnnotation)annotations[i]).getReferenceSequence();
					res.append( annotationDelimiter + " " + annotations[i].getIdentifier() + keyValueDelimiter + rs.toString( rs.getAlphabetContainer().getDelim(), 0, rs.getLength() ) );
				}else{
					res.append( annotationDelimiter + " " + annotations[i].getType() + keyValueDelimiter + annotations[i].getIdentifier() );
				}
			}
		}
		return res.toString();
	}

}
