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
package de.jstacs.data.sequences.annotation;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.Sequence;


/**
 * This class implements an {@link SequenceAnnotationParser} that parses a {@link ReferenceSequenceAnnotation} from the comment lines of a sequence.
 * 
 * @author Jan Grau
 */
public class ReferenceSequenceAnnotationParser extends SplitSequenceAnnotationParser {

	private String key;
	private AlphabetContainer alphabet;
	private String delim;
	
	/**
	 * Creates a new {@link ReferenceSequenceAnnotationParser} with user-specified delimiters.
	 * 
	 * @param key the key for the {@link ReferenceSequenceAnnotation}
	 * @param alphabet the {@link AlphabetContainer} of the reference sequence
	 * @param keyValueDelimiter the delimiter between key and corresponding value
	 * @param annotationDelimiter the delimiter between different {@link SequenceAnnotation}s
	 * @param delim the delimiter between symbols a the reference sequence
	 * 
	 * @throws IllegalArgumentException if the delimiters are identical
	 */
	public ReferenceSequenceAnnotationParser( String key, AlphabetContainer alphabet, String keyValueDelimiter, String annotationDelimiter, String delim ) throws IllegalArgumentException {
		super( keyValueDelimiter, annotationDelimiter );
		this.key = key;
		this.alphabet = alphabet;
		this.delim = delim;
	}
	
	/**
	 * Creates a new {@link ReferenceSequenceAnnotationParser} with user-specified delimiters.
	 * 
	 * @param key the key for the {@link ReferenceSequenceAnnotation}
	 * @param alphabet the {@link AlphabetContainer} of the reference sequence
	 * @param keyValueDelimiter the delimiter between key and corresponding value
	 * @param annotationDelimiter the delimiter between different {@link SequenceAnnotation}s
	 * 
	 * @throws IllegalArgumentException if the delimiters are identical
	 */
	public ReferenceSequenceAnnotationParser( String key, AlphabetContainer alphabet, String keyValueDelimiter, String annotationDelimiter ) throws IllegalArgumentException {
		this( key, alphabet, keyValueDelimiter, annotationDelimiter, alphabet.getDelim() );
	}

	protected void add( String type, String identifier ) {
		if( type.equalsIgnoreCase( this.key ) ) {
			try{
				annot.add( new ReferenceSequenceAnnotation( key, Sequence.create( alphabet, identifier, delim ) ) );
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
				res.append( annotations[0].getIdentifier() + keyValueDelimiter + ((ReferenceSequenceAnnotation)annotations[0]).getReferenceSequence() );
			}else{
				res.append( annotations[0].getType() + keyValueDelimiter + annotations[0].getIdentifier() );
			}
			for( int i = 1; i < annotations.length; i++ ) {
				if(annotations[i] instanceof ReferenceSequenceAnnotation){
					res.append( annotationDelimiter + " " + annotations[i].getIdentifier() + keyValueDelimiter + ((ReferenceSequenceAnnotation)annotations[i]).getReferenceSequence() );
				}else{
					res.append( annotationDelimiter + " " + annotations[i].getType() + keyValueDelimiter + annotations[i].getIdentifier() );
				}
			}
		}
		return res.toString();
	}

}
