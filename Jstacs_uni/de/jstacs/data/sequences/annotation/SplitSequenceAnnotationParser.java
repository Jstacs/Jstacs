/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.data.sequences.annotation;

import java.util.LinkedList;

/**
 * This class implements a simple {@link SequenceAnnotationParser} which simply splits the comments by specific delimiters.
 * 
 * @author Jens Keilwagen
 */
public class SplitSequenceAnnotationParser implements SequenceAnnotationParser {

	/**
	 * The delimiter between key and value
	 */
	protected String keyValueDelimiter;
	
	/**
	 * The delimiter between different annotations
	 */
	protected String annotationDelimiter;
	/**
	 * The internal list of current {@link SequenceAnnotation}s.
	 */
	protected LinkedList<SequenceAnnotation> annot;
	private static final SequenceAnnotation[] empty = new SequenceAnnotation[0];	
	
	/**
	 * Creates a new {@link SplitSequenceAnnotationParser} with specific delimiters, i.e., key value
	 * delimiter &quot;=&quot; and annotation delimiter &quot;;&quot;.
	 *  
	 * @see SplitSequenceAnnotationParser#SplitSequenceAnnotationParser(String, String)
	 */
	public SplitSequenceAnnotationParser() {
		this( "=", ";" );
	}
	
	/**
	 * Creates a new {@link SplitSequenceAnnotationParser} with user-specified delimiters.
	 * 
	 * @param keyValueDelimiter the delimiter between key and corresponding value
	 * @param annotationDelimiter the delimiter between different {@link SequenceAnnotation}s
	 * 
	 * @throws IllegalArgumentException if the delimiters are identical
	 */
	public SplitSequenceAnnotationParser( String keyValueDelimiter, String annotationDelimiter ) {
		if( keyValueDelimiter.equals( annotationDelimiter ) ) {
			throw new IllegalArgumentException( "The delimiters have to be different." );
		}
		this.annotationDelimiter = annotationDelimiter;
		this.keyValueDelimiter = keyValueDelimiter;
		annot = new LinkedList<SequenceAnnotation>();
	}
	
	public void addToAnnotation( String unparsed ) {
		String[] split = unparsed.substring(1).split( annotationDelimiter );
		for( String current : split ) {
			int idx = current.indexOf( keyValueDelimiter );
			if( idx >= 0 ) {
				add( current.substring(0,idx).trim(), current.substring(idx+keyValueDelimiter.length()).trim() );
			}
		}
	}
	
	/**
	 * This method actually adds a {@link SequenceAnnotation} to the internal list.
	 * 
	 * @param type the type of the {@link SequenceAnnotation}
	 * @param identifier the identifier of the {@link SequenceAnnotation}
	 */
	protected void add( String type, String identifier ) {
		annot.add( new SequenceAnnotation( type, identifier ) );
	}

	public void clearAnnotation() {
		annot.clear();
	}

	public SequenceAnnotation[] getCurrentAnnotation() {
		return annot.size() == 0 ? null : annot.toArray( empty );
	}

	public String parseAnnotationToComment( char commentChar, SequenceAnnotation... annotations ) {
		StringBuffer res = new StringBuffer();
		res.append(  commentChar );
		if( annotations != null && annotations.length > 0 ) {
			res.append( annotations[0].getType() + keyValueDelimiter + annotations[0].getIdentifier() );
			for( int i = 1; i < annotations.length; i++ ) {
				res.append( annotationDelimiter + " " + annotations[i].getType() + keyValueDelimiter + annotations[i].getIdentifier() );
			}
		}
		return res.toString();
	}
}
