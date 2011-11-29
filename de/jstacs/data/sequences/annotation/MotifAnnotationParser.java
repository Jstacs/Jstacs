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

import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;

/**
 * Annotation parser for {@link MotifAnnotation}s.
 * 
 * @author Jan Grau
 *
 */
public class MotifAnnotationParser extends SplitSequenceAnnotationParser {

	/**
	 * Creates a new {@link MotifAnnotationParser} with default delimiters.
	 * @see SplitSequenceAnnotationParser#SplitSequenceAnnotationParser()
	 */
	public MotifAnnotationParser() {
		super();
	}

	/**
	 * Creates a new {@link MotifAnnotationParser} with the supplied delimiters
	 * @param keyValueDelimiter the delimiter between key and value
	 * @param annotationDelimiter the delimiter between different annotations
	 */
	public MotifAnnotationParser( String keyValueDelimiter, String annotationDelimiter ) {
		super( keyValueDelimiter, annotationDelimiter );
	}

	@Override
	public void addToAnnotation( String unparsed ) {
		String[] split = unparsed.substring(1).split( annotationDelimiter );
		for( String current : split ) {
			int idx = current.indexOf( keyValueDelimiter );
			if( idx >= 0 ) {
				String type = current.substring(0,idx).trim();
				String identifier = current.substring(idx+keyValueDelimiter.length()).trim();
				int idxl, idxr;
				if("Motif".equals( type ) &&  (idxl = identifier.indexOf( "(" ) ) >=0 && (idxr = identifier.indexOf( ")") ) >=0 ){
					
					String ann = identifier.substring( idxl+1, idxr ).trim();
					String[] parts = ann.split( ", " );
					int position = Integer.parseInt( parts[0] );
					int length = Integer.parseInt( parts[1] );
					Strand strandedness = Strand.valueOf( parts[2].trim() );
					
					add( identifier.substring( 0, idxl ).trim(), position, length, strandedness );
				}else{
					add( type, identifier );
				}
			}
		}
	}
	
	/**
	 * Adds the motif with identifier <code>identifier</code> at position <code>position</code>
	 * with length <code>length</code> and {@link Strand} <code>strandedness</code>
	 * @param identifier the identifier of the motif
	 * @param position the position of the motif
	 * @param length the length of the motif
	 * @param strandedness the strand of the motif
	 */
	protected void add( String identifier, int position, int length, Strand strandedness ) {
		annot.add( new MotifAnnotation( identifier, position, length, strandedness ) );
	}
	
	public String parseAnnotationToComment( char commentChar, SequenceAnnotation... annotations ) {
		StringBuffer res = new StringBuffer();
		res.append(  commentChar );
		MotifAnnotation ma = null;
		if( annotations != null && annotations.length > 0 ) {
			res.append( annotations[0].getType() + keyValueDelimiter + annotations[0].getIdentifier() );
			if( annotations[0] instanceof MotifAnnotation ){
				ma = (MotifAnnotation) annotations[0];
				res.append( "( "+ma.getPosition()+", "+ma.getLength()+", "+ma.getStrandedness().name()+" )" );
			}
			for( int i = 1; i < annotations.length; i++ ) {
				res.append( annotationDelimiter + " " + annotations[i].getType() + keyValueDelimiter + annotations[i].getIdentifier() );
				if( annotations[i] instanceof MotifAnnotation ){
					ma = (MotifAnnotation) annotations[i];
					res.append( " ( "+ma.getPosition()+", "+ma.getLength()+", "+ma.getStrandedness().name()+" )" );
				}
			}
		}
		return res.toString();
	}

}
