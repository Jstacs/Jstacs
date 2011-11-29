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

import de.jstacs.results.CategoricalResult;

/**
 * This class implements a naive {@link SequenceAnnotationParser} which simply paste the comments into {@link SequenceAnnotation}.
 * 
 * @author Jens Keilwagen
 */
public class SimpleSequenceAnnotationParser implements SequenceAnnotationParser {

	private LinkedList<SequenceAnnotation> annot;
	private static final SequenceAnnotation[] empty = new SequenceAnnotation[0];
	private static final String TYPE = "unparsed comment line";
	
	/**
	 * The constructor of a {@link SimpleSequenceAnnotationParser} which simply paste the comments into {@link SequenceAnnotation}.
	 */
	public SimpleSequenceAnnotationParser() {
		annot = new LinkedList<SequenceAnnotation>();
	}
	
	public void addToAnnotation( String unparsed ) {
		annot.add( new SequenceAnnotation( TYPE ,annot.size()+"", new CategoricalResult( "unparsed comment", "", unparsed.substring(1) ) ) );
	}

	public void clearAnnotation() {
		annot.clear();
	}

	public SequenceAnnotation[] getCurrentAnnotation() {
		return annot.size() == 0 ? null : annot.toArray( empty );
	}

	public String parseAnnotationToComment( char commentChar, SequenceAnnotation... annotations ) {
		String simple = "" + commentChar;
		if( annotations == null || annotations.length == 0 ) {
			return simple;
		} else {
			StringBuffer res = new StringBuffer();
			String start = "";
			for( int i = 0; i < annotations.length; i++ ) {
				if( annotations[i].getType().equals( TYPE )
						&& annotations[i].getNumberOfResults() == 1
						&& annotations[i].getResultAt( 0 ) instanceof CategoricalResult ) {
					res.append( start + commentChar +  annotations[i].getResultAt( 0 ).getResult().toString() );
					start = "\n";
				}
			}
			if( res.length() == 0 ) {
				return simple;
			} else {
				return res.toString();
			}
		}
	}
}
