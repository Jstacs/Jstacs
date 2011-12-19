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

import java.util.Collection;

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.Result;

/**
 * Class for a {@link SequenceAnnotation} that has a position on the sequence
 * and a length, e.g. for donor splice sites, exons or genes.
 * 
 * @author Jan Grau
 */
public class LocatedSequenceAnnotationWithLength extends LocatedSequenceAnnotation {

	private int length;

	/**
	 * Creates a new {@link LocatedSequenceAnnotationWithLength} of type
	 * <code>type</code> with identifier <code>identifier</code> and additional
	 * annotation (that does not fit the {@link SequenceAnnotation} definitions)
	 * given as an array of {@link Result}s <code>result</code>.
	 * 
	 * @param position
	 *            the position of the
	 *            {@link LocatedSequenceAnnotationWithLength} on the sequence
	 * @param length
	 *            the length of the {@link LocatedSequenceAnnotationWithLength}
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param results
	 *            the additional annotation
	 * 
	 * @see LocatedSequenceAnnotation#LocatedSequenceAnnotation(int, String,
	 *      String, Result...)
	 */
	public LocatedSequenceAnnotationWithLength( int position, int length, String type, String identifier, Result... results ) {
		super( position, type, identifier, results );
		this.length = length;
	}

	/**
	 * Creates a new {@link LocatedSequenceAnnotationWithLength} of type
	 * <code>type</code> with identifier <code>identifier</code> and additional
	 * annotation (that does not fit the {@link SequenceAnnotation} definitions)
	 * given as a {@link Collection} of {@link Result}s <code>result</code>.
	 * 
	 * @param position
	 *            the position of the
	 *            {@link LocatedSequenceAnnotationWithLength} on the sequence
	 * @param length
	 *            the length of the {@link LocatedSequenceAnnotationWithLength}
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param results
	 *            the additional annotation
	 * 
	 * @see LocatedSequenceAnnotation#LocatedSequenceAnnotation(int, String,
	 *      String, Collection)
	 */
	public LocatedSequenceAnnotationWithLength( int position, int length, String type, String identifier, Collection<Result> results ) {
		super( position, type, identifier, results );
		this.length = length;
	}

	/**
	 * Creates a new {@link LocatedSequenceAnnotationWithLength} of type
	 * <code>type</code> with identifier <code>identifier</code>, additional
	 * annotation (that does not fit the {@link SequenceAnnotation} definitions)
	 * given as an array of {@link Result}s <code>additionalAnnotations</code>
	 * and sub-annotations <code>annotations</code>.
	 * 
	 * @param position
	 *            the position of the
	 *            {@link LocatedSequenceAnnotationWithLength} on the sequence
	 * @param length
	 *            the length of the {@link LocatedSequenceAnnotationWithLength}
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param annotations
	 *            the sub-annotations
	 * @param additionalAnnotations
	 *            the additional annotation
	 * 
	 * @see LocatedSequenceAnnotation#LocatedSequenceAnnotation(int, String,
	 *      String, SequenceAnnotation[], Result...)
	 */
	public LocatedSequenceAnnotationWithLength( int position, int length, String type, String identifier, SequenceAnnotation[] annotations,
												Result... additionalAnnotations ) {
		super( position, type, identifier, annotations, additionalAnnotations );
		this.length = length;
	}

	/**
	 * Creates a new {@link LocatedSequenceAnnotationWithLength} of type
	 * <code>type</code> with identifier <code>identifier</code>, additional
	 * annotation (that does not fit the {@link SequenceAnnotation} definitions)
	 * given as an array of {@link Result}s <code>additionalAnnotations</code>
	 * and sub-annotations <code>annotations</code>. The position of the new
	 * {@link LocatedSequenceAnnotationWithLength} is the minimal position of
	 * all positions of <code>annotations</code> and the length is determined so
	 * that it is the maximum of these positions and (if applicable) the
	 * corresponding values of {@link #getEnd()}.
	 * 
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param annotations
	 *            the sub-annotations
	 * @param additionalAnnotations
	 *            the additional annotation
	 * 
	 * @see LocatedSequenceAnnotation#LocatedSequenceAnnotation(String, String,
	 *      LocatedSequenceAnnotation[], Result...)
	 * @see LocatedSequenceAnnotationWithLength#getEnd()
	 */
	public LocatedSequenceAnnotationWithLength( String type, String identifier, LocatedSequenceAnnotation[] annotations,
												Result... additionalAnnotations ) {
		super( type, identifier, annotations, additionalAnnotations );
		int end = -Integer.MAX_VALUE;
		for( int i = 0; i < annotations.length; i++ ) {
			if( annotations[i] instanceof LocatedSequenceAnnotationWithLength ) {
				if( ( (LocatedSequenceAnnotationWithLength)annotations[i] ).getEnd() > end ) {
					end = ( (LocatedSequenceAnnotationWithLength)annotations[i] ).getEnd();
				}
			} else {
				if( annotations[i].getPosition() + 1 > end ) {
					end = annotations[i].getPosition() + 1;
				}
			}
		}
		this.length = end - getPosition();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link LocatedSequenceAnnotationWithLength} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link LocatedSequenceAnnotationWithLength} could not
	 *             be reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see LocatedSequenceAnnotation#LocatedSequenceAnnotation(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public LocatedSequenceAnnotationWithLength( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.sequences.annotation.LocatedSequenceAnnotation#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML( StringBuffer representation ) throws NonParsableException {
		representation = XMLParser.extractForTag( representation, "locatedSequenceAnnotationWithLength" );
		super.fromXML( XMLParser.extractForTag( representation, "locatedSequenceAnnotation" ) );
		length = XMLParser.extractObjectForTags( representation, "length", int.class );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.sequences.annotation.LocatedSequenceAnnotation#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags( buf, "locatedSequenceAnnotation" );
		XMLParser.appendObjectWithTags( buf, length, "length" );
		XMLParser.addTags( buf, "locatedSequenceAnnotationWithLength" );
		return buf;
	}

	/**
	 * Returns the length of this {@link LocatedSequenceAnnotationWithLength} as
	 * given in the constructor.
	 * 
	 * @return the length of this {@link LocatedSequenceAnnotationWithLength}
	 */
	public int getLength() {
		return length;
	}

	/**
	 * Returns the end of this {@link LocatedSequenceAnnotationWithLength}, i.e.
	 * {@link #getPosition()} + {@link #getLength()}.
	 * 
	 * @return the end of this {@link LocatedSequenceAnnotationWithLength}
	 * 
	 * @see LocatedSequenceAnnotationWithLength#getPosition()
	 * @see LocatedSequenceAnnotationWithLength#getLength()
	 */
	public int getEnd() {
		return getPosition() + length;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.sequences.annotation.LocatedSequenceAnnotation#toString()
	 */
	@Override
	public String toString() {
		StringBuffer buf = new StringBuffer();
		buf.append( super.toString() );
		buf.append( "length: " );
		buf.append( length );
		buf.append( "\n" );
		return buf.toString();
	}
	
	/**
	 * Returns <code>true</code> if this {@link LocatedSequenceAnnotationWithLength} overlaps
	 * with the location of <code>second</code>.
	 * @param second the other {@link LocatedSequenceAnnotationWithLength}
	 * @return if both annotations overlap
	 */
	public boolean overlaps(LocatedSequenceAnnotationWithLength second){
		return ( (this.getPosition() <= second.getPosition() && this.getEnd() > second.getPosition() ) ||
				(this.getPosition() >= second.getPosition() && second.getEnd() > this.getPosition()) );
	}

}
