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

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.Result;

/**
 * Class for a {@link SequenceAnnotation} that has a position on the sequence,
 * e.g for transcription start sites or intron-exon junctions.
 * 
 * @author Jan Grau
 */
public class LocatedSequenceAnnotation extends SequenceAnnotation {

	private int position;

	/**
	 * Creates a new {@link LocatedSequenceAnnotation} of type <code>type</code>
	 * with identifier <code>identifier</code> and additional annotation (that
	 * does not fit the {@link SequenceAnnotation} definitions) given as an
	 * array of {@link Result}s <code>result</code>.
	 * 
	 * @param position
	 *            the position of the {@link LocatedSequenceAnnotation} on the
	 *            sequence
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param results
	 *            the additional annotation
	 * 
	 * @see SequenceAnnotation#SequenceAnnotation(String, String, Result[][])
	 */
	public LocatedSequenceAnnotation( int position, String type, String identifier, Result... results ) {
		super( type, identifier, results );
		this.position = position;
	}

	/**
	 * Creates a new {@link LocatedSequenceAnnotation} of type <code>type</code>
	 * with identifier <code>identifier</code> and additional annotation (that
	 * does not fit the {@link SequenceAnnotation} definitions) given as a
	 * {@link Collection} of {@link Result}s <code>result</code>.
	 * 
	 * @param position
	 *            the position of the {@link LocatedSequenceAnnotation} on the
	 *            sequence
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param results
	 *            the additional annotation
	 * 
	 * @see SequenceAnnotation#SequenceAnnotation(String, String, Collection)
	 */
	public LocatedSequenceAnnotation( int position, String type, String identifier, Collection<Result> results ) {
		super( type, identifier, results );
		this.position = position;
	}

	/**
	 * Creates a new {@link LocatedSequenceAnnotation} of type <code>type</code>
	 * with identifier <code>identifier</code>, additional annotation (that does
	 * not fit the {@link SequenceAnnotation} definitions) given as an array of
	 * {@link Result}s <code>additionalAnnotation</code> and sub-annotations
	 * <code>annotations</code>.
	 * 
	 * @param position
	 *            the position of the {@link LocatedSequenceAnnotation} on the
	 *            sequence
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param annotations
	 *            the sub-annotations
	 * @param additionalAnnotation
	 *            the additional annotation
	 * 
	 * @see SequenceAnnotation#SequenceAnnotation(String, String,
	 *      SequenceAnnotation[], Result...)
	 */
	public LocatedSequenceAnnotation( int position, String type, String identifier, SequenceAnnotation[] annotations,
										Result... additionalAnnotation ) {
		super( type, identifier, annotations, additionalAnnotation );
		this.position = position;
	}

	/**
	 * Creates a new {@link LocatedSequenceAnnotation} of type <code>type</code>
	 * with identifier <code>identifier</code>, additional annotation (that does
	 * not fit the {@link SequenceAnnotation} definitions) given as an array of
	 * {@link Result}s <code>additionalAnnotation</code> and sub-annotations
	 * <code>annotations</code>. The position of this
	 * {@link LocatedSequenceAnnotation} is the minimum of the positions of
	 * <code>annotations</code>.
	 * 
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param annotations
	 *            the sub-annotations
	 * @param additionalAnnotation
	 *            the additional annotation
	 * 
	 * @see SequenceAnnotation#SequenceAnnotation(String, String,
	 *      SequenceAnnotation[], Result...)
	 */
	public LocatedSequenceAnnotation( String type, String identifier, LocatedSequenceAnnotation[] annotations,
										Result... additionalAnnotation ) {
		super( type, identifier, annotations, additionalAnnotation );
		this.position = Integer.MAX_VALUE;
		for( int i = 0; i < annotations.length; i++ ) {
			if( ( (LocatedSequenceAnnotation)annotations[i] ).getPosition() < position ) {
				position = ( (LocatedSequenceAnnotation)annotations[i] ).getPosition();
			}
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link LocatedSequenceAnnotation} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link LocatedSequenceAnnotation} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see SequenceAnnotation#SequenceAnnotation(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public LocatedSequenceAnnotation( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.sequences.annotation.SequenceAnnotation#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML( StringBuffer representation ) throws NonParsableException {
		representation = XMLParser.extractForTag( representation, "locatedSequenceAnnotation" );
		super.fromXML( XMLParser.extractForTag( representation, "annotation" ) );
		position = XMLParser.extractObjectForTags( representation, "position", int.class );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.sequences.annotation.SequenceAnnotation#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags( buf, "annotation" );
		XMLParser.appendObjectWithTags( buf, position, "position" );
		XMLParser.addTags( buf, "locatedSequenceAnnotation" );
		return buf;
	}

	/**
	 * Returns the position of this {@link LocatedSequenceAnnotation} on the
	 * sequence.
	 * 
	 * @return the position on the sequence
	 */
	public int getPosition() {
		return position;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.sequences.annotation.SequenceAnnotation#toString()
	 */
	@Override
	public String toString() {
		StringBuffer buf = new StringBuffer();
		buf.append( super.toString() );
		buf.append( "position: " );
		buf.append( position );
		buf.append( "\n" );
		return buf.toString();
	}

}
