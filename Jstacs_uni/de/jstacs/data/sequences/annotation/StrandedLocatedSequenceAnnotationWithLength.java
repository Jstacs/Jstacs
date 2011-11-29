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
 * Class for a {@link SequenceAnnotation} that has a position, a length and an
 * orientation on the strand of a {@link de.jstacs.data.Sequence}.
 * 
 * @author Jan Grau
 */
public class StrandedLocatedSequenceAnnotationWithLength extends LocatedSequenceAnnotationWithLength {

	/**
	 * This enum defines possible orientations on the strands.
	 * 
	 * @author Jan Grau
	 */
	public enum Strand {
		/**
		 * The annotation is located on the forward strand.
		 */
		FORWARD( "Forward strand" ),
		/**
		 * The annotation is located on the reverse strand.
		 */
		REVERSE( "Reverse strand" ),
		/**
		 * The orientation of the annotation is not known.
		 */
		UNKNOWN( "Unknown strand" ), ;

		private final String strandedness;

		Strand( String strandedness ) {
			this.strandedness = strandedness;
		}

		/**
		 * Returns the strandedness, i.e. the orientation on the strand of the
		 * sequence as a {@link String}.
		 * 
		 * @return the strandedness, i.e. the orientation on the strand of the
		 *         sequence
		 */
		public String strandedness() {
			return strandedness;
		}
	};

	private Strand strandedness;

	/**
	 * Creates a new {@link StrandedLocatedSequenceAnnotationWithLength} of type
	 * <code>type</code> with identifier <code>identifier</code> and additional
	 * annotation (that does not fit the {@link SequenceAnnotation} definitions)
	 * given as an array of {@link Result}s <code>results</code>.
	 * 
	 * @param position
	 *            the position of the
	 *            {@link StrandedLocatedSequenceAnnotationWithLength} on the
	 *            sequence
	 * @param length
	 *            the length of the
	 *            {@link StrandedLocatedSequenceAnnotationWithLength}
	 * @param strandedness
	 *            the orientation on the strand
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param results
	 *            the additional annotation
	 * 
	 * @see LocatedSequenceAnnotationWithLength#LocatedSequenceAnnotationWithLength(int,
	 *      int, String, String, Result...)
	 */
	public StrandedLocatedSequenceAnnotationWithLength( int position, int length, Strand strandedness, String type, String identifier,
														Result... results ) {
		super( position, length, type, identifier, results );
		this.strandedness = strandedness;
	}

	/**
	 * Creates a new {@link StrandedLocatedSequenceAnnotationWithLength} of type
	 * <code>type</code> with identifier <code>identifier</code> and additional
	 * annotation (that does not fit the {@link SequenceAnnotation} definitions)
	 * given as a {@link Collection} of {@link Result}s <code>results</code>.
	 * 
	 * @param position
	 *            the position of the
	 *            {@link StrandedLocatedSequenceAnnotationWithLength} on the
	 *            sequence
	 * @param length
	 *            the length of the
	 *            {@link StrandedLocatedSequenceAnnotationWithLength}
	 * @param strandedness
	 *            the orientation on the strand
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param results
	 *            the additional annotation
	 * 
	 * @see LocatedSequenceAnnotationWithLength#LocatedSequenceAnnotationWithLength(int,
	 *      int, String, String, Collection)
	 */
	public StrandedLocatedSequenceAnnotationWithLength( int position, int length, Strand strandedness, String type, String identifier,
														Collection<Result> results ) {
		super( position, length, type, identifier, results );
		this.strandedness = strandedness;
	}

	/**
	 * Creates a new {@link StrandedLocatedSequenceAnnotationWithLength} of type
	 * <code>type</code> with identifier <code>identifier</code>, additional
	 * annotation (that does not fit the {@link SequenceAnnotation} definitions)
	 * given as an array of {@link Result}s <code>additionalAnnotations</code>
	 * and sub-annotations <code>annotations</code>.
	 * 
	 * @param position
	 *            the position of the
	 *            {@link StrandedLocatedSequenceAnnotationWithLength} on the
	 *            sequence
	 * @param length
	 *            the length of the
	 *            {@link StrandedLocatedSequenceAnnotationWithLength}
	 * @param strandedness
	 *            the orientation on the strand
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param annotations
	 *            the sub-annotations
	 * @param additionalAnnotations
	 *            the additional annotation
	 * 
	 * @see LocatedSequenceAnnotationWithLength#LocatedSequenceAnnotationWithLength(int,
	 *      int, String, String, SequenceAnnotation[], Result...)
	 */
	public StrandedLocatedSequenceAnnotationWithLength( int position, int length, Strand strandedness, String type, String identifier,
														SequenceAnnotation[] annotations, Result... additionalAnnotations ) {
		super( position, length, type, identifier, annotations, additionalAnnotations );
		this.strandedness = strandedness;
	}

	/**
	 * Creates a new {@link StrandedLocatedSequenceAnnotationWithLength} of type
	 * <code>type</code> with identifier <code>identifier</code>, additional
	 * annotation (that does not fit the {@link SequenceAnnotation} definitions)
	 * given as an array of {@link Result}s <code>additionalAnnotations</code>
	 * and sub-annotations <code>annotations</code>. The position of the new
	 * {@link LocatedSequenceAnnotationWithLength} is the minimal position of
	 * all positions of <code>annotations</code> and the length is determined so
	 * that it is the maximum of these positions and (if applicable) the
	 * corresponding values of {@link #getEnd()}.
	 * 
	 * @param strandedness
	 *            the orientation on the strand
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param annotations
	 *            the sub-annotations
	 * @param additionalAnnotations
	 *            the additional annotation
	 * 
	 * @see LocatedSequenceAnnotationWithLength#LocatedSequenceAnnotationWithLength(String,
	 *      String, LocatedSequenceAnnotation[], Result...)
	 * @see StrandedLocatedSequenceAnnotationWithLength#getEnd()
	 */
	public StrandedLocatedSequenceAnnotationWithLength( String type, String identifier, Strand strandedness,
														LocatedSequenceAnnotation[] annotations, Result... additionalAnnotations ) {
		super( type, identifier, annotations, additionalAnnotations );
		this.strandedness = strandedness;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link StrandedLocatedSequenceAnnotationWithLength} out of
	 * its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StrandedLocatedSequenceAnnotationWithLength}
	 *             could not be reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see LocatedSequenceAnnotationWithLength#LocatedSequenceAnnotationWithLength(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public StrandedLocatedSequenceAnnotationWithLength( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * Returns the strandedness, i.e the orientation of this annotation.
	 * 
	 * @return the strandedness, i.e. the orientation
	 */
	public Strand getStrandedness() {
		return strandedness;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.sequences.annotation.LocatedSequenceAnnotationWithLength#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML( StringBuffer representation ) throws NonParsableException {
		representation = XMLParser.extractForTag( representation, "strandedAnnotation" );
		super.fromXML( XMLParser.extractForTag( representation, "locatedSequenceAnnotationWithLength" ) );
		strandedness = XMLParser.extractObjectForTags( representation, "strandedness", Strand.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.sequences.annotation.LocatedSequenceAnnotationWithLength#toString()
	 */
	@Override
	public String toString() {
		StringBuffer buf = new StringBuffer( super.toString() );
		buf.append( "strand: " );
		buf.append( strandedness.strandedness() );
		buf.append( "\n" );
		return buf.toString();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.sequences.annotation.LocatedSequenceAnnotationWithLength#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags( buf, "locatedSequenceAnnotationWithLength" );
		XMLParser.appendObjectWithTags( buf, strandedness, "strandedness" );
		XMLParser.addTags( buf, "strandedAnnotation" );
		return buf;
	}

}
