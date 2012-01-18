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

import de.jstacs.io.NonParsableException;
import de.jstacs.results.Result;

/**
 * Class for some annotations that consist mainly of one position on a sequence.
 * This includes transcription start sites, translation start sites, splice
 * donor sites and splice acceptor sites.
 * 
 * @author Jan Grau
 */
public class SinglePositionSequenceAnnotation extends LocatedSequenceAnnotation {

	/**
	 * This <code>enum</code> defines possible types of a
	 * {@link SinglePositionSequenceAnnotation}. May be extended in the future.
	 */
	public enum Type {

		/**
		 * Transcription start site.
		 */
		TSS( "Transcription start" ),
		/**
		 * Translation start site.
		 */
		TLS( "Translation start" ),
		/**
		 * Splice donor site.
		 */
		SDS( "Splice donor site" ),
		/**
		 * Splice acceptor site.
		 */
		SAS( "Splice acceptor site" );

		private String identifier;

		Type( String identifier ) {
			this.identifier = identifier;
		}

		/* (non-Javadoc)
		 * @see java.lang.Enum#toString()
		 */
		@Override
		public String toString() {
			return identifier;
		}

	}

	/**
	 * Creates a new {@link SinglePositionSequenceAnnotation} of type
	 * <code>type</code> with identifier <code>identifier</code> and position
	 * <code>position</code>.
	 * 
	 * @param position
	 *            the position of the {@link SinglePositionSequenceAnnotation}
	 *            on the sequence
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * 
	 * @see LocatedSequenceAnnotation#LocatedSequenceAnnotation(int, String,
	 *      String, Result...)
	 */
	public SinglePositionSequenceAnnotation( Type type, String identifier, int position ) {
		super( position, type.toString(), identifier );
	}

	/**
	 * Creates a new {@link SinglePositionSequenceAnnotation} of type
	 * <code>type</code> with identifier <code>identifier</code>, position
	 * <code>position</code> and additional annotations
	 * <code>additionalAnnotation</code>.
	 * 
	 * @param position
	 *            the position of the {@link SinglePositionSequenceAnnotation}
	 *            on the sequence
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param additionalAnnotation
	 *            the additional annotation
	 * 
	 * @see LocatedSequenceAnnotation#LocatedSequenceAnnotation(int, String,
	 *      String, Result...)
	 */
	public SinglePositionSequenceAnnotation( Type type, String identifier, int position, Result... additionalAnnotation ) {
		super( position, type.toString(), identifier, additionalAnnotation );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SinglePositionSequenceAnnotation} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SinglePositionSequenceAnnotation} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see LocatedSequenceAnnotation#LocatedSequenceAnnotation(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public SinglePositionSequenceAnnotation( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

}
