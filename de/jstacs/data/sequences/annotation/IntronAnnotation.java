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
 * Annotation class for an intron as defined by a donor and an acceptor splice
 * site.
 * 
 * @author Jan Grau
 */
public class IntronAnnotation extends LocatedSequenceAnnotationWithLength {

	/**
	 * Creates a new {@link IntronAnnotation} from a donor
	 * {@link SinglePositionSequenceAnnotation} and an acceptor
	 * {@link SinglePositionSequenceAnnotation} and a set of additional
	 * annotations.
	 * 
	 * @param identifier
	 *            the identifier of this annotation
	 * @param donor
	 *            the donor annotation
	 * @param acceptor
	 *            the acceptor annotation
	 * @param additionalAnnotation
	 *            the additional annotations
	 * 
	 * @throws Exception
	 *             if <code>donor</code> or <code>acceptor</code> are not of the
	 *             correct
	 *             {@link de.jstacs.data.sequences.annotation.SinglePositionSequenceAnnotation.Type}
	 * 
	 * @see SinglePositionSequenceAnnotation
	 * @see SinglePositionSequenceAnnotation.Type
	 */
	public IntronAnnotation( String identifier, SinglePositionSequenceAnnotation donor, SinglePositionSequenceAnnotation acceptor,
								Result... additionalAnnotation ) throws Exception {
		super( "Intron", identifier, new SinglePositionSequenceAnnotation[]{ donor, acceptor }, additionalAnnotation );
		if( !donor.getType().equals( SinglePositionSequenceAnnotation.Type.SDS.toString() ) || !acceptor.getType()
					.equals( SinglePositionSequenceAnnotation.Type.SAS.toString() ) ) {
			throw new Exception( "An intron is defined by donor and acceptor site." );
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link IntronAnnotation} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link IntronAnnotation} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer}
	 *             <code>representation</code> could not be parsed)
	 * 
	 * @see LocatedSequenceAnnotationWithLength#LocatedSequenceAnnotationWithLength(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public IntronAnnotation( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

}
