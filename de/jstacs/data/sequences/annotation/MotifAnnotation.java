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

import de.jstacs.NonParsableException;
import de.jstacs.results.Result;

/**
 * Class for a {@link StrandedLocatedSequenceAnnotationWithLength} that is a
 * motif. The usefulness of this class amounts to having a defined type for all
 * motifs.
 * 
 * @author Jan Grau
 */
public class MotifAnnotation extends StrandedLocatedSequenceAnnotationWithLength {

	/**
	 * Creates a new {@link MotifAnnotation} of type <code>type</code> with
	 * identifier <code>identifier</code> and additional annotation (that does
	 * not fit the {@link SequenceAnnotation} definitions) given as an array of
	 * {@link Result}s <code>additionalAnnotation</code>. The orientation of the
	 * motif on the strand is also considered by the value
	 * <code>strandedness</code> of the <code>enum</code> {@link StrandedLocatedSequenceAnnotationWithLength.Strand}.
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
	 * @param identifier
	 *            the identifier of the annotation
	 * @param additionalAnnotation
	 *            the additional annotation
	 * 
	 * @see StrandedLocatedSequenceAnnotationWithLength.Strand
	 * @see StrandedLocatedSequenceAnnotationWithLength#StrandedLocatedSequenceAnnotationWithLength(int,
	 *      int, StrandedLocatedSequenceAnnotationWithLength.Strand, String, String, Result[])
	 */
	public MotifAnnotation( String identifier, int position, int length, Strand strandedness, Result... additionalAnnotation ) {
		super( position, length, strandedness, "Motif", identifier, additionalAnnotation );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link MotifAnnotation} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MotifAnnotation} could not be reconstructed out
	 *             of the XML representation (the {@link StringBuffer}
	 *             <code>representation</code> could not be parsed)
	 * 
	 * @see StrandedLocatedSequenceAnnotationWithLength#StrandedLocatedSequenceAnnotationWithLength(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public MotifAnnotation( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

}
