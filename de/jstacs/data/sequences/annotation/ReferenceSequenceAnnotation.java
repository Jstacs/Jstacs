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
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;
import de.jstacs.io.XMLParser;
import de.jstacs.results.Result;

/**
 * This class implements a {@link SequenceAnnotation} that contains a reference
 * sequence. This sequence can be obtained using the method
 * {@link #getReferenceSequence()}.
 * 
 * @author Jan Grau
 */
public class ReferenceSequenceAnnotation extends SequenceAnnotation {

	private Sequence ref;
	
	/**
	 * The type of a {@link ReferenceSequenceAnnotation}.
	 * 
	 * @see SequenceAnnotation#getType()
	 */
	public static final String TYPE = "reference"; 

	/**
	 * Creates a new {@link ReferenceSequenceAnnotation} with identifier
	 * <code>identifier</code>, reference sequence <code>ref</code>, and
	 * additional annotation (that does not fit the {@link SequenceAnnotation}
	 * definitions) given as a {@link Result} <code>result</code>.
	 * 
	 * @param identifier
	 *            the identifier of the annotation
	 * @param ref
	 *            the reference sequence
	 * @param results
	 *            the additional annotation
	 * 
	 * @see de.jstacs.results.ResultSet#ResultSet(Result)
	 */
	public ReferenceSequenceAnnotation( String identifier, Sequence ref, Result... results ) {
		super( TYPE, identifier, results );
		this.ref = ref;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ReferenceSequenceAnnotation} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ReferenceSequenceAnnotation} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer}
	 *             <code>representation</code> could not be parsed)
	 * 
	 * @see SequenceAnnotation#SequenceAnnotation(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public ReferenceSequenceAnnotation( StringBuffer representation ) throws NonParsableException {
		super( XMLParser.extractForTag( representation, "super" ) );
		representation = XMLParser.extractForTag( representation, getClass().getSimpleName() );
		AlphabetContainer cont = XMLParser.extractObjectForTags( representation, "alphabet", AlphabetContainer.class );
		try {
			ref = Sequence.create( cont, XMLParser.extractObjectForTags( representation, "ref", String.class ) );
		} catch ( Exception e ) {
			throw new NonParsableException();
		}
	}

	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags( buf, "super" );
		XMLParser.appendObjectWithTags( buf, ref.getAlphabetContainer(), "alphabet" );
		XMLParser.appendObjectWithTags( buf, ref.toString(), "ref" );
		XMLParser.addTags( buf, getClass().getSimpleName() );
		return buf;
	}

	public String toString() {
		return super.toString() + ref;
	}

	/**
	 * Returns the reference sequence.
	 * 
	 * @return the reference sequence
	 */
	public Sequence getReferenceSequence() {
		return ref;
	}
}