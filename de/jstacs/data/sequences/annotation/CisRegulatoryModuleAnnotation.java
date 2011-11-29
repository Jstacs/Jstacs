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
 * Annotation for a cis-regulatory module as defined by a set of
 * {@link MotifAnnotation}s of the motifs in the module.
 * 
 * @author Jan Grau
 * 
 * @see MotifAnnotation
 */
public class CisRegulatoryModuleAnnotation extends SequenceAnnotation {

	/**
	 * Creates a new {@link CisRegulatoryModuleAnnotation} from a set of motifs
	 * and possibly additional annotations.
	 * 
	 * @param identifier
	 *            the identifier of this annotation
	 * @param motifs
	 *            the motifs in the module
	 * @param additionalAnnotation
	 *            the additional annotations
	 * 
	 * @see SequenceAnnotation#SequenceAnnotation(String, String,
	 *      SequenceAnnotation[], Result...)
	 */
	public CisRegulatoryModuleAnnotation( String identifier, MotifAnnotation[] motifs, Result... additionalAnnotation ) {
		super( "Cis-regulatory module", identifier, motifs, additionalAnnotation );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link CisRegulatoryModuleAnnotation} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link CisRegulatoryModuleAnnotation} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see SequenceAnnotation#SequenceAnnotation(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public CisRegulatoryModuleAnnotation( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

}
