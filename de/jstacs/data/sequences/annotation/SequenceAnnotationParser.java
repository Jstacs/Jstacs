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

/**
 * This interface declares the methods which are used by
 * {@link de.jstacs.io.AbstractStringExtractor} to annotate a {@link String}
 * which will be parsed to a {@link de.jstacs.data.Sequence}.
 * 
 * @author Jens Keilwagen
 */
public interface SequenceAnnotationParser {

	/**
	 * This method returns the current {@link SequenceAnnotation}.
	 * 
	 * @return the current {@link SequenceAnnotation}
	 */
	SequenceAnnotation[] getCurrentAnnotation();

	/**
	 * This method adds the <code>unparsed</code> {@link String} in some way to
	 * the {@link SequenceAnnotation}.
	 * 
	 * @param unparsed
	 *            the String containing the annotation
	 */
	void addToAnnotation( String unparsed );

	/**
	 * This method reset the current {@link SequenceAnnotation}.
	 */
	void clearAnnotation();

	/**
	 * This method returns a {@link String} representation of the given
	 * {@link SequenceAnnotation}s that can be used as comment line in a file.
	 * <br>
	 * Normally, the method should do the opposite of {@link SequenceAnnotationParser#addToAnnotation(String)}.
	 * 
	 * @param commentChar
	 *            the char at the beginning of comment lines
	 * @param annotations
	 *            the {@link SequenceAnnotation}s to be parsed
	 * 
	 * @return a String representing the given {@link SequenceAnnotation}s
	 * 
	 * @see de.jstacs.data.DataSet#save(java.io.OutputStream, char, SequenceAnnotationParser)
	 */
	String parseAnnotationToComment( char commentChar, SequenceAnnotation... annotations );
}
