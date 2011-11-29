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

package de.jstacs.io;

import java.util.Enumeration;
import java.util.regex.Pattern;

import de.jstacs.data.sequences.annotation.SequenceAnnotation;

/**
 * This class implements the reader that extracts strings. The class ignores
 * lines starting with a given character, since those lines are treated as
 * comments. If the user does not specify this character, it is set to
 * &quot;#&quot; internally. If the user specifies this character as
 * &quot;&gt;&quot;, the file or {@link String} will be treated as in
 * FastA-format, i.e. lines beginning with &quot;&gt;&quot; will be stripped and
 * the lines between two &quot;&gt;&quot; (or until the end of the file) will be
 * appended to form a new {@link String}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public abstract class AbstractStringExtractor implements Enumeration<String> {

	/**
	 * The usual comment character is &quot;#&quot;. Lines beginning with this
	 * sign will be ignored.
	 */
	public final static char USUALLY = '#';

	/**
	 * The comment character for FastA-formatted files is &quot;>&quot;. If
	 * &quot;>&quot; is specified as the comment character, the file or
	 * {@link String} will be interpreted as in FastA format.
	 */
	public final static char FASTA = '>';

	/**
	 * The annotation of the source.
	 */
	protected String annotation;
	
	/**
	 * The pattern for ignoring comment lines.
	 */
	protected Pattern ignorePattern;

	/**
	 * The internal comment character.
	 */
	protected char ignore;

	/**
	 * Creates a new {@link AbstractStringExtractor} with the specified
	 * character as start of each comment line.
	 * 
	 * @param ignore
	 *            the comment character
	 * 
	 * @see Pattern
	 */
	protected AbstractStringExtractor( char ignore ) {
		this.ignore = ignore;
		this.ignorePattern = Pattern.compile( "^\\s*" + ignore + ".*" );
	}

	/**
	 * Returns the annotation of the source.
	 * 
	 * @return the annotation
	 */
	public final String getAnnotation() {
		return annotation;
	}
	
	/**
	 * Returns the {@link SequenceAnnotation} or <code>null</code> if no {@link SequenceAnnotation} is available.
	 * 
	 * @return the {@link SequenceAnnotation} or <code>null</code> if no {@link SequenceAnnotation} is available.
	 */
	public SequenceAnnotation[] getCurrentSequenceAnnotations() {
		return null;
	}
}
