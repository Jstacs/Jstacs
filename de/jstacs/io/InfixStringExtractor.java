/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.io;

import de.jstacs.data.sequences.annotation.SequenceAnnotation;

/**
 * This class implements an {@link AbstractStringExtractor} that can be seen as a filter.
 * Each instance only returns an infix instead of the full {@link String}. 
 * 
 * @author Jens Keilwagen
 */
public class InfixStringExtractor extends AbstractStringExtractor {

	private AbstractStringExtractor se;
	private int start, end;
	
	/**
	 * This constructor creates an instance that uses only a infix of the string returned by <code>se</code>.

	 * @param se the main {@link AbstractStringExtractor} returned the full {@link String}s
	 * @param start the start position of the infix
	 * @param length the length of the infix 
	 */
	public InfixStringExtractor( AbstractStringExtractor se, int start, int length ) {
		super( se.ignore );
		this.se = se;
		this.start = start;
		this.end = start+length;
		annotation = "infix data set (start=" + start +", length=" + length + ") of " + se.getAnnotation();
	}

	public boolean hasMoreElements() {
		return se.hasMoreElements();
	}

	public String nextElement() {
		return se.nextElement().substring( start, end );
	}
	
	public SequenceAnnotation[] getCurrentSequenceAnnotations() {
		return se.getCurrentSequenceAnnotations();
	}
}
