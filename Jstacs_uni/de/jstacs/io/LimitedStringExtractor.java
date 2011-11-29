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
 * This {@link AbstractStringExtractor} allows to read only a limited amount of data.
 * 
 * @author Jens Keilwagen
 */
public class LimitedStringExtractor extends AbstractStringExtractor {

	private AbstractStringExtractor se;
	private int num, current;
	
	/**
	 * This constructor creates a new instance based on a given {@link AbstractStringExtractor} and a maximal number of Strings to be read.
	 *  
	 * @param se the underlying {@link AbstractStringExtractor}
	 * @param num the maximal number of Strings to be read 
	 */
	public LimitedStringExtractor( AbstractStringExtractor se, int num ) {
		super( se.ignore );
		this.se = se;
		if( num < 0 ) {
			throw new IllegalArgumentException( "The number has to be non-negative." );
		}
		this.num = num;
		current = 0;
		annotation = "limited sample (max=" + num + ") of " + se.getAnnotation();
	}

	public boolean hasMoreElements() {
		return current < num && se.hasMoreElements();
	}

	public String nextElement() {
		String res = se.nextElement();
		current++;
		return res;
	}
	
	public SequenceAnnotation[] getCurrentSequenceAnnotations() {
		return se.getCurrentSequenceAnnotations();
	}
}
