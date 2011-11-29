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

package de.jstacs.data;

import java.util.Collection;

/**
 * This class implements a {@link RecyclableSequenceEnumerator} on user-specified {@link Sequence}s.
 * 
 * @author Jens Keilwagen
 */
public class SequenceEnumeration implements RecyclableSequenceEnumerator {

	private Sequence[] content;
	private int index;

	/**
	 * This constructor creates an instance based on the user-specified {@link Sequence}s <code>sequences</code>.
	 * 
	 * @param sequences the user-specified {@link Sequence}s
	 */
	public SequenceEnumeration( Sequence... sequences ) {
		content = sequences.clone();
		reset();
	}
	
	/**
	 * This constructor creates an instance based on the user-specified {@link Collection} of {@link Sequence}s <code>sequences</code>.
	 * 
	 * @param sequences the user-specified {@link Collection} of {@link Sequence}s
	 */
	public SequenceEnumeration( Collection<Sequence> sequences ) {
		content = sequences.toArray( new Sequence[0] );
		reset();
	}
	
	public void reset() {
		index = 0;
	}

	public boolean hasMoreElements() {
		return index < content.length;
	}

	public Sequence nextElement() {
		return content[index++];
	}
}
