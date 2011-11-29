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

package de.jstacs.utils;

/**
 * This class is a container for any objects that have to be compared. It is
 * enables you e.g. to use the method <code>Arrays.sort(...)</code> without
 * writing a comparator. Furthermore it helps to speed up the sorting since the
 * &quot;weight&quot; can be precomputed and does not have to be recomputed in
 * every comparison.
 * 
 * @author Jens Keilwagen
 * 
 * @param <E>
 *            any {@link Class}
 * @param <C>
 * 			  an implementation of {@link Comparable} that is used to compare instance of {@link ComparableElement}
 */
public class ComparableElement<E, C extends Comparable> implements Comparable<ComparableElement<E, C>> {

	/**
	 * The element.
	 */
	private E o;

	/**
	 * The weight.
	 */
	private C w;

	/**
	 * Creates a new {@link ComparableElement}.
	 * 
	 * @param o
	 *            the element
	 * @param w
	 *            the weight
	 */
	public ComparableElement( E o, C w ) {
		this.o = o;
		this.w = w;
	}

	/**
	 * This method returns the element.
	 * 
	 * @return the element
	 */
	public E getElement() {
		return o;
	}

	/**
	 * This method returns the weight of the element.
	 * 
	 * @return the weight of the element
	 */
	public C getWeight() {
		return w;
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo( ComparableElement<E, C> o ) throws ClassCastException {
		return w.compareTo( o.w );
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return o + "\nweight: " + w;
	}
}
