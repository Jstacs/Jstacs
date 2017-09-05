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
 * A very simple container.
 * 
 * @author Jens Keilwagen
 *
 * @param <E1> the type of the first element
 * @param <E2> the type of the second element
 */
public class Pair<E1,E2> {
	private E1 element1;
	private E2 element2;
	
	/**
	 * This constructor creates a pair of <code>element1</code> and <code>element2</code>.
	 * 
	 * @param element1 the first element
	 * @param element2 the second element
	 */
	public Pair( E1 element1, E2 element2 ) {
		this.element1 = element1;
		this.element2 = element2;
	}
	
	/**
	 * This method returns the first element.
	 * @return the first element
	 */
	public E1 getFirstElement() {
		return element1;
	}
	
	/**
	 * This method returns the second element.
	 * @return the second element
	 */
	public E2 getSecondElement() {
		return element2;
	}
	
	
	public String toString(){
		return "["+element1+",\n"+element2+"]";
	}
}
