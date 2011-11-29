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

import java.util.Arrays;

/**
 * A simple list of primitive type <code>int</code>.
 * 
 * @author Jens Keilwagen
 */
public final class IntList implements Cloneable {

	private int size;

	private int[] array;

	/**
	 * This is the default constructor that creates an {@link IntList} with
	 * initial length 10.
	 */
	public IntList() {
		this( 10 );
	}

	/**
	 * This is the default constructor that creates an {@link IntList} with
	 * initial length <code>size</code>.
	 * 
	 * @param size
	 *            the initial size, has to be positive
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>size</code> is less than 1
	 */
	public IntList( int size ) throws IllegalArgumentException {
		super();
		if( size <= 0 ) {
			throw new IllegalArgumentException( "The size has to be positive." );
		}
		this.size = 0;
		array = new int[size];
	}
	
	public IntList clone() throws CloneNotSupportedException {
		IntList clone = (IntList) super.clone();
		clone.array = array.clone();
		return clone;
	}

	/**
	 * Adds the element <code>val</code> at the end of the list.
	 * 
	 * @param val
	 *            the element that should be inserted
	 */
	public void add( int val ) {
		if( array.length == size ) {
			int[] help = new int[2 * array.length];
			System.arraycopy( array, 0, help, 0, size );
			array = help;
		}
		array[size++] = val;
	}
	
	/**
	 * Adds <code>val</code> to the list only if
	 * it is not already contained in the list.
	 * @param val the value to be added
	 * @return the index of <code>val</code> if it is already contained in the list and -1 otherwise
	 */
	public int addConditional( int val ) {
		int b = contains( val );
		if(b<0){
			add(val);
		}
		return b;
	}

	/**
	 * Removes all elements from the list.
	 */
	public void clear() {
		size = 0;
	}

	/**
	 * Returns the element with the specified <code>index</code>.
	 * 
	 * @param index
	 *            the specified index of the element to return
	 * 
	 * @return the corresponding element
	 */
	public int get( int index ) {
		if( index < size ) {
			return array[index];
		} else {
			throw new ArrayIndexOutOfBoundsException( index );
		}
	}
	
	/**
	 * Checks if <code>val</code> is already returned in the list.
	 * @param val the value
	 * @return the index of <code>val</code> if it is contained in the list and -1 otherwise
	 */
	public int contains( int val ){
		for(int i=0;i<size;i++){
			if(array[i] == val){
				return i;
			}
		}
		return -1;
	}

	/**
	 * Returns the number of inserted elements.
	 * 
	 * @return the number of inserted elements
	 */
	public int length() {
		return size;
	}

	/**
	 * This method returns an <code>int</code> array containing all elements of
	 * the list.
	 * 
	 * @return an <code>int</code> array containing all elements of the list
	 */
	public int[] toArray() {
		int[] res = new int[size];
		System.arraycopy( array, 0, res, 0, size );
		return res;
	}
	
	/**
	 * This method reverses the list, i.e., after invoking this method the formerly
	 * first element is than the last, the formerly second element is than the last but one, ... 
	 */
	public void reverse() {
		int i = 0, n = size/2, help;
		while( i < n ) {
			help = array[i];
			array[i] = array[size-1-i];
			array[size-1-i++] = help;
		}
	}
	
	/**
	 * This method sorts the elements of the list.
	 * 
	 * @see Arrays#sort(int[], int, int)
	 */
	public void sort() {
		Arrays.sort( array, 0, size );
	}

	public boolean equals(Object o) {
		if( this == o ) {
		    return true;
		}
		if( o instanceof IntList ) {
			IntList list = (IntList)o;
		    if( size == list.size ) {
		    	int n = 0;
		    	while( n < size && array[n] != list.array[n] ) {
		    		n++;
		    	}
				return n == size;
			}
			return false;
		}
		return false;
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		for( int i = 0; i < size; i++) {
			sb.append( i==0 ? "[" : ", " );
			sb.append( array[i] );
		}
		sb.append( "]" );
		return sb.toString();
	}
}
