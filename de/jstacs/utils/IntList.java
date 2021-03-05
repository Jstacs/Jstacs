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
	 * Adds a constant to all internal values between start and end
	 * @param start the first index
	 * @param end the last index (exclusive)
	 * @param offset the constant to be added
	 */
	public void addToValues(int start, int end, int offset){
		for(int i=start;i<end;i++){
			array[i] += offset;
		}
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
	 * Returns the last element and removes it from the list.
	 * 
	 * @return the last element
	 */
	public int pop(){
		int temp = array[size-1];
		if(size > 0){
			size--;
		}
		return temp;
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

	
	/**
	 * Sorts the elements of the list and removes duplicate entries.
	 */
	public void sortAndMakeUnique() {
		sort();
		int i=0,j=0;
		while(i < size) {
			while(i < size-1 && array[i] == array[i+1]) {
				i++;
			}
			array[j] = array[i];
			j++;
			i++;
		}
		size = j;
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
		sb.append( "[" );
		for( int i = 0; i < size; i++) {
			sb.append( i==0 ? "" : ", " );
			sb.append( array[i] );
		}
		sb.append( "]" );
		return sb.toString();
	}

	/**
	 * Performs a binary search for element <code>key</code> by calling
	 * {@link Arrays#binarySearch(int[], int, int, int)} on the internal array.
	 * @param key the key to search for
	 * @param fromIndex the first element considered
	 * @param toIndex the first element not considered
	 * 
	 * @return see {@link Arrays#binarySearch(int[], int, int, int)}
	 */
	public int binarySearch(int key, int fromIndex, int toIndex ) {
		if(fromIndex<0 || toIndex > size){
			throw new ArrayIndexOutOfBoundsException();
		}
		return Arrays.binarySearch( array, fromIndex, toIndex, key );
	}
	
	/**
	 * Performs an interpolation search for element <code>key</code> on the internal array.
	 * @param key the key to search for
	 * @param fromIndex the first element considered
	 * @param toIndex the first element not considered
	 * 
	 * @return see {@link Arrays#binarySearch(int[], int, int, int)}
	 */
	public int interpolationSearch(int key, int fromIndex, int toIndex) {
		if(fromIndex<0 || toIndex > size){
			throw new ArrayIndexOutOfBoundsException();
		}
		return interpolationSearch(key, fromIndex, toIndex, array[fromIndex], array[toIndex-1]);
	}
	
	private int interpolationSearch(int idx2, int fromIndex, int toIndex, int first, int last) {
		int idx = (int)((idx2-first)/(double)(last-first)*(toIndex-fromIndex-1)) + fromIndex;
		int mid = array[idx];

		while(mid != idx2){
			if(mid>idx2){
				toIndex = idx;
				last = mid;
			}else{
				fromIndex = idx+1;
				first = array[idx+1];
			}
			
			idx = (int)((idx2-first)/(double)(last-first)*(toIndex-fromIndex-1)) + fromIndex;
			mid = array[idx];
		}
		return idx;
	}
	
}
