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
import java.util.Map;
import java.util.TreeMap;

import de.jstacs.Storable;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * A simple list of primitive type <code>double</code>.
 * 
 * @author Jens Keilwagen
 */
public final class DoubleList implements Storable, Cloneable {

	private int size;

	private double[] array, help;

	/**
	 * This is the default constructor that creates a {@link DoubleList} with
	 * initial length 10.
	 */
	public DoubleList() {
		this( 10 );
	}

	/**
	 * This is the default constructor that creates a {@link DoubleList} with
	 * initial length <code>size</code>.
	 * 
	 * @param size
	 *            the initial size, has to be positive
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>size</code> is less than 1
	 */
	public DoubleList( int size ) throws IllegalArgumentException {
		super();
		if( size <= 0 ) {
			throw new IllegalArgumentException( "The size has to be positive." );
		}
		this.size = 0;
		array = new double[size];
	}
	
	public DoubleList clone() throws CloneNotSupportedException {
		DoubleList clone = (DoubleList) super.clone();
		clone.array = array.clone();
		clone.help = null;
		return clone;
	}

	/**
	 * This is the constructor for the interface {@link Storable}. Constructs a
	 * {@link DoubleList} out of an XML representation.
	 * 
	 * @param rep
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public DoubleList( StringBuffer rep ) throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag( rep, getClass().getSimpleName() );
		size = XMLParser.extractObjectForTags( xml, "size", int.class );
		array = new double[size];
		Map<String,String> filter = new TreeMap<String, String>();
		for( int i = 0; i < size; i++ ) { 
			filter.clear();
			filter.put( "val", "" + i );
			array[i] = XMLParser.extractObjectAndAttributesForTags( xml, "pos", null, filter, double.class );
		}
	}

	/**
	 * This method adds all elements of {@link DoubleList} <code>list2</code> to the current list.
	 * 
	 * @param list2 a list of elements to be added
	 */
	public void addAll(DoubleList list2){
		if(array.length <= size+list2.length()){
			double[] help = new double[2* (array.length+list2.length())];
			System.arraycopy( array, 0, help, 0, size );
			array = help;
		}
		System.arraycopy( list2.array, 0, array, size, list2.length() );
		size += list2.length();
	}
	
	/**
	 * Adds the element <code>val</code> at the end of the list.
	 * 
	 * @param val
	 *            the element that should be inserted
	 */
	public final void add( double val ) {
		if( array.length == size ) {
			help = new double[2 * array.length];
			System.arraycopy( array, 0, help, 0, size );
			array = help;
		}
		array[size++] = val;
	}

	/**
	 * Adds the element <code>val</code> from <code>fromIndex</code> to
	 * <code>toIndex</code> (exclusive).
	 * 
	 * @param val
	 *            the element that should be inserted
	 * @param fromIndex
	 *            the start index (inclusive)
	 * @param toIndex
	 *            the end index (exclusive)
	 */
	public void add( double val, int fromIndex, int toIndex ) {
		if( toIndex > array.length ) {
			double[] help = new double[toIndex];
			System.arraycopy( array, 0, help, 0, size );
			array = help;
		}
		while( fromIndex < toIndex ) {
			array[fromIndex++] = val;
		}
		if( size < toIndex ) {
			size = toIndex;
		}
	}

	/**
	 * Removes all elements from the list.
	 */
	public final void clear() {
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
	public final double get( int index ) {
		if( index < size ) {
			return array[index];
		} else {
			throw new ArrayIndexOutOfBoundsException( index );
		}
	}

	/**
	 * Returns the number of inserted elements.
	 * 
	 * @return the number of inserted elements
	 */
	public final int length() {
		return size;
	}

	/**
	 * This method returns a <code>double</code> array containing all elements
	 * of the list.
	 * 
	 * @return a <code>double</code> array containing all elements of the list
	 */
	public double[] toArray() {
		double[] erg = new double[size];
		System.arraycopy( array, 0, erg, 0, size );
		return erg;
	}

	/**
	 * Multiplies all values in the list from index <code>start</code> to
	 * <code>end</code> with the value <code>factor</code>.
	 * 
	 * @param start
	 *            the start index (inclusive)
	 * @param end
	 *            the end index (exclusive)
	 * @param factor
	 *            the factor with which the values in the list are multiplied
	 */
	public void multiply( int start, int end, double factor ) {
		while( start < end ) {
			array[start++] *= factor;
		}
	}
	
	/**
	 * Adds to all values in the list from index <code>start</code> to
	 * <code>end</code> the value <code>summand</code>.
	 * 
	 * @param start
	 *            the start index (inclusive)
	 * @param end
	 *            the end index (exclusive)
	 * @param summand
	 *            the summand which is added to the values in the list
	 */
	public void addTo( int start, int end, double summand ) {
		while( start < end ) {
			array[start++] += summand;
		}
	}

	/**
	 * This method computes the mean of a part of the list. Please note that
	 * <code>start</code> has to be smaller than <code>end</code> and
	 * <code>end</code> has to be smaller than {@link DoubleList#length()}.
	 * 
	 * @param start
	 *            the start index (inclusive)
	 * @param end
	 *            the end index (exclusive)
	 * 
	 * @return the mean of the part of the list
	 */
	public double mean( int start, int end ) {
		if( start >= end ) {
			throw new IllegalArgumentException( "The start index has to be smaller than the end index." );
		}
		double sum = 0;
		for( int k = start; k < end; k++ ) {
			sum += array[k];
		}
		return sum / (double)( end - start );
	}
	
	/**
	 * This method computes the minimum of a part of the list. Please note that
	 * <code>start</code> has to be smaller than <code>end</code> and
	 * <code>end</code> has to be smaller than {@link DoubleList#length()}.
	 * 
	 * @param start
	 *            the start index (inclusive)
	 * @param end
	 *            the end index (exclusive)
	 * 
	 * @return the minimum of the part of the list
	 */
	public double min( int start, int end ) {
		if( start >= end ) {
			throw new IllegalArgumentException( "The start index has to be smaller than the end index." );
		}
		double min = array[start];
		for( int k = start+1; k < end; k++ ) {
			if( array[k] < min ) {
				min = array[k];
			}
		}
		return min;
	}
	
	/**
	 * This method computes the median of a part of the list. Please note that
	 * <code>start</code> has to be smaller than <code>end</code> and
	 * <code>end</code> has to be smaller than {@link DoubleList#length()}.
	 * 
	 * @param start
	 *            the start index (inclusive)
	 * @param end
	 *            the end index (exclusive)
	 * 
	 * @return the median of the part of the list
	 */
	public double median( int start, int end ) {
		if( start >= end ) {
			throw new IllegalArgumentException( "The start index has to be smaller than the end index." );
		}
		double[] part = new double[end-start];
		System.arraycopy( array, start, part, 0, part.length );
		Arrays.sort( part );
		if(part.length % 2 == 0){
			return (part[part.length/2] + part[part.length/2 - 1])/2.0;
		}else{
			return part[(part.length-1)/2];
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( size * 30 );
		XMLParser.appendObjectWithTags( xml, size, "size" );
		for( int i = 0; i < size; i++ ) {
			XMLParser.appendObjectWithTagsAndAttributes( xml, array[i], "pos", "val=\"" + i + "\"" );
		}
		XMLParser.addTags( xml, getClass().getSimpleName() );
		return xml;
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
}
