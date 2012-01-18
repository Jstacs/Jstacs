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

package de.jstacs.motifDiscovery.history;

import java.util.Arrays;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * This class implements a history that allows operations (i.e. a pair of
 * <code>int</code>), that are not a priorily forbidden and that are done before
 * less than a specified threshold. This enables the user for instance to allow
 * each operation exactly once. Furthermore, this class has a switch that allows
 * to forbid an operation for which the reverse operation has been performed in the
 * last step. 
 * 
 * <br><br>
 * 
 * If shrink operations are allowed by the user then a shrink operation is allow anytime.
 * This helps to keep the motif short.
 * 
 * @author Jens Keilwagen
 */
public class RestrictedRepeatHistory implements History {

	private int threshold;

	private int[] anz, last;

	private Hashtable<String, int[]> hash;

	/**
	 * Switches for a priori forbidden operations.
	 */
	private boolean allowShift, allowShrink, allowExpand, allowReverse;

	/**
	 * This constructor creates an instance that allows to shift, shrink, and
	 * expand the motif, to do the reverse operation, and to do each operation at most once.
	 */
	public RestrictedRepeatHistory() {
		this( true, true, true, true, 1 );
	}
	
	/**
	 * This constructor creates an instance that allows to shift, shrink, and
	 * expand the motif, and allows to do each operation at most once.
	 * 
	 * @param allowReverse
	 *            whether it is allowed to do the reverse operation as next operation or not
	 */
	public RestrictedRepeatHistory( boolean allowReverse ) {
		this( true, true, true, allowReverse, 1 );
	}

	/**
	 * This constructor creates an instance that allows to shift shrink and
	 * expand the motif.
	 * 
	 * @param threshold
	 *            the number of times each operation can be done
	 * 
	 * @throws IllegalArgumentException
	 *             if the <code>threshold</code> is below 1
	 */
	public RestrictedRepeatHistory( int threshold ) {
		this( true, true, true, true, threshold );
	}

	/**
	 * This constructor creates an instance that allows to do each operation at
	 * most once.
	 * 
	 * @param allowShift
	 *            whether it is allowed to shift the motif
	 * @param allowShrink
	 *            whether it is allowed to shrink the motif
	 * @param allowExpand
	 *            whether it is allowed to expand the motif
	 */
	public RestrictedRepeatHistory( boolean allowShift, boolean allowShrink, boolean allowExpand ) {
		this( allowShift, allowShrink, allowExpand, true, 1 );
	}

	/**
	 * This constructor creates an instance with user specified allowed
	 * operations and <code>threshold</code>.
	 * 
	 * @param allowShift
	 *            whether it is allowed to shift the motif
	 * @param allowShrink
	 *            whether it is allowed to shrink the motif
	 * @param allowExpand
	 *            whether it is allowed to expand the motif
	 * @param allowReverse
	 *            whether it is allowed to do the reverse operation as next operation or not
	 * @param threshold
	 *            the number of times each operation can be done
	 * 
	 * @throws IllegalArgumentException
	 *             if the <code>threshold</code> is below 1
	 */
	public RestrictedRepeatHistory( boolean allowShift, boolean allowShrink, boolean allowExpand, boolean allowReverse, int threshold ) throws IllegalArgumentException {
		this.allowShift = allowShift;
		this.allowShrink = allowShrink;
		this.allowExpand = allowExpand;
		this.allowReverse = allowReverse;
		if( threshold < 1 ) {
			throw new IllegalArgumentException( "The threshold has to be at least 1." );
		}
		this.threshold = threshold;
		hash = new Hashtable<String, int[]>();
		last = new int[2];
		clear();
	}
	
	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public RestrictedRepeatHistory( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getXMLTag() );
		last = (int[]) XMLParser.extractObjectForTags( xml, "last" );
		anz = (int[]) XMLParser.extractObjectForTags( xml, "anz" );
		
		String[] opList = (String[]) XMLParser.extractObjectForTags( xml, "opList" );
		int[][] anzList = (int[][]) XMLParser.extractObjectForTags( xml, "anzList" );
		hash = new Hashtable<String, int[]>();
		for( int i = 0; i < opList.length; i++ ) {
			hash.put( opList[i], anzList[i] );
		}
		
		threshold = (Integer) XMLParser.extractObjectForTags( xml, "threshold" );
		allowExpand = (Boolean) XMLParser.extractObjectForTags( xml, "allowExpand" );
		allowShift = (Boolean) XMLParser.extractObjectForTags( xml, "allowShift" );
		allowShrink = (Boolean) XMLParser.extractObjectForTags( xml, "allowShrink" );
		allowReverse = (Boolean) XMLParser.extractObjectForTags( xml, "allowReverse" );
	}
	
	private String getXMLTag() {
		return getClass().getSimpleName();
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, last, "last" );
		XMLParser.appendObjectWithTags( xml, anz, "anz" );
		
		Iterator<Entry<String,int[]>> it = hash.entrySet().iterator();
		Entry<String,int[]> e;
		LinkedList<String> opList = new LinkedList<String>();
		LinkedList<int[]> anzList = new LinkedList<int[]>();
		while( it.hasNext() ) {
			e = it.next();
			opList.add( e.getKey() );
			anzList.add( e.getValue() );
		}
		XMLParser.appendObjectWithTags( xml, opList.toArray( new String[0] ), "opList" );
		XMLParser.appendObjectWithTags( xml, anzList.toArray( new int[0][] ), "anzList" );
		
		XMLParser.appendObjectWithTags( xml, threshold, "threshold" );
		XMLParser.appendObjectWithTags( xml, allowExpand, "allowExpand" );
		XMLParser.appendObjectWithTags( xml, allowShift, "allowShift" );
		XMLParser.appendObjectWithTags( xml, allowShrink, "allowShrink" );
		XMLParser.appendObjectWithTags( xml, allowReverse, "allowReverse" );
		XMLParser.addTags( xml, getXMLTag() );
		return xml;
	}

	public RestrictedRepeatHistory clone() throws CloneNotSupportedException {
		RestrictedRepeatHistory clone = (RestrictedRepeatHistory)super.clone();
		clone.anz = anz==null?null:anz.clone();
		clone.hash = new Hashtable<String, int[]>();
		Iterator<Entry<String, int[]>> it = hash.entrySet().iterator();
		Entry<String, int[]> e;
		while( it.hasNext() ) {
			e = it.next();
			clone.hash.put( e.getKey(), e.getValue().clone() );
		}
		clone.last = last.clone();
		return clone;
	}

	public void clear() {
		hash.clear();
	}

	public boolean operationAllowed( int... op ) {
		if( ( op.length != 2 ) // not correct length 
			|| ( !allowReverse && op[0] == -last[0] && op[1] == -last[1] ) ) { // reverse operation
			return false;
		}
		int sum = op[1] - op[0];
		if( ( !allowShift && sum == 0 ) || ( !allowShrink && sum < 0 ) || ( !allowExpand && sum > 0 ) ) {
			return false;
		} else if ( allowShrink && sum < 0 ) {
			return true;
		} else {
			anz = hash.get( Arrays.toString( op ) );
			return( anz == null || anz[0] < threshold );
		}
	}

	public void operationPerfomed( int... op ) {
		String s = Arrays.toString( op );
		anz = hash.get( s );
		if( anz == null ) {
			anz = new int[]{ 1 };
		} else {
			anz[0]++;
		}
		hash.put( s, anz );
		last[0] = op[0];
		last[1] = op[1];
	}
}