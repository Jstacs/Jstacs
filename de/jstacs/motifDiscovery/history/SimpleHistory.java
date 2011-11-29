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

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * This class implements a simple history that has a limited memory that will be
 * used cyclicly. Only operations are allowed, that are not a priorily forbidden
 * and the complementary operation has not been performed before (in the local
 * history) except shrink operations. The last point helps to keep the motifs short
 * that will be found. 
 * 
 * <br>
 * 
 * <b>This history may not cause an termination of the motif discovery.</b>
 * 
 * @author Jens Keilwagen
 * 
 * @see CappedHistory
 */
public class SimpleHistory implements History {

	/**
	 * Local history.
	 */
	private int[][] memory;

	/**
	 * The pointer.
	 */
	private int index;

	/**
	 * Switches for a priori forbidden operations.
	 */
	private boolean allowShift, allowShrink, allowExpand;

	/**
	 * This constructor creates a simple history with limited memory.
	 * All operation are a priori allowed.
	 * 
	 * @param slots
	 *            the number of memory slots
	 *            
	 * @see SimpleHistory#SimpleHistory(int, boolean, boolean, boolean)
	 */
	public SimpleHistory( int slots ) {
		this( slots, true, true, true );
	}

	/**
	 * This constructor creates a simple history with limited memory.
	 * 
	 * @param slots
	 *            the number of memory slots
	 * @param allowShift
	 *            if <tt>true</tt> shifts are allow
	 * @param allowShrink
	 *            if <tt>true</tt> shrinks are allow
	 * @param allowExpand
	 *            if <tt>true</tt> expands are allow
	 */
	public SimpleHistory( int slots, boolean allowShift, boolean allowShrink, boolean allowExpand ) {
		memory = new int[slots][];
		index = 0;
		this.allowShift = allowShift;
		this.allowShrink = allowShrink;
		this.allowExpand = allowExpand;
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
	public SimpleHistory( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getXMLTag() );
		memory = (int[][]) XMLParser.extractObjectForTags( xml, "memory" );
		index = (Integer) XMLParser.extractObjectForTags( xml, "index" );
		allowExpand = (Boolean) XMLParser.extractObjectForTags( xml, "allowExpand" );
		allowShift = (Boolean) XMLParser.extractObjectForTags( xml, "allowShift" );
		allowShrink = (Boolean) XMLParser.extractObjectForTags( xml, "allowShrink" );
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
		XMLParser.appendObjectWithTags( xml, memory, "memory" );
		XMLParser.appendObjectWithTags( xml, index, "index" );
		XMLParser.appendObjectWithTags( xml, allowExpand, "allowExpand" );
		XMLParser.appendObjectWithTags( xml, allowShift, "allowShift" );
		XMLParser.appendObjectWithTags( xml, allowShrink, "allowShrink" );
		XMLParser.addTags( xml, getXMLTag() );
		return xml;
	}

	public boolean operationAllowed( int... operation ) {
		if( operation.length != 2 ) {
			return false;
		}
		int sum = operation[1] - operation[0];
		if( sum < 0 ) { //shrink operation
			return allowShrink;
		} else if( (!allowShift && sum == 0) //shift forbidden, but operation is a shift
				|| (!allowExpand && sum > 0) //expand forbidden, but operation is an expand
			)
		{
			return false;
		}
		for( int j, i = 0; i < memory.length; i++ ) {
			if( memory[i] != null && operation.length == memory[i].length ) {
				j = 0;
				while( j < memory[i].length && memory[i][j] == -operation[j] ) {
					j++;
				}
				if( j == memory[i].length ) {
					return false;
				}
			}
		}
		return true;
	}

	public void operationPerfomed( int... operation ) {
		if( memory.length != 0 ) {
			memory[index] = operation;
			index++;
			index %= memory.length;
		}
	}

	public void clear() {
		Arrays.fill( memory, null );
		index = 0;
	}

	public SimpleHistory clone() throws CloneNotSupportedException {
		SimpleHistory clone = (SimpleHistory)super.clone();
		clone.memory = new int[memory.length][];
		for( int i = 0; i < memory.length; i++ ) {
			if( memory[i] != null ) {
				clone.memory[i] = memory[i].clone();
			}
		}
		return clone;
	}
}
