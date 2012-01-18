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

import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * This class implements a history that allows operations, that are not a
 * priorily forbidden and do not create a configuration that has already be
 * considered. If not a priorily forbidden shrink operations are always allowed.
 * 
 * <br>
 * 
 * <b>This history may not cause an termination of the motif discovery.</b>
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see CappedHistory
 */
public class NoRevertHistory implements History {

	private int[] current;

	private LinkedList<int[]> list;

	/**
	 * Switches for a priori forbidden operations.
	 */
	private boolean allowShift, allowShrink, allowExpand;

	//private boolean next;
	
	/**
	 * This constructor creates an instance that allows to shift shrink and expand the motif.
	 */
	public NoRevertHistory() {
		this( true, true, true );
	}

	/**
	 * This constructor creates an instance with user specified allowed operations.
	 * 
	 * @param allowShift whether it is allowed to shift the motif
	 * @param allowShrink whether it is allowed to shrink the motif
	 * @param allowExpand whether it is allowed to expand the motif
	 */
	public NoRevertHistory( boolean allowShift, boolean allowShrink, boolean allowExpand ) {
		this.allowShift = allowShift;
		this.allowShrink = allowShrink;
		this.allowExpand = allowExpand;
		list = new LinkedList<int[]>();
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
	public NoRevertHistory( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getXMLTag() );
		current = (int[]) XMLParser.extractObjectForTags( xml, "current" );
		int[][] help = (int[][]) XMLParser.extractObjectForTags( xml, "list" );
		list = new LinkedList<int[]>();
		if( help != null ) {
			for( int i = 0; i < help.length; i++ ) {
				list.add( help[i] );
			}
		}
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
		XMLParser.appendObjectWithTags( xml, current, "current" );
		XMLParser.appendObjectWithTags( xml, list.toArray( new int[0][] ), "list" );
		XMLParser.appendObjectWithTags( xml, allowExpand, "allowExpand" );
		XMLParser.appendObjectWithTags( xml, allowShift, "allowShift" );
		XMLParser.appendObjectWithTags( xml, allowShrink, "allowShrink" );
		XMLParser.addTags( xml, getXMLTag() );
		return xml;
	}

	public NoRevertHistory clone() throws CloneNotSupportedException {
		NoRevertHistory clone = (NoRevertHistory)super.clone();
		clone.current = current.clone();
		clone.list = new LinkedList<int[]>();
		Iterator<int[]> it = list.iterator();
		while( it.hasNext() ) {
			clone.list.add( it.next().clone() );
		}
		return clone;
	}

	public void clear() {
		list.clear();
		current = new int[]{ 0, 0 };
		list.add( current.clone() );
		//next = true;
	}

	public boolean operationAllowed( int... op ) {
		if( op.length != 2 ) {
			return false;
		}
		int sum = op[1] - op[0];
		if( sum < 0 ) { //shrink operation
			return allowShrink;
		} else if( (!allowShift && sum == 0) //shift forbidden, but operation is a shift
				|| (!allowExpand && sum > 0) //expand forbidden, but operation is an expand
				//|| (!next && sum >= 0)
			)
		{
			return false;
		}
		Iterator<int[]> it = list.iterator();
		int[] temp = { current[0] + op[0], current[1] + op[1] };
		while( it.hasNext() ) {
			int[] temp2 = it.next();
			if( temp[0] == temp2[0] && temp[1] == temp2[1] ) {
				//we have been in this situation before ...
				/*
				if( sum < 0 ) {
					//next = false;
					return true;
				} else {
					return false;
				}*/
				return false;
			}
		}
		return true;
	}

	public void operationPerfomed( int... op ) {
		current[0] += op[0];
		current[1] += op[1];
		//System.out.println( Arrays.toString( current ) );
		list.add( current.clone() );
	}
}