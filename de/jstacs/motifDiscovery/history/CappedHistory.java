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

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * This class combines a threshold on the number of steps which can be performed with any other {@link History}.
 * Using a threshold ensures to converge using any other history. 
 * 
 * @author Jens Keilwagen
 * 
 */
public class CappedHistory implements History {

	private int threshold, number;
	private History h;
	
	/**
	 * This constructor creates an instance that allows at most <code>t</code> steps using the history <code>h</code>.
	 * 
	 * @param t the maximal number of steps/operations
	 * @param h the underlying history
	 * 
	 * @throws CloneNotSupportedException is the history <code>h</code> could not be cloned
	 */
	public CappedHistory( int t, History h ) throws CloneNotSupportedException {
		super();
		if( t < 0 ) {
			throw new IllegalArgumentException( "The threshold has to be non-negative." );
		}
		threshold = t;
		this.h = h.clone();
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
	public CappedHistory( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getXMLTag() );
		threshold = (Integer) XMLParser.extractObjectForTags( xml, "threshold" );
		number = (Integer) XMLParser.extractObjectForTags( xml, "number" );
		h = (History) XMLParser.extractObjectForTags( xml, "history" );
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
		XMLParser.appendObjectWithTags( xml, threshold, "threshold" );
		XMLParser.appendObjectWithTags( xml, number, "number" );
		XMLParser.appendObjectWithTags( xml, h, "history" );
		XMLParser.addTags( xml, getXMLTag() );
		return xml;
	}
	
	public CappedHistory clone() throws CloneNotSupportedException {
		CappedHistory clone = (CappedHistory) super.clone();
		clone.h = h.clone();
		return clone;
	}

	public void clear() {
		number = 0;
		h.clear();		
	}
	
	public boolean operationAllowed( int... op ) {
		if( number < threshold ) {
			return h.operationAllowed( op );
		} else {
			return false;
		}
	}

	public void operationPerfomed( int... op ) {
		number++;
		if( number <= threshold ) {
			h.operationPerfomed( op );
		}
	}
}