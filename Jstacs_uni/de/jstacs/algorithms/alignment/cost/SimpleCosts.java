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
package de.jstacs.algorithms.alignment.cost;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * Class for simple costs with costs 
 * <ul>
 * <li><code>match</code> for a match,</li>
 * <li><code>mismatch</code> for a mismatch, and</li>
 * <li><code>gap</code> for a gap (of length 1).</li>
 * </ul>
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class SimpleCosts implements Costs {

	private double match;
	private double mismatch;
	private double gap;

	/**
	 * Creates a new instance of simple costs with costs
	 * <code>match</code> for a match,
	 * <code>mismatch</code> for a mismatch, and
	 * <code>gap</code> for a gap (of length 1).
	 *
	 * @param match
	 *            the match costs
	 * @param mismatch
	 *            the mismatch costs
	 * @param gap
	 *            the costs for a gap
	 */
	public SimpleCosts( double match, double mismatch, double gap ) {
		this.gap = gap;
		if( mismatch <= 0 ) {
			throw new IllegalArgumentException( "Problem: mismatch <= 0" );
		}
		this.mismatch = mismatch;
		if( match > 0 ) {
			throw new IllegalArgumentException( "Problem: match > 0" );
		}
		this.match = match;
	}
	
	/**
	 * Restores {@link SimpleCosts} object from its XML representation.
	 * @param xml the XML representation
	 * @throws NonParsableException if the XML could not be parsed
	 */
	public SimpleCosts( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, "SimpleCosts" );
		gap = (Double)XMLParser.extractObjectForTags( xml, "gap" );
		match = (Double)XMLParser.extractObjectForTags( xml, "match" );
		mismatch = (Double)XMLParser.extractObjectForTags( xml, "mismatch" );
		
	}

	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, gap, "gap" );
		XMLParser.appendObjectWithTags( xml, match, "match" );
		XMLParser.appendObjectWithTags( xml, mismatch, "mismatch" );
		XMLParser.addTags( xml, "SimpleCosts" );
		return xml;
	}
	
	public double getCostFor( Sequence s1, Sequence s2, int i, int j ) {
		return s1.discreteVal( i - 1 ) != s2.discreteVal( j - 1 ) ? mismatch : match;
	}

	public double getGapCosts() {
		return gap;
	}
}
