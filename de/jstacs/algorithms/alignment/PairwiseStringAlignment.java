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
package de.jstacs.algorithms.alignment;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * Class for the representation of an alignment of two {@link String}s. It
 * contains the two {@link String}s that were aligned and expanded by
 * gap-symbols, and the edit-costs according to the employed
 * {@link de.jstacs.algorithms.alignment.cost.Costs} instance.
 * 
 * @author Jan Grau
 */
public class PairwiseStringAlignment extends StringAlignment {
	
	private int start, end, nummatches;

	/**
	 * Creates the instance for the two (extended) {@link String}s and the
	 * edit-costs.
	 * 
	 * @param r1
	 *            the first {@link String}
	 * @param r2
	 *            the second {@link String}
	 * @param cost
	 *            the edit-costs
	 * @param startPos
	 *            the start position of the aligned block in the first {@link String}
	 * @param endPos
	 *            the end position of the aligned block in the first {@link String}
	 * @param numMatches the number of matches in the alignment
	 */
	protected PairwiseStringAlignment( String r1, String r2, double cost, int startPos, int endPos, int numMatches ) {
		super( cost, r1, r2 );
		this.start = startPos;
		this.end = endPos;
		this.nummatches = numMatches;
	}
	
	
	/**
	 * @param xml
	 * @throws NonParsableException
	 */
	public PairwiseStringAlignment( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}


	/**
	 * Returns the number of matches in this alignment.
	 * @return the number of matches
	 */
	public int getNumberOfMatches(){
		return nummatches;
	}
	
	/**
	 * This method returns the start index of the alignment in the first sequence.
	 * 
	 * @return the start index of the alignment in the first sequence
	 */
	public int getStartIndexOfAlignmentForFirst() {
		return start;
	}
	
	/**
	 * This method returns the end index of the alignment in the first sequence.
	 * 
	 * @return the end index of the alignment in the first sequence
	 */
	public int getEndIndexOfAlignmentForFirst() {
		return end;
	}
	
	protected void fromXML(StringBuffer xml) throws NonParsableException{
		xml = XMLParser.extractForTag( xml, "PairwiseStringAlignment" );
		super.fromXML( xml );
		end = (Integer)XMLParser.extractObjectForTags( xml, "end" );
		nummatches = (Integer)XMLParser.extractObjectForTags( xml, "nummatches" );
		start = (Integer)XMLParser.extractObjectForTags( xml, "start" );
		
	}
	
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		super.toXML();
		XMLParser.appendObjectWithTags( xml, end, "end" );
		XMLParser.appendObjectWithTags( xml, nummatches, "nummatches" );
		XMLParser.appendObjectWithTags( xml, start, "start" );
		XMLParser.addTags( xml, "PairwiseStringAlignment" );
		return xml;
	}
	
}
