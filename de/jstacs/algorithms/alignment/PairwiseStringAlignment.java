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

/**
 * Class for the representation of an alignment of two {@link String}s. It
 * contains the two {@link String}s that were aligned and expanded by
 * gap-symbols, and the edit-costs according to the employed
 * {@link de.jstacs.algorithms.alignment.cost.Costs} instance.
 * 
 * @author Jan Grau
 */
public class PairwiseStringAlignment extends StringAlignment {
	
	private int start, end;

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
	 */
	protected PairwiseStringAlignment( String r1, String r2, double cost, int startPos, int endPos ) {
		super( cost, r1, r2 );
		this.start = startPos;
		this.end = endPos;
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
}
