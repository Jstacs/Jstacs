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
 * This class implements affine gap costs, i.e., the costs for starting a new gap are given by <code>start</code>, and
 * the costs for elongating a gap by one position are given by <code>elong</code>.
 * In this implementation, we may differentiate between the costs of opening an insert gap, i.e., a gap in the first string, 
 * and a delete gap, i.e., a gap in the second string.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class AffineCosts implements Costs {
	private double startInsert, startDelete;
	private double elongInsert, elongDelete;
	private boolean equal;
	private Costs c;
	
	/**
	 * This constructor creates a new instance of cost using affine gap costs.
	 * The costs for match, mismatch and gap elongation are defined by <code>c</code>,
	 * while the costs for a start of a gap are given by <code>start</code>.
	 * 
	 * @param start the costs for starting a gap
	 * @param c the cost for match, mismatch and gap elongation
	 */
	public AffineCosts( double start, Costs c ) {
		this(start, start, c);
	}
	
	/**
	 * This constructor creates a new instance of cost using affine gap costs.
	 * The costs for match, mismatch and gap elongation are defined by <code>c</code>,
	 * while the costs for a start of a gap are given by <code>startInsert</code> and <code>startDelete</code>, respectively.
	 * 
	 * @param startInsert the costs for starting an insert gap
	 * @param startDelete the costs for starting a delete gap
	 * @param c the cost for match, mismatch and gap elongation
	 */
	public AffineCosts( double startInsert, double startDelete, Costs c ) {
		this.c = c;
		this.elongInsert = c.getInsertCosts();
		this.elongDelete = c.getDeleteCosts();
		if( (startInsert < 0 && -startInsert > elongInsert)
				|| ( startDelete < 0 && -startDelete > elongDelete) ) {
			throw new IllegalArgumentException( "Problem: start < 0 && -start > elong" );
		}
		this.startInsert = startInsert;
		this.startDelete = startDelete;
		this.equal = elongInsert == elongDelete && startInsert == startDelete;
	}
	
	/**
	 * Restores {@link AffineCosts} object from its XML representation.
	 * @param xml the XML representation
	 * @throws NonParsableException if the XML could not be parsed
	 */
	public AffineCosts(StringBuffer xml) throws NonParsableException{
		xml = XMLParser.extractForTag( xml, "AffineCosts" );
		c = (Costs)XMLParser.extractObjectForTags( xml, "c" );
		try{
		elongInsert = (Double)XMLParser.extractObjectForTags( xml, "elongInsert" );
		startInsert = (Double)XMLParser.extractObjectForTags( xml, "startInsert" );
		elongDelete = (Double)XMLParser.extractObjectForTags( xml, "elongDelete" );
		startDelete = (Double)XMLParser.extractObjectForTags( xml, "startDelete" );
		}catch(NonParsableException ex){
			elongInsert = elongDelete = (Double)XMLParser.extractObjectForTags( xml, "elong" );
			startInsert = startDelete = (Double)XMLParser.extractObjectForTags( xml, "start" );
		}
	}
	
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, c, "c" );
		XMLParser.appendObjectWithTags( xml, elongInsert, "elongInsert" );
		XMLParser.appendObjectWithTags( xml, startInsert, "startInsert" );
		XMLParser.appendObjectWithTags( xml, elongDelete, "elongDelete" );
		XMLParser.appendObjectWithTags( xml, startDelete, "startDelete" );
		XMLParser.addTags( xml, "AffineCosts" );
		return xml;
	}
	
	/**
	 * Returns the internal costs (supplied as <code>c</code> to {@link #AffineCosts(double, Costs)}) used for matches and mismatches.
	 * @return the internal costs
	 */
	public Costs getInternalCosts(){
		return c;
	}
	
	/**
	 * Returns the costs for elongating an insert gap by one position.
	 * 
	 * @return the corresponding costs
	 */
	public double getElongateInsertCosts() {
		return elongInsert;
	}
	
	/**
	 * Returns the costs for elongating a delete gap by one position.
	 * 
	 * @return the corresponding costs
	 */
	public double getElongateDeleteCosts() {
		return elongDelete;
	}

	/**
	 * Returns the costs for an insert gap of length <code>length</code>.
	 * 
	 * @param length
	 *            the length of the gap
	 * 
	 * @return the corresponding costs
	 */
	public double getInsertCostsFor( int length ) {
		return startInsert + ( length * elongInsert );
	}
	
	/**
	 * Returns the costs for a delete gap of length <code>length</code>.
	 * 
	 * @param length
	 *            the length of the gap
	 * 
	 * @return the corresponding costs
	 */
	public double getDeleteCostsFor( int length ) {
		return startDelete + ( length * elongDelete );
	}
	
	/**
	 * Returns the costs for a gap of length <code>length</code>.
	 * 
	 * 
	 * @param length
	 *            the length of the gap
	 * 
	 * @return the corresponding costs
	 */
	public double getGapCostsFor( int length ) {
		if(equal){
			return startDelete + ( length * elongDelete );
		}else{
			throw new IllegalArgumentException("Costs for delete and insert not equal.");
		}
	}
	
	
	
	@Override
	public double getCostFor(Sequence s1, Sequence s2, int i, int j) {
		return c.getCostFor(s1, s2, i, j);
	}

	@Override
	public double getInsertCosts() {
		return c.getInsertCosts();
	}
	
	@Override
	public double getDeleteCosts() {
		return c.getDeleteCosts();
	}
}