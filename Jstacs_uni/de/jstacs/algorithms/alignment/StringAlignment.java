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

import java.util.List;

import de.jstacs.Storable;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.Result;

/**
 * Class for the representation of an alignment of {@link String}s. It
 * contains the {@link String}s that were aligned and expanded by
 * gap-symbols, and the edit-costs according to the employed
 * {@link de.jstacs.algorithms.alignment.cost.Costs} instance.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class StringAlignment implements Comparable<StringAlignment>, Storable{
	private String[] r;
	private double cost;
	private Result res;

	/**
	 * Restores {@link StringAlignment} object from its XML representation.
	 * @param xml the XML representation
	 * @throws NonParsableException if the XML could not be parsed
	 */
	public StringAlignment(StringBuffer xml) throws NonParsableException{
		fromXML( xml );
	}
	
	/**
	 * This constructor creates an instance storing the aligned Strings and the costs of the alignment.
	 * 
	 * @param cost the cost of the alignment
	 * @param strings the aligned and expanded Strings
	 */
	public StringAlignment(double cost, String... strings) {
		this(cost, strings, (Result)null);
	}
	
	/**
	 * This constructor creates an instance storing the aligned Strings and the costs of the alignment.
	 * 
	 * @param cost the cost of the alignment
	 * @param strings the aligned and expanded Strings
	 * @param res a result describing the alignment
	 */
	public StringAlignment(double cost, String[] strings, Result res) {
		this.cost = cost;
		r = strings.clone();
		this.res = res;
	}
	
	/**
	 * Returns the annotation result.
	 * 
	 * @return the annotation result
	 */
	public Result getAnnotationResult() {
		return res;
	}

	/**
	 * Returns the costs.
	 * 
	 * @return the costs
	 */
	public double getCost() {
		return cost;
	}

	/**
	 * Returns the number of sequences in this alignment.
	 * 
	 * @return the number of sequences in this alignment.
	 */
	public int getNumberOfAlignedSequences() {
		return r.length;
	}
	
	/**
	 * Returns the aligned {@link String} with index <code>index</code>.
	 * 
	 * @param index
	 *            the index of the {@link String}
	 * 
	 * @return the aligned {@link String} with index <code>index</code>
	 */
	public String getAlignedString(int index) {
		return r[index];
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return toString( Integer.MAX_VALUE, -1, null );
	}
	
	/**
	 * This method returns a String representation of the alignment with a given chunk size.
	 * I.e. all aligned sequence will be split each <code>chunkSize</code> symbols.
	 * 
	 * @param chunkSize the size of the chunks
	 * @param simplifyIdx if non-negative the alignment is given relative to the sequence with the corresponding index; i.e. match is represented by &quot;.&quot; and mismatch by the corresponding symbol
	 * @param annot a list of annotations for the aligned Strings, can be <code>null</code>
	 * 
	 * @return a String representation of the alignment
	 */
	public String toString( int chunkSize, int simplifyIdx, List<String> annot ) {
		StringBuffer buf = new StringBuffer();
		
		String[] d;
		if( simplifyIdx >= 0 ) {
			StringBuffer[] data = new StringBuffer[r.length];
			int[] start = new int[r.length], end = new int[r.length];
			for( int i = 0; i < r.length; i++ ) {
				data[i] = new StringBuffer();
				
				int l = 0; 
				while( l < r[i].length() && r[i].charAt( l ) == '-' ) {
					l++;
				}
				start[i] = l;
				
				l = r[i].length()-1;
				while( l >= 0 && r[i].charAt( l ) == '-' ) {
					l--;
				}
				end[i] = l;
			}
			
			char same = '.', border = ' ';
			for( int l = 0; l < r[simplifyIdx].length(); l++ ) {
				for( int i = 0; i < r.length; i++ ) {
					data[i].append( i==simplifyIdx ? r[i].charAt( l ) : (start[i] <= l && l <= end[i] ? (r[i].charAt(l)==r[simplifyIdx].charAt(l) ? same : r[i].charAt(l)) : border) );
				}
			}
			
			d = new String[r.length];
			for( int i = 0; i < r.length; i++ ) {
				d[i] = data[i].toString();
			}
		} else {
			d = r;
		}
		
		for( int end, start = 0; start < d[0].length(); start += chunkSize ) {
			end = Math.min( start+chunkSize, d[0].length() );
			for (int i = 0; i < r.length; i++) {
				buf.append( d[i].substring( start, end ) );
				if( annot != null && i < annot.size() ) {
					buf.append( "\t" );
					buf.append( annot.get(i) );
				}
				buf.append("\n");
			}
			buf.append("\n");
		}
		buf.append("costs: ");
		buf.append(cost);
		return buf.toString();
	}

	/**
	 * This method return the length of the alignment.
	 * 
	 * @return the length of the alignment
	 */
	public int getLength() {
		return r[0].length();
	}

	@Override
	public int compareTo(StringAlignment o) {
		if( this.getNumberOfAlignedSequences() == o.getNumberOfAlignedSequences() ) {
			int n = getNumberOfAlignedSequences(), c = 0, i = 0;
			while( i < n && (c=r[i].compareTo(o.r[i])) == 0 ) {
				i++;
			}
			return c;
		} else {
			return this.getNumberOfAlignedSequences() - o.getNumberOfAlignedSequences();
		}
	}

	/**
	 * Parses the XML representation.
	 * @param xml the XML representation
	 * @throws NonParsableException if XML could not be parsed
	 */
	protected void fromXML(StringBuffer xml) throws NonParsableException{
		xml = XMLParser.extractForTag( xml, "StringAlignment" );
		
		cost = (Double)XMLParser.extractObjectForTags( xml, "cost" );
		r = (String[])XMLParser.extractObjectForTags( xml, "r" );
		res = (Result)XMLParser.extractObjectForTags( xml, "res" );
		
	}
	
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, cost, "cost" );
		XMLParser.appendObjectWithTags( xml, r, "r" );
		XMLParser.appendObjectWithTags( xml, res, "res" );
		XMLParser.addTags( xml, "StringAlignment" );
		return xml;
	}
}