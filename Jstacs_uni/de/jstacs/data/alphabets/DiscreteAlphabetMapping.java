/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.data.alphabets;

import de.jstacs.Storable;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * This class implements the transformation of discrete values to other discrete values
 * which define a {@link DiscreteAlphabet}.
 * 
 * @author Jens Keilwagen
 */
public class DiscreteAlphabetMapping implements Storable {

	/**
	 * The assignment of the old value (index) to the new value (entry of the array). 
	 */
	private int[] newValues;
	
	private double[] logNumSimilarSymbols;
	
	/**
	 * The new {@link de.jstacs.data.alphabets.Alphabet}.
	 */
	private DiscreteAlphabet newAlphabet;
	
	/**
	 * The main constructor creating a {@link DiscreteAlphabetMapping}.
	 * 
	 * @param newValues the assignment of the old value (index) to the new value (entry of the array)
	 * @param newAlphabet the new {@link de.jstacs.data.alphabets.Alphabet}
	 */
	public DiscreteAlphabetMapping( int[] newValues, DiscreteAlphabet newAlphabet ) {
		this.newValues = newValues.clone();
		this.newAlphabet = newAlphabet;
		precompute();
	}
	
	private static final String XML_TAG = "DiscreteAlphabetMapping";
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link DiscreteAlphabetMapping} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DiscreteAlphabetMapping} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public DiscreteAlphabetMapping( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, XML_TAG );
		newAlphabet = (DiscreteAlphabet) XMLParser.extractObjectForTags( xml, "newAlphabet" );
		newValues = XMLParser.extractObjectForTags( xml, "newValues", int[].class );
		precompute();
	}
	
	private void precompute() {
		int max = 0;
		for( int i = 0; i < newValues.length; i++ ) {
			max = Math.max( max, newValues[i] );
		}
		if( (max+1) != newAlphabet.length() ) {
			throw new IllegalArgumentException( "Please check the assignment from the original symbol/values to the new symbols/values." );
		}
		double[] originalSymbols = new double[max+1];
		for( int i = 0; i < newValues.length; i++ ) {
			originalSymbols[newValues[i]]++;
		}
		logNumSimilarSymbols = new double[newValues.length];
		for( int i = 0; i < newValues.length; i++ ) {
			logNumSimilarSymbols[i] = Math.log(originalSymbols[newValues[i]]);
		}
	}
	
	public StringBuffer toXML() {
		StringBuffer sb = new StringBuffer();
		XMLParser.appendObjectWithTags( sb, newAlphabet, "newAlphabet" );
		XMLParser.appendObjectWithTags( sb, newValues, "newValues" );
		XMLParser.addTags( sb, XML_TAG );
		return sb;
	}
	
	/**
	 * Returns the new Alphabet that is used for mapping.
	 * 
	 * @return the new Alphabet that is used for mapping.
	 */
	public DiscreteAlphabet getNewAlphabet() {
		return newAlphabet;
	}
	
	/**
	 * This method implements the main transformation of the values.
	 * 
	 * @param oldValue the old value
	 * 
	 * @return the new value
	 */
	public int getNewDiscreteValue( int oldValue ) {
		return newValues[oldValue];
	}
	
	/**
	 * This method returns the logarithm of the number of old values that yield the same new value.
	 * 
	 * @param oldValue the old value
	 * 
	 * @return the logarithm of the number of old values that yield the same new value
	 * 
	 * @see de.jstacs.data.sequences.MappedDiscreteSequence#getLogNumberOfPossibleOriginalSequences()
	 */
	public double getLogNumberOfSimilarSymbols( int oldValue ) {
		return logNumSimilarSymbols[oldValue];
	}
}
