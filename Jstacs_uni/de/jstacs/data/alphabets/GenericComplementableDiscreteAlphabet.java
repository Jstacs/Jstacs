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
package de.jstacs.data.alphabets;

import java.util.StringTokenizer;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.SimpleParameter;

/**
 * This class implements an generic complementable discrete alphabet. The complement is determined in the constructor, where the user has to specify for each symbol the index of the complement symbol.
 * 
 * @author Jens Keilwagen
 */
public class GenericComplementableDiscreteAlphabet extends
		ComplementableDiscreteAlphabet {

	private int[] comp;
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link GenericComplementableDiscreteAlphabet} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link GenericComplementableDiscreteAlphabet} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public GenericComplementableDiscreteAlphabet( StringBuffer representation ) throws NonParsableException {
		super(representation);
	}

	/**
	 * This constructor creates a {@link GenericComplementableDiscreteAlphabet} from a parameter set.
	 * 
	 * @param parameters the set of parameters needed to create a {@link GenericComplementableDiscreteAlphabet}
	 * 
	 * @throws IllegalArgumentException
	 *             if space or tab will be used as symbols
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 */
	public GenericComplementableDiscreteAlphabet( GenericComplementableDiscreteAlphabetParameterSet parameters ) throws IllegalArgumentException, DoubleSymbolException {
		super( parameters );
		StringTokenizer tok = new StringTokenizer( (String)parameters.getParameterAt( 2 ).getValue(), " " );
		comp = new int[tok.countTokens()];
		int i = 0;
		while( tok.hasMoreTokens() ) {
			comp[i++] = Integer.parseInt( tok.nextToken() );
		}
	}
	
	/**
	 * The main constructor.
	 *  
	 * @param alphabet
	 *            the array of symbols
	 * @param caseInsensitive
	 *            indicates if there will be no difference between lowercase and
	 *            uppercase letters/symbols in the alphabet (no case sensitivity)
	 * @param revComp the array is used to provide the complementarity, <code>revComp[i]</code> contains the index of the complement of the i-th symbol
	 * 
	 * @throws IllegalArgumentException
	 *             if space or tab will be used as symbols
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 */
	public GenericComplementableDiscreteAlphabet( boolean caseInsensitive, String[] alphabet, int[] revComp ) throws DoubleSymbolException,
			IllegalArgumentException {
		super( caseInsensitive, alphabet );
		this.comp = revComp.clone();
	}

	@Override
	public int getComplementaryCode(int code) {
		return comp[code];
	}

	@Override
	public AlphabetParameterSet getCurrentParameterSet() throws Exception {
		if( parameters != null ) {
			return parameters.clone();
		} else {
			return new GenericComplementableDiscreteAlphabetParameterSet( alphabet.clone(), caseInsensitive, comp.clone() );
		}
	}
	
	/**
	 * This class is used as container for the parameters of a {@link GenericComplementableDiscreteAlphabet}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class GenericComplementableDiscreteAlphabetParameterSet extends DiscreteAlphabetParameterSet {
		
		/**
		 * This constructor creates an empty parameter set the has to be filled before it can be used to create a {@link GenericComplementableDiscreteAlphabet}.
		 */
		public GenericComplementableDiscreteAlphabetParameterSet() {
			super( GenericComplementableDiscreteAlphabet.class );
		}
		
		/**
		 * The main constructor.
		 *  
		 * @param alphabet
		 *            the array of symbols
		 * @param caseInsensitive
		 *            indicates if there will be no difference between lowercase and
		 *            uppercase letters/symbols in the alphabet (no case sensitivity)
		 * @param revComp the array is used to provide the complementarity, <code>revComp[i]</code> contains the index of the complement of the i-th symbol
		 * 
		 * @throws Exception if some parameter could not be set.
		 */
		public GenericComplementableDiscreteAlphabetParameterSet( String[] alphabet, boolean caseInsensitive, int[] revComp ) throws Exception {
			this();
			
			loadParameters();
			String alphString = "" + alphabet[0], revString = "" + revComp[0];
			for( int i = 1; i < alphabet.length; i++ ) {
				alphString += " " + alphabet[i];
				revString += " " + revComp[i];
			}
			parameters.get( 0 ).setValue( alphString );
			parameters.get( 1 ).setValue( new String( caseInsensitive ? "Case insensitive" : "Case sensitive" ) );
			parameters.get( 2 ).setValue( revString );
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Creates a new {@link GenericComplementableDiscreteAlphabet.GenericComplementableDiscreteAlphabetParameterSet} out of its XML
		 * representation.
		 * 
		 * @param representation
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link GenericComplementableDiscreteAlphabet.GenericComplementableDiscreteAlphabetParameterSet} could not be
		 *             reconstructed out of the XML representation (the
		 *             {@link StringBuffer} <code>representation</code> could
		 *             not be parsed)
		 * 
		 * @see de.jstacs.data.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(StringBuffer) Alphabet.AlphabetParameterSet#AlphabetParameterSet(StringBuffer)
		 * @see de.jstacs.Storable
		 */
		public GenericComplementableDiscreteAlphabetParameterSet( StringBuffer representation ) throws NonParsableException {
			super( representation );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.ParameterSet#loadParameters()
		 */
		@Override
		protected void loadParameters() throws Exception {
			super.loadParameters();
			parameters.add( new SimpleParameter( DataType.STRING,
					"Values of the index for computings the reverse complement",
					"",
					true ) );
		}
	}
}