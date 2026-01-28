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

import java.util.Arrays;
import java.util.Hashtable;
import java.util.StringTokenizer;

import de.jstacs.DataType;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;

/**
 * Class for an alphabet that consists of arbitrary {@link String}s. For DNA
 * alphabets, the class {@link DNAAlphabet} should be used.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class DiscreteAlphabet extends Alphabet {

	private static final String XML_TAG = "DiscreteAlphabet";

	/**
	 * The alphabet as {@link String} array.
	 */
	protected String[] alphabet;

	/**
	 * For encoding.
	 */
	private Hashtable<String, Integer> hash;

	/**
	 * Switch to decide whether the input should be treated case sensitive or insensitive.
	 */
	protected boolean caseInsensitive;

	private int longestCharacter;

	/**
	 * The parameter set describing this {@link DiscreteAlphabet} .
	 */
	protected AlphabetParameterSet parameters;

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link DiscreteAlphabet} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DiscreteAlphabet} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer}
	 *             <code>representation</code> could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public DiscreteAlphabet( StringBuffer representation ) throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag( representation, XML_TAG );
		try {
			init( XMLParser.extractObjectForTags( xml, "symbols", String[].class ), XMLParser.extractObjectForTags( xml, "caseInsensitive", boolean.class ) );
			parseFromXML(xml);
		} catch ( Exception e ) {
			NonParsableException n = new NonParsableException( e.getMessage() );
			n.setStackTrace( e.getStackTrace() );
			throw n;
		}
	}

	/**
	 * The constructor for the {@link de.jstacs.InstantiableFromParameterSet}
	 * interface. Creates a new {@link DiscreteAlphabet} from a given set of
	 * parameters.
	 * 
	 * @param parameters
	 *            the parameter set for the {@link DiscreteAlphabet}
	 * 
	 * @throws IllegalArgumentException
	 *             if space or tab will be used as symbols
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 * 
	 * @see de.jstacs.InstantiableFromParameterSet
	 */
	public DiscreteAlphabet( DiscreteAlphabetParameterSet parameters ) throws IllegalArgumentException, DoubleSymbolException {
		try {
			this.parameters = parameters.clone();
			String alphValue = (String)parameters.getParameterAt( 0 ).getValue();
			boolean caseInsensitive = (Boolean)parameters.getParameterAt( 1 ).getValue();

			StringTokenizer tok = new StringTokenizer( alphValue, " " );
			String[] vals = new String[tok.countTokens()];
			int i = 0;
			while( tok.hasMoreTokens() ) {
				vals[i++] = tok.nextToken();
			}
			init( vals, caseInsensitive );
		} catch ( CloneNotSupportedException e ) {
			throw new IllegalArgumentException( e.getCause().getMessage() );
		}
	}
	
	protected void parseFromXML(StringBuffer xml) throws NonParsableException{
		
	}

	/* (non-Javadoc)
	 * @see de.jstacs.InstantiableFromParameterSet#getCurrentParameterSet()
	 */
	public AlphabetParameterSet getCurrentParameterSet() throws Exception {
		if( parameters != null ) {
			return parameters.clone();
		} else {
			return new DiscreteAlphabetParameterSet( alphabet.clone(), caseInsensitive );
		}
	}

	private void init( String[] alphabet, boolean caseInsensitive ) throws IllegalArgumentException, DoubleSymbolException {
		hash = new Hashtable<String, Integer>( alphabet.length, 1.0f );
		this.alphabet = new String[alphabet.length];
		this.caseInsensitive = caseInsensitive;
		longestCharacter = 0;
		for( int i = 0; i < alphabet.length; i++ ) {
			if( alphabet[i].length() == 0 ) {
				throw new IllegalArgumentException( "\"\" can not be a symbol/character." );
			}
			if( alphabet[i].indexOf( " " ) >= 0 || alphabet[i].indexOf( "\t" ) >= 0 ) {
				throw new IllegalArgumentException( "blanks and tabs can not be part of a symbol/character: \'" + alphabet[i] + "\'" );
			}
			if( hash.containsKey( caseInsensitive ? alphabet[i].toUpperCase() : alphabet[i] ) ) {
				throw new DoubleSymbolException( caseInsensitive ? alphabet[i].toUpperCase() : alphabet[i] );
			}
			hash.put( caseInsensitive ? alphabet[i].toUpperCase() : alphabet[i], i );
			this.alphabet[i] = alphabet[i];
			if( longestCharacter < alphabet[i].length() ) {
				longestCharacter = alphabet[i].length();
			}
		}
	}

	/**
	 * Creates a new {@link DiscreteAlphabet} from a minimal and a maximal
	 * value, i.e. in <code>[min,max]</code>.
	 * 
	 * @param min
	 *            the minimal value (inclusive)
	 * @param max
	 *            the maximal value (inclusive)
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>min &gt; max</code>
	 */
	public DiscreteAlphabet( int min, int max ) throws IllegalArgumentException {
		if( min > max ) {
			throw new IllegalArgumentException( "The maximal value has to be equal or greater than the minimal value." );
		}
		String[] alphabet = new String[max - min + 1];
		for( int i = min; i <= max; i++ ) {
			alphabet[i - min] = "" + i;
		}

		try {
			init( alphabet, false );
		} catch ( Exception doesNotHappen ) {
			IllegalArgumentException i = new IllegalArgumentException( doesNotHappen.getCause() );
			i.setStackTrace( doesNotHappen.getStackTrace() );
			throw i;
		}
	}

	/**
	 * Creates a new {@link DiscreteAlphabet} from a given alphabet as a
	 * {@link String} array. This {@link String} array is cloned internally.
	 * 
	 * @param alphabet
	 *            the alphabet as {@link String} array
	 * @param caseInsensitive
	 *            indicates if there will be no difference between lowercase and
	 *            uppercase letters/symbols in the alphabet (no case sensitivity)
	 * 
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 * @throws IllegalArgumentException
	 *             if one of the symbols is either empty or a white-space
	 *             character
	 */
	public DiscreteAlphabet( boolean caseInsensitive, String... alphabet ) throws DoubleSymbolException, IllegalArgumentException {
		init( alphabet, caseInsensitive );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 200 );
		XMLParser.appendObjectWithTags( xml, alphabet, "symbols" );
		XMLParser.appendObjectWithTags( xml, caseInsensitive, "caseInsensitive" );
		addToXML(xml);
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}
	
	protected void addToXML(StringBuffer xml){
		
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo( Alphabet b ) {
		//if( !getClass().equals( b.getClass() ) ) { /TODO FIXME problem with DNAAlphabet vs {A,C,G,T}
//			return getClass().getName().compareTo( b.getClass().getName() );
//		} else {
			if( b instanceof ContinuousAlphabet ) {//TODO not perfect
				return -1;
			}
			
			if( b == this ) {
				return 0;
			}
			DiscreteAlphabet a = (DiscreteAlphabet)b;
			if( a.alphabet.length != alphabet.length ) {
				return a.alphabet.length - alphabet.length;
			}
			if( a.caseInsensitive != caseInsensitive ) {
				return a.caseInsensitive ? -1 : 1;
			}
			int i = 0;
			if( caseInsensitive ) {
				for( ; i < alphabet.length; i++ ) {
					if( !alphabet[i].equalsIgnoreCase( a.getSymbolAt( i ) ) ) {
						break;
					}
				}
			} else {
				for( ; i < alphabet.length; i++ ) {
					if( !alphabet[i].equals( a.getSymbolAt( i ) ) ) {
						break;
					}
				}
			}
			if( i < alphabet.length ) {
				return alphabet[i].compareTo( a.getSymbolAt( i ) );
			} else {
				return 0;
			}
//		}
	}

	/**
	 * Returns the code of a given symbol.
	 * 
	 * @param symbol
	 *            the given symbol
	 * 
	 * @return the code of a given symbol
	 * 
	 * @throws WrongAlphabetException
	 *             if the symbol is not defined in the alphabet
	 */
	public int getCode( String symbol ) throws WrongAlphabetException {
		if( caseInsensitive ) {
			symbol = symbol.toUpperCase();
		}
		Integer i = hash.get( symbol );
		if( i == null ) {
			throw new WrongAlphabetException( "Symbol \"" + symbol + "\" from input not defined in alphabet: " + Arrays.toString(alphabet) );
		}
		return i;
	}

	/**
	 * Returns the length of the longest &quot;symbol&quot; in the alphabet.
	 * 
	 * @return the length of the longest &quot;symbol&quot;
	 */
	public final int getMaximalSymbolLength() {
		return longestCharacter;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.Alphabet#getMin()
	 */
	@Override
	public double getMin() {
		return 0;
	}

	/**
	 * Returns the symbol at position <code>i</code> in the alphabet.
	 * 
	 * @param i
	 *            the position in the alphabet
	 * 
	 * @return the symbol at position <code>i</code>
	 */
	public final String getSymbolAt( int i ) {
		return alphabet[i];
	}

	/**
	 * Indicates if the alphabet ignores the case.
	 * 
	 * @return <code>true</code> if the alphabet ignores the case,
	 *         <code>false</code> otherwise
	 */
	public final boolean ignoresCase() {
		return caseInsensitive;
	}

	/**
	 * Indicates if <code>candidate</code> is an element of the internal
	 * interval.
	 * 
	 * @param candidate
	 *            the value to be tested
	 * 
	 * @return <code>true</code> if <code>candidate</code> is an element of the
	 *         internal interval, <code>false</code> otherwise
	 */
	public final boolean isEncodedSymbol( int candidate ) {
		return ( 0 <= candidate ) && ( candidate < alphabet.length );
	}

	/**
	 * Tests if a given symbol is contained in the alphabet.
	 * 
	 * @param candidat
	 *            the candidat symbol
	 * 
	 * @return <code>true</code> if the <code>candidat</code> is a symbol of the
	 *         alphabet, <code>false</code> otherwise
	 */
	public final boolean isSymbol( String candidat ) {
		return hash.containsKey( caseInsensitive ? candidat.toUpperCase() : candidat );
	}

	/**
	 * Returns the number of symbols in the calling alphabet.
	 * 
	 * @return the number of symbols in the alphabet
	 */
	@Override
	public final double length() {
		return alphabet.length;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.Alphabet#toString()
	 */
	@Override
	public String toString() {
		String erg = "{\"" + alphabet[0];
		for( int i = 1; i < alphabet.length; i++ ) {
			erg += "\", \"" + alphabet[i];
		}
		return erg + "\"}";
	}

	/**
	 * Class for the {@link de.jstacs.parameters.ParameterSet} of a
	 * {@link DiscreteAlphabet}.
	 * 
	 * @author Jan Grau
	 */
	public static class DiscreteAlphabetParameterSet extends AlphabetParameterSet {
		
		/**
		 * This constructor should only be used for parameter sets that are intended to created subclasses of {@link DiscreteAlphabet}. 
		 * 
		 * @param clazz the class the should be created with this parameter set
		 */
		protected DiscreteAlphabetParameterSet( Class<? extends DiscreteAlphabet> clazz ) throws ParameterException {
			super( clazz );
			parameters.add( new SimpleParameter( DataType.STRING,
					"Values of the alphabet",
					"The possible values of the discrete alphabet." + "If the alphabet consists of single characters, e.g. A, C, G, and T,"
							+ " the values may be set as a single string, e.g. \"ACGT\"."
							+ "If the alphabet consists of multi-character symbols, e.g. Gly, Asp, Ser,"
							+ "the symbols must be separated by spaces.",
					true ) );
			parameters.add( new SelectionParameter( DataType.BOOLEAN,
					new String[]{ "Case insensitive", "Case sensitive" },
					new Boolean[]{ true, false },
					"Case insensitive",
					"Use the alphabet case insensitive",
					true ) );
		}
		
		/**
		 * Creates a new {@link DiscreteAlphabetParameterSet} with empty values.
		 * @throws ParameterException if the parameters could not be created
		 * 
		 * @see de.jstacs.data.alphabets.Alphabet.AlphabetParameterSet#AlphabetParameterSet(Class)
		 */
		public DiscreteAlphabetParameterSet() throws ParameterException {
			this( DiscreteAlphabet.class );
			
		}

		/**
		 * Creates a new {@link DiscreteAlphabetParameterSet} from an alphabet
		 * given as a {@link String} array.
		 * 
		 * @param alphabet
		 *            the alphabet as {@link String} array
		 * @param caseInsensitive
		 *            indicates if the {@link DiscreteAlphabet} shall be case
		 *            insensitive
		 * 
		 * @throws Exception
		 *             if the parameter set could not be created
		 * 
		 * @see de.jstacs.data.alphabets.DiscreteAlphabet.DiscreteAlphabetParameterSet#DiscreteAlphabetParameterSet()
		 */
		public DiscreteAlphabetParameterSet( String[] alphabet, boolean caseInsensitive ) throws Exception {
			this();
			String alphString = "" + alphabet[0];
			for( int i = 1; i < alphabet.length; i++ ) {
				alphString += " " + alphabet[i];
			}
			parameters.get( 0 ).setValue( alphString );
			parameters.get( 1 ).setValue( new String( caseInsensitive ? "Case insensitive" : "Case sensitive" ) );
		}

		/**
		 * Creates a new {@link DiscreteAlphabetParameterSet} from an alphabet
		 * of symbols given as a <code>char</code> array.
		 * 
		 * @param alphabet
		 *            the array of symbols
		 * @param caseInsensitive
		 *            indicates if the {@link DiscreteAlphabet} shall be case
		 *            insensitive
		 * 
		 * @throws Exception
		 *             if the parameter set could not be created
		 * 
		 * @see  DiscreteAlphabet.DiscreteAlphabetParameterSet#DiscreteAlphabetParameterSet()
		 */
		public DiscreteAlphabetParameterSet( char[] alphabet, boolean caseInsensitive ) throws Exception {
			this();
			String alphString = "" + alphabet[0];
			for( int i = 1; i < alphabet.length; i++ ) {
				alphString += " " + alphabet[i];
			}
			parameters.get( 0 ).setValue( alphString );
			parameters.get( 1 ).setValue( new String( caseInsensitive ? "Case insensitive" : "Case sensitive" ) );
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link DiscreteAlphabetParameterSet} out of its XML
		 * representation.
		 * 
		 * @param representation
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link DiscreteAlphabetParameterSet} could not be
		 *             reconstructed out of the XML representation (the
		 *             {@link StringBuffer} <code>representation</code> could
		 *             not be parsed)
		 * 
		 * @see de.jstacs.data.alphabets.Alphabet.AlphabetParameterSet#AlphabetParameterSet(StringBuffer)
		 * @see de.jstacs.Storable
		 */
		public DiscreteAlphabetParameterSet( StringBuffer representation ) throws NonParsableException {
			super( representation );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "An alphabet that consists of discrete values.";
		}
	}
}