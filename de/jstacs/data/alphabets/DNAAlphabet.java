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

import de.jstacs.NonParsableException;

/**
 * This class implements the discrete alphabet that is used for DNA.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public final class DNAAlphabet extends ComplementableDiscreteAlphabet {

	public static DNAAlphabet get() {
		DNAAlphabet res = null;
		try {
			res = new DNAAlphabet();
		} catch (Exception doesNotHappen) {
			System.err.println( "Unexpected Error:" );
			doesNotHappen.printStackTrace();
			System.exit(1);
		}
		return res;
	}
	
	
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link DNAAlphabet} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DNAAlphabet} could not be reconstructed out of
	 *             the XML representation (the {@link StringBuffer}
	 *             <code>representation</code> could not be parsed)
	 * 
	 * @see ComplementableDiscreteAlphabet#ComplementableDiscreteAlphabet(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public DNAAlphabet( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * The constructor for the {@link de.jstacs.InstantiableFromParameterSet}
	 * interface. Creates a new {@link DNAAlphabet} from a given parameter set.
	 * 
	 * @param parameters
	 *            the given set of parameters
	 * 
	 * @throws IllegalArgumentException
	 *             if space or tab will be used as symbols
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 * 
	 * @see DNAAlphabet#DNAAlphabet()
	 */

	public DNAAlphabet( DNAAlphabetParameterSet parameters ) throws IllegalArgumentException, DoubleSymbolException {
		this();
		this.parameters = parameters;
	}

	/**
	 * The main constructor. Creates a new {@link DNAAlphabet} with the standard
	 * DNA-alphabet (&quot;A&quot;, &quot;C&quot;, &quot;G&quot;, &quot;T&quot;).
	 * 
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 * @throws IllegalArgumentException
	 *             if space or tab will be used as symbols
	 * 
	 * @see ComplementableDiscreteAlphabet#ComplementableDiscreteAlphabet(boolean, String...)
	 */
	public DNAAlphabet() throws DoubleSymbolException, IllegalArgumentException {
		super( true, DNAAlphabetParameterSet.DNA );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.ComplementableDiscreteAlphabet#getComplementaryCode(int)
	 */
	@Override
	public int getComplementaryCode( int code ) {
		return 3 - code;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.ComplementableDiscreteAlphabet#getCurrentParameterSet()
	 */
	@Override
	public AlphabetParameterSet getCurrentParameterSet() throws Exception {
		if( this.parameters == null ) {
			return new DNAAlphabetParameterSet();
		} else {
			return this.parameters;
		}
	}

	/**
	 * The parameter set for a {@link DNAAlphabet}.
	 * 
	 * @author Jan Grau, Jens Keilwagen
	 */
	public static final class DNAAlphabetParameterSet extends AlphabetParameterSet {

		private static final String[] DNA = { "A", "C", "G", "T" };

		/**
		 * Creates a new {@link DNAAlphabetParameterSet}.
		 * 
		 * @throws Exception
		 *             if an error occurred
		 * 
		 * @see de.jstacs.data.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(Class) Alphabet.AlphabetParameterSet#AlphabetParameterSet(Class)
		 */
		public DNAAlphabetParameterSet() throws Exception {
			super( DNAAlphabet.class );
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link DNAAlphabetParameterSet} out of its XML
		 * representation.
		 * 
		 * @param representation
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link DNAAlphabetParameterSet} could not be
		 *             reconstructed out of the XML representation (the
		 *             {@link StringBuffer} <code>representation</code> could
		 *             not be parsed)
		 * 
		 * @see de.jstacs.data.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(StringBuffer) Alphabet.AlphabetParameterSet#AlphabetParameterSet(StringBuffer)
		 * @see de.jstacs.Storable
		 */
		public DNAAlphabetParameterSet( StringBuffer representation ) throws NonParsableException {
			super( representation );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "An alphabet for DNA.";
		}
	}
}
