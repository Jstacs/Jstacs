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
 * This class implements the discrete alphabet that is used for proteins (one letter code).
 * 
 * @author Jens Keilwagen
 */
public final class ProteinAlphabet extends DiscreteAlphabet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ProteinAlphabet} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ProteinAlphabet} could not be reconstructed out of
	 *             the XML representation (the {@link StringBuffer}
	 *             <code>representation</code> could not be parsed)
	 * 
	 * @see DiscreteAlphabet#DiscreteAlphabet(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public ProteinAlphabet( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * The constructor for the {@link de.jstacs.InstantiableFromParameterSet}
	 * interface. Creates a new {@link ProteinAlphabet} from a given parameter set.
	 * 
	 * @param parameters
	 *            the given set of parameters
	 * 
	 * @throws IllegalArgumentException
	 *             if space or tab will be used as symbols
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 * 
	 * @see ProteinAlphabet#ProteinAlphabet()
	 */

	public ProteinAlphabet( ProteinAlphabetParameterSet parameters ) throws IllegalArgumentException, DoubleSymbolException {
		this();
		this.parameters = parameters;
	}

	/**
	 * The main constructor. Creates a new {@link ProteinAlphabet} with the one letter code
	 * Protein-alphabet.
	 * 
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 * @throws IllegalArgumentException
	 *             if space or tab will be used as symbols
	 * 
	 * @see DiscreteAlphabet#DiscreteAlphabet(boolean, String...)
	 */
	public ProteinAlphabet() throws DoubleSymbolException, IllegalArgumentException {
		super( true, ProteinAlphabetParameterSet.AMINOACID );
	}


	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.ComplementableDiscreteAlphabet#getCurrentParameterSet()
	 */
	@Override
	public AlphabetParameterSet getCurrentParameterSet() throws Exception {
		if( this.parameters == null ) {
			return new ProteinAlphabetParameterSet();
		} else {
			return this.parameters;
		}
	}

	/**
	 * The parameter set for a {@link ProteinAlphabet}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static final class ProteinAlphabetParameterSet extends AlphabetParameterSet {

		private static final String[] AMINOACID = { "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V" };

		/**
		 * Creates a new {@link ProteinAlphabetParameterSet}.
		 * 
		 * @throws Exception
		 *             if an error occurred
		 * 
		 * @see de.jstacs.data.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(Class) Alphabet.AlphabetParameterSet#AlphabetParameterSet(Class)
		 */
		public ProteinAlphabetParameterSet() throws Exception {
			super( ProteinAlphabet.class );
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link ProteinAlphabetParameterSet} out of its XML
		 * representation.
		 * 
		 * @param representation
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link ProteinAlphabetParameterSet} could not be
		 *             reconstructed out of the XML representation (the
		 *             {@link StringBuffer} <code>representation</code> could
		 *             not be parsed)
		 * 
		 * @see de.jstacs.data.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(StringBuffer) Alphabet.AlphabetParameterSet#AlphabetParameterSet(StringBuffer)
		 * @see de.jstacs.Storable
		 */
		public ProteinAlphabetParameterSet( StringBuffer representation ) throws NonParsableException {
			super( representation );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "The (one letter code) alphabet for protein.";
		}
	}
}
