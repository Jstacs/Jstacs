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
 * This abstract class indicates that an alphabet can be used to compute the
 * complement. The most important method is {@link #getComplementaryCode(int)}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public abstract class ComplementableDiscreteAlphabet extends DiscreteAlphabet {

	/**
	 * The constructor for the {@link de.jstacs.InstantiableFromParameterSet}
	 * interface. Creates a new {@link ComplementableDiscreteAlphabet} from a given set of
	 * parameters.
	 * 
	 * @param parameters
	 *            the parameter set for the {@link ComplementableDiscreteAlphabet}
	 * 
	 * @throws IllegalArgumentException
	 *             if space or tab will be used as symbols
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 * 
	 * @see de.jstacs.InstantiableFromParameterSet
	 */
	protected ComplementableDiscreteAlphabet( DiscreteAlphabetParameterSet parameters ) throws IllegalArgumentException, DoubleSymbolException {
		super( parameters );
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ComplementableDiscreteAlphabet} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ComplementableDiscreteAlphabet} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see DiscreteAlphabet#DiscreteAlphabet(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	protected ComplementableDiscreteAlphabet( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * Creates a new {@link ComplementableDiscreteAlphabet} from a given array
	 * of symbols.
	 * 
	 * @param alphabet
	 *            the array of symbols
	 * @param caseInsensitive
	 *            indicates if there will be no difference between lowercase and
	 *            uppercase letters/symbols in the alphabet (no case sensitivity)
	 * 
	 * @throws DoubleSymbolException
	 *             if a symbol occurs more than once in the alphabet
	 * @throws IllegalArgumentException
	 *             if one of the symbols is either empty or a white-space
	 *             character
	 * 
	 * @see DiscreteAlphabet#DiscreteAlphabet(boolean, String...)
	 */
	protected ComplementableDiscreteAlphabet( boolean caseInsensitive, String... alphabet ) throws DoubleSymbolException,
																							IllegalArgumentException {
		super( caseInsensitive, alphabet );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.DiscreteAlphabet#getCurrentParameterSet()
	 */
	@Override
	public abstract AlphabetParameterSet getCurrentParameterSet() throws Exception;

	/**
	 * This method returns the code of the symbol that is the complement of the
	 * symbol encoded by <code>code</code>.
	 * 
	 * @param code
	 *            the encoded symbol
	 * 
	 * @return the code of the complement
	 */
	public abstract int getComplementaryCode( int code );
}
