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

import de.jstacs.Singleton;

/**
 * This class implements the discrete alphabet that is used for proteins (one letter code).
 * 
 * @author Jens Keilwagen
 */
public final class ProteinAlphabet extends DiscreteAlphabet implements Singleton {

	public static final ProteinAlphabet SINGLETON = get();
	
	private static ProteinAlphabet get() {
		ProteinAlphabet res = null;
		try {
			res = new ProteinAlphabet();
		} catch (Exception doesNotHappen) {
			throw new RuntimeException( doesNotHappen.getMessage() );
		}
		return res;
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
	private ProteinAlphabet() throws DoubleSymbolException, IllegalArgumentException {
		super( true, ProteinAlphabetParameterSet.AMINOACID );
		this.parameters = ProteinAlphabetParameterSet.SINGLETON;
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
	public static final class ProteinAlphabetParameterSet extends AlphabetParameterSet implements Singleton {

		public static final ProteinAlphabetParameterSet SINGLETON = get();
		
		private static final String[] AMINOACID = { "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V" };

		private static ProteinAlphabetParameterSet get() {
			ProteinAlphabetParameterSet res = null;
			try {
				res = new ProteinAlphabetParameterSet();
			} catch (Exception doesNotHappen) {
				throw new RuntimeException( doesNotHappen.getMessage() );
			}
			return res;
		}
		
		/**
		 * Creates a new {@link ProteinAlphabetParameterSet}.
		 * 
		 * @throws Exception
		 *             if an error occurred
		 * 
		 * @see de.jstacs.data.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(Class) Alphabet.AlphabetParameterSet#AlphabetParameterSet(Class)
		 */
		private ProteinAlphabetParameterSet() throws Exception {
			super( ProteinAlphabet.class );
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