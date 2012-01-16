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
 * This class implements the discrete alphabet that is used for DNA.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public final class DNAAlphabet extends ComplementableDiscreteAlphabet implements Singleton {

	/**
	 * The only instance of this class.
	 * 
	 * @see Singleton
	 */
	public final static DNAAlphabet SINGLETON = get(); 
	
	private static DNAAlphabet get() {
		DNAAlphabet res = null;
		try {
			res = new DNAAlphabet();
		} catch (Exception doesNotHappen) {
			throw new RuntimeException( doesNotHappen.getMessage() );
		}
		return res;
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
	private DNAAlphabet() throws DoubleSymbolException, IllegalArgumentException {
		super( true, DNAAlphabetParameterSet.DNA );
		this.parameters = DNAAlphabetParameterSet.SINGLETON;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.ComplementableDiscreteAlphabet#getComplementaryCode(int)
	 */
	@Override
	public int getComplementaryCode( int code ) {
		return 3 - code;
	}

	/**
	 * The parameter set for a {@link DNAAlphabet}.
	 * 
	 * @author Jan Grau, Jens Keilwagen
	 */
	public static final class DNAAlphabetParameterSet extends AlphabetParameterSet<DNAAlphabet> implements Singleton {

		/**
		 * The only instance of this class.
		 * 
		 * @see Singleton
		 */
		public final static DNAAlphabetParameterSet SINGLETON = get(); 
		
		private static DNAAlphabetParameterSet get() {
			DNAAlphabetParameterSet res = null;
			try {
				res = new DNAAlphabetParameterSet();
			} catch (Exception doesNotHappen) {
				throw new RuntimeException( doesNotHappen.getMessage() );
			}
			return res;
		}
		
		private static final String[] DNA = { "A", "C", "G", "T" };

		/**
		 * Creates a new {@link DNAAlphabetParameterSet}.
		 * 
		 * @throws Exception
		 *             if an error occurred
		 * 
		 * @see de.jstacs.data.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(Class) Alphabet.AlphabetParameterSet#AlphabetParameterSet(Class)
		 */
		private DNAAlphabetParameterSet() throws Exception {
			super( DNAAlphabet.class );
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
