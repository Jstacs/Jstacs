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
import de.jstacs.data.WrongAlphabetException;

/**
 * This class implements a discrete alphabet for the IUPAC DNA code.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public final class UIPACDNAAlphabet extends ComplementableDiscreteAlphabet implements Singleton {

	/**
	 * The only instance of this class.
	 * 
	 * @see Singleton
	 */
	public final static UIPACDNAAlphabet SINGLETON = get(); 
	
	private static UIPACDNAAlphabet get() {
		UIPACDNAAlphabet res = null;
		try {
			res = new UIPACDNAAlphabet();
		} catch (Exception doesNotHappen) {
			throw new RuntimeException( doesNotHappen.getMessage() );
		}
		return res;
	}

	/**
	 * The main constructor. Creates a new {@link UIPACDNAAlphabet} with the standard
	 * UIPAC DNA-alphabet.
	 * 
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 * @throws IllegalArgumentException
	 *             if space or tab will be used as symbols
	 * 
	 * @see ComplementableDiscreteAlphabet#ComplementableDiscreteAlphabet(boolean, String...)
	 */
	private UIPACDNAAlphabet() throws DoubleSymbolException, IllegalArgumentException {
		super( true, UIPACDNAAlphabetParameterSet.DNA );
		this.parameters = UIPACDNAAlphabetParameterSet.SINGLETON;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.ComplementableDiscreteAlphabet#getComplementaryCode(int)
	 */
	@Override
	public int getComplementaryCode( int code ) {
		if(code < 12){
			return 11-code;
		}else{
			return code;
		}
	}
	
	/**
	 * The parameter set for a {@link UIPACDNAAlphabet}.
	 * 
	 * @author Jan Grau, Jens Keilwagen
	 */
	public static final class UIPACDNAAlphabetParameterSet extends AlphabetParameterSet<UIPACDNAAlphabet> implements Singleton {

		/**
		 * The only instance of this class.
		 * 
		 * @see Singleton
		 */
		public final static UIPACDNAAlphabetParameterSet SINGLETON = get(); 
		
		private static UIPACDNAAlphabetParameterSet get() {
			UIPACDNAAlphabetParameterSet res = null;
			try {
				res = new UIPACDNAAlphabetParameterSet();
			} catch (Exception doesNotHappen) {
				throw new RuntimeException( doesNotHappen.getMessage() );
			}
			return res;
		}
		
		private static final String[] DNA = { "A", "C", "R", "K", "B", "D", "H", "V", "M", "Y", "G", "T", "S", "W", "N"};
		

		/**
		 * Creates a new {@link UIPACDNAAlphabetParameterSet}.
		 * 
		 * @throws Exception
		 *             if an error occurred
		 * 
		 * @see de.jstacs.data.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(Class) Alphabet.AlphabetParameterSet#AlphabetParameterSet(Class)
		 */
		private UIPACDNAAlphabetParameterSet() throws Exception {
			super( UIPACDNAAlphabet.class );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "A UIPAC alphabet for DNA.";
		}
	}
}
