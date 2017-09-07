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
public final class IUPACDNAAlphabet extends ComplementableDiscreteAlphabet implements Singleton {

	/**
	 * The only instance of this class.
	 * 
	 * @see Singleton
	 */
	public final static IUPACDNAAlphabet SINGLETON = get(); 
	
	private static IUPACDNAAlphabet get() {
		IUPACDNAAlphabet res = null;
		try {
			res = new IUPACDNAAlphabet();
		} catch (Exception doesNotHappen) {
			throw new RuntimeException( doesNotHappen.getMessage() );
		}
		return res;
	}

	/**
	 * The main constructor. Creates a new {@link IUPACDNAAlphabet} with the standard
	 * UIPAC DNA-alphabet.
	 * 
	 * @throws DoubleSymbolException
	 *             if one of the symbols occurred more than once
	 * @throws IllegalArgumentException
	 *             if space or tab will be used as symbols
	 * 
	 * @see ComplementableDiscreteAlphabet#ComplementableDiscreteAlphabet(boolean, String...)
	 */
	private IUPACDNAAlphabet() throws DoubleSymbolException, IllegalArgumentException {
		super( true, IUPACDNAAlphabetParameterSet.DNA );
		this.parameters = IUPACDNAAlphabetParameterSet.SINGLETON;
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
	 * Indicates if <code>query</code> is contained in <code>code</code>
	 * according to the IUPAC DNA alphabet.
	 * @param query the query symbol
	 * @param code the reference symbol
	 * @return if <code>query</code> is contained in <code>code</code>
	 * @throws WrongAlphabetException if either of the arguments is not defined in the IUPAC alphabet
	 */
	public boolean isPart( String query, String code ) throws WrongAlphabetException {
		return isPart( getCode(query), getCode(code) );
	}
	
	/**
	 * Indicates if <code>query</code> is contained in <code>code</code>
	 * according to the IUPAC DNA alphabet.
	 * @param query the code of the query symbol
	 * @param code the code of the reference symbol
	 * @return if <code>query</code> is contained in <code>code</code>
	 */
	public boolean isPart( int query, int code ) {
		return IUPACDNAAlphabetParameterSet.subset[query][code];
	}
	
	/**
	 * The parameter set for a {@link IUPACDNAAlphabet}.
	 * 
	 * @author Jan Grau, Jens Keilwagen
	 */
	public static final class IUPACDNAAlphabetParameterSet extends AlphabetParameterSet<IUPACDNAAlphabet> implements Singleton {

		/**
		 * The only instance of this class.
		 * 
		 * @see Singleton
		 */
		public final static IUPACDNAAlphabetParameterSet SINGLETON = get(); 
		
		private static IUPACDNAAlphabetParameterSet get() {
			IUPACDNAAlphabetParameterSet res = null;
			try {
				res = new IUPACDNAAlphabetParameterSet();
			} catch (Exception doesNotHappen) {
				throw new RuntimeException( doesNotHappen.getMessage() );
			}
			return res;
		}
		
		private static final String[] DNA = { "A", "C", "R", "K", "B", "D", "H", "V", "M", "Y", "G", "T", "S", "W", "N"};
		private static final boolean[][] subset = new boolean[DNA.length][DNA.length];
		static {
			for( int i = 0; i < DNA.length; i++ ) {
				subset[i][i] = true;
			}
			
			try {
				set( 2, 0, 10 );
				set( 3, 10, 11 );
				set( 4, 1, 10, 11,     3, 9, 12 );
				set( 5, 0, 10, 11,     3, 13, 2 );
				set( 6, 0, 1, 11,      9, 13, 8 );
				set( 7, 0, 1, 10,      12, 2, 8 );
				set( 8, 0, 1 );
				set( 9, 1, 11 );
				set( 12, 1, 10 );
				set( 13, 0, 11 );
				set( 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 );
			} catch( WrongAlphabetException doesNotHappen ) {
				throw new RuntimeException(doesNotHappen.getCause());
			}
		}
		
		private static void set( int ambigious, int... subsets ) throws WrongAlphabetException {
			for( int i = 0; i < subsets.length; i++ ) {
				subset[subsets[i]][ambigious] = true;
			}
		}

		/**
		 * Creates a new {@link UIPACDNAAlphabetParameterSet}.
		 * 
		 * @throws Exception
		 *             if an error occurred
		 * 
		 * @see de.jstacs.data.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(Class) Alphabet.AlphabetParameterSet#AlphabetParameterSet(Class)
		 */
		private IUPACDNAAlphabetParameterSet() throws Exception {
			super( IUPACDNAAlphabet.class );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "A IUPAC alphabet for DNA.";
		}
	}
}
