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
import de.jstacs.data.AlphabetContainer;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;

/**
 * This class implements a singleton for an {@link AlphabetContainer} that can be used for DNA. 
 * 
 * @see Singleton
 * @see DNAAlphabet
 * @see de.jstacs.data.DNADataSet
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public final class DNAAlphabetContainer extends AlphabetContainer implements Singleton {
	
	/**
	 * The only instance of this class.
	 * 
	 * @see Singleton
	 */
	public static final DNAAlphabetContainer SINGLETON = get();
	
	private static DNAAlphabetContainer get() {
		DNAAlphabetContainer res = null;
		try {
			res = new DNAAlphabetContainer();
		} catch (Exception doesNotHappen) {
			throw new RuntimeException( doesNotHappen.getMessage() );
		}
		return res;
	}
	
	private DNAAlphabetContainer() throws IllegalArgumentException, NotInstantiableException, DoubleSymbolException {
		super( DNAAlphabet.SINGLETON );
		parameters = DNAAlphabetContainerParameterSet.SINGLETON;
	}

	/**
	 * This class implements a singleton for the {@link de.jstacs.parameters.ParameterSet} of a {@link DNAAlphabetContainer}. 
	 * 
	 * @author Jens Keilwagen
	 */
	public final static class DNAAlphabetContainerParameterSet extends AbstractAlphabetContainerParameterSet<DNAAlphabetContainer> implements Singleton {

		/**
		 * The only instance of this class.
		 * 
		 * @see Singleton
		 */
		public static final DNAAlphabetContainerParameterSet SINGLETON = get();
		
		private static DNAAlphabetContainerParameterSet get() {
			DNAAlphabetContainerParameterSet res = null;
			try {
				res = new DNAAlphabetContainerParameterSet();
			} catch (Exception doesNotHappen) {
				throw new RuntimeException( doesNotHappen.getMessage() );
			}
			return res;
		}
		
		private DNAAlphabetContainerParameterSet() throws Exception {
			super( DNAAlphabetContainer.class );
		}
		
		public DNAAlphabetContainerParameterSet clone() {
			return this;
		}
		
		@Override
		public String getInstanceComment() {
			return "An alphabet container for DNA.";
		}

		@Override
		public String getInstanceName() {
			return "DNA alphabet container";
		}

		@Override
		public int getPossibleLength() {
			return 0;
		}

		@Override
		public boolean isDiscrete() {
			return true;
		}

		@Override
		public boolean isSimple() {
			return true;
		}
	}
}
