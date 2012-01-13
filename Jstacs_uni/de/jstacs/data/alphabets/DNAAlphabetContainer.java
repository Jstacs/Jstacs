package de.jstacs.data.alphabets;

import de.jstacs.Singleton;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;
import de.jstacs.parameters.InstanceParameterSet;

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
	
	public final static class DNAAlphabetContainerParameterSet extends InstanceParameterSet<DNAAlphabetContainer> implements Singleton {

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
			// TODO Auto-generated method stub
			return null;
		}
	}
}
