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