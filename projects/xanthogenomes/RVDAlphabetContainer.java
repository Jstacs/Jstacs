package projects.xanthogenomes;

import de.jstacs.Singleton;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.DoubleSymbolException;
import projects.xanthogenomes.Tools.Translator;

public final class RVDAlphabetContainer extends AlphabetContainer implements Singleton {
	
	/**
	 * The only instance of this class.
	 * 
	 * @see Singleton
	 */
	public static final RVDAlphabetContainer SINGLETON = get();

	private RVDAlphabetContainer() throws IllegalArgumentException, DoubleSymbolException{
		super(new DiscreteAlphabet( true, getRVDs() ));
	}
	
	private static String[] getRVDs(){
		AlphabetContainer prot = Translator.DEFAULT.getProteinAlphabet();
		String[] rvds = new String[(int)(prot.getAlphabetLengthAt( 0 )*prot.getAlphabetLengthAt( 0 ))];
		int k=0;
		for(int i=0;i<prot.getAlphabetLengthAt( 0 );i++){
			for(int j=0;j<prot.getAlphabetLengthAt( 0 );j++,k++){
				rvds[k] = prot.getSymbol( 0, i )+prot.getSymbol( 0, j );
			}
		}
		return rvds;
	}
	
	private static RVDAlphabetContainer get(){
		
		try {
			RVDAlphabetContainer RVDAlphabet = new RVDAlphabetContainer();
			return RVDAlphabet;
		} catch ( Exception e ) { 
			return null;
		}
	}
	
}