package supplementary.cookbook.recipes;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.GenericComplementableDiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;

/**
 * This class exemplarily shows how two create a user-specified alphabet.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class AlphabetCreation {

	public static void main(String[] args) throws Exception {
		String[] symbols = {"A", "C", "G", "T", "-"};
		//new alphabet
		DiscreteAlphabet abc = new DiscreteAlphabet( true, symbols );
		
		//new alphabet that allows to build the reverse complement of a sequence
		int[] revComp = new int[symbols.length];
		revComp[0] = 3; //symbols[0]^rc = symbols[3]
		revComp[1] = 2; //symbols[1]^rc = symbols[2]
		revComp[2] = 1; //symbols[2]^rc = symbols[1]
		revComp[3] = 0; //symbols[3]^rc = symbols[0]
		revComp[4] = 4; //symbols[4]^rc = symbols[4]
		GenericComplementableDiscreteAlphabet abc2 = new GenericComplementableDiscreteAlphabet( true, symbols, revComp );
		
		Sequence seq = Sequence.create( new AlphabetContainer( abc2 ), "ACGT-" );
		Sequence rc = seq.reverseComplement();
		System.out.println( seq );
		System.out.println( rc );
	}
}