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
package supplementary.cookbook.recipes;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.GenericComplementableDiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;

/**
 * This class exemplarily shows how to create a user-specified alphabet.
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