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

package de.jstacs.data;

import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.CyclicSequenceAdaptor;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.utils.IntList;

/**
 * Generates De Buijn sequences using the algorithm from Frank Ruskey's Combinatorial Generation.
 * @author Jan Grau
 *
 */
public class DeBruijnSequenceGenerator {
	
	/**
	 * Generates a De Bruijn sequence of length {@latex.inline $|A|^n$}, where A denotes the alphabet.
	 * @param alphabet the alphabet
	 * @param n the exponent of length length, corresponds to the length of n-mers covered exactly once
	 * @return the sequence (wrapped in an array)
	 * @throws WrongAlphabetException forwarded from {@link IntSequence}, should not happen
	 * @throws WrongSequenceTypeException forwarded from {@link IntSequence}, should not happen
	 */
	public static CyclicSequenceAdaptor[] generate(DiscreteAlphabet alphabet, int n) throws WrongAlphabetException, WrongSequenceTypeException {
		
		return new CyclicSequenceAdaptor[]{generate(alphabet, n, 0)};
	}
	
	/**
	 * Generates a De Bruijn sequence using the supplied alphabet and the given alphabet shift, i.e., for a cyclic shift of the symbols 
	 * of the alphabet.
	 * @param alphabet the alphabet
	 * @param n the length of the covered n-mers
	 * @param alphabetShift the alphabet shift (0 equals no shift)
	 * @return the De Bruijn sequence
	 * @throws WrongAlphabetException forwarded from {@link IntSequence}, should not happen
	 * @throws WrongSequenceTypeException forwarded from {@link IntSequence}, should not happen
	 */
	public static CyclicSequenceAdaptor generate(DiscreteAlphabet alphabet, int n, int alphabetShift) throws WrongAlphabetException, WrongSequenceTypeException{
		int k = (int)alphabet.length();
		if(alphabetShift >= k || alphabetShift < 0){
			throw new IllegalArgumentException( "alphabet shift greater than alphabet size" );
		}
		
		IntList seq = new IntList();
		int[] a = new int[k*n];
		db(1,1,n,k,seq,a);
		int[] seqa = seq.toArray();
		
		if(alphabetShift != 0){
			for(int i=0;i<seqa.length;i++){
				seqa[i] = (seqa[i] + alphabetShift) % k;
			}
		}
		
		return new CyclicSequenceAdaptor( new IntSequence( new AlphabetContainer( alphabet ), seqa ),seqa.length );
	}
	
	private static void db(int t, int p, int n, int k, IntList seq, int[] a){
		if(t > n){
			if(n % p == 0){
				for(int j=1;j<p+1;j++){
					seq.add( a[j] );
				}
			}
		}else{
			a[t] = a[t-p];
			db(t+1,p,n,k,seq,a);
			for(int j=a[t-p]+1;j<k;j++){
				a[t] = j;
				db(t+1,t,n,k,seq,a);
			}
		}
	}

}
