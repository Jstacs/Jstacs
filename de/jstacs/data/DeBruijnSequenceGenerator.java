package de.jstacs.data;

import java.util.HashMap;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.alphabets.ComplementableDiscreteAlphabet;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.CyclicSequenceAdaptor;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModelFactory;
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
	 * @throws WrongAlphabetException if the alphabet is 
	 * @throws WrongSequenceTypeException
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
	 * @throws WrongAlphabetException 
	 * @throws WrongSequenceTypeException 
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