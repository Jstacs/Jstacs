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


public class DeBruijnSequenceGenerator {
	
	public static void main(String[] args) throws Exception {
		
		int k=6;
		int m=3;
		
		CyclicSequenceAdaptor seq = generate( DNAAlphabet.SINGLETON, k, 0 );
		System.out.println(seq);
		HashMap<Sequence, Integer> map = new HashMap<Sequence, Integer>();
		DataSet ds = new DataSet( new DataSet( "", seq.getSuperSequence( seq.getLength()+k-1 ) ), k );
		for(int i=0;i<ds.getNumberOfElements();i++){
			map.put( ds.getElementAt( i ), 1 );
		}
		
		DiscreteSequenceEnumerator en = new DiscreteSequenceEnumerator( seq.getAlphabetContainer(), k, false );		
		
		while(en.hasMoreElements()){
			Sequence temp = en.nextElement();
			if(map.get( temp ) == null){
				throw new Exception( temp.toString() );
			}
		}
		System.out.println("all");
		
	}
	
	
	
	public static CyclicSequenceAdaptor[] generate(DiscreteAlphabet alphabet, int n) throws WrongAlphabetException, WrongSequenceTypeException, OperationNotSupportedException{
		/*CyclicSequenceAdaptor[] seqs = new CyclicSequenceAdaptor[(int)alphabet.length()];
		for(int i=0;i<seqs.length;i++){
			seqs[i] = generate( DNAAlphabet.SINGLETON, n, i );
		}
		
		return seqs;*/
		return new CyclicSequenceAdaptor[]{generate(alphabet, n, 0)};
	}
	
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
