package de.jstacs.clustering.distances;

import java.util.Random;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.ByteSequence;
import de.jstacs.data.sequences.CyclicSequenceAdaptor;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.sequenceScores.statisticalModels.StatisticalModel;

/**
 * Class for a distance metric between {@link StatisticalModel}s based on the correlation of score
 * profiles on random sequences.
 * 
 * For two {@link StatisticalModel}s {@latex.inline $M_1$} and {@latex.inline $M_2$}, we compute the score profiles
 * {@latex.inline $s_1(x,M_1)$} and {@latex.inline $s_2(x,M_2)$} on a random sequence {@latex.inline x} of length {@latex.inline $4^n$}. 
 * The distance is then defined based on the Pearson correlation as {@latex.inline $1 - cor( s_1(x,M_1), s_2(x,M_2) )$} between these score profiles,
 * maximizing over suitable shifts of the score profiles and both strand orientations. 
 * 
 * @author Jan Grau
 *
 */
public class RandomSequenceScoreDistance extends SequenceScoreDistance {

	private static Random r = new Random(117);

	/**
	 * Creates a distance using a random sequence of length {@latex.inline $|A|^n$}.
	 * @param alphabet the alphabet
	 * @param n the exponent of the sequence length
	 * @param exp if exponential profiles should be computed
	 * @throws WrongAlphabetException if the alphabet does not match
	 * @throws WrongSequenceTypeException if the sequence could not be created for this alphabet
	 */
	public RandomSequenceScoreDistance(DiscreteAlphabet alphabet, int n, boolean exp) throws WrongAlphabetException, WrongSequenceTypeException {
		super(createSequences(alphabet,n), exp);
	}
	
	@Override
	public double[][] getProfile(StatisticalModel o, boolean rc) throws Exception{
		return DeBruijnMotifComparison.getProfilesForMotif( seqs, o, rc, exp );
	}

	/**
	 * Creates a new random sequence of the given length and alphabet.
	 * @param alphabet the alphabet
	 * @param n the exponent of the sequence length (see above)
	 * @return the sequence
	 * @throws WrongAlphabetException if the alphabet does not match
	 * @throws WrongSequenceTypeException if the sequence could not be created for this alphabet
	 */
	public static CyclicSequenceAdaptor[] createSequences(
			DiscreteAlphabet alphabet, int n) throws WrongAlphabetException, WrongSequenceTypeException {
		int l = (int)Math.pow(alphabet.length(), n);
		int al = (int) alphabet.length();
		byte[] seq = new byte[l];
		for(int i=0;i<l;i++){
			seq[i] = (byte) r.nextInt(al);
		}
		
		AlphabetContainer ac = null;
		if(alphabet instanceof DNAAlphabet){
			ac = DNAAlphabetContainer.SINGLETON;
		}else{
			ac = new AlphabetContainer(alphabet);
		}
		
		return new CyclicSequenceAdaptor[]{
				new CyclicSequenceAdaptor(new ByteSequence(ac, seq))
		};
	}

}
