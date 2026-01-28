package de.jstacs.sequenceScores;

/**
 * Interface for {@link SequenceScore} that provide additional methods for computing scores of infix sequences
 * and filtering infix sequences.
 * 
 * @author Jan Grau
 *
 */
public interface QuickScanningSequenceScore extends SequenceScore {

	/**
	 * Computes the position-wise scores of an infix of the sequence <code>seq</code> (which must be encoded by the
	 * same alphabet as this {@link QuickScanningSequenceScore}) beginning at <code>start</code> and extending for <code>length</code>
	 * positions. The scores are computed per position and filled into the provided array <code>scores</code>, which is of the same
	 * length as <code>seq</code>. 
	 * 
	 * This must be implemented such that {@link de.jstacs.utils.ToolBox#sum(double...)} applied to <code>scores</code>
	 * computed from <code>0</code> to <code>seq.length</code> returns the same value as {@link #getLogScoreFor(de.jstacs.data.sequences.Sequence)} called
	 * on the {@link de.jstacs.data.sequences.IntSequence} created from <code>seq</code>.
	 * 
	 * @param seq the sequence
	 * @param start the start of the infix
	 * @param length the length of the infix
	 * @param scores the array of scores to be (partly) filled
	 */
	public void fillInfixScore(int[] seq, int start, int length, double[] scores);
	
	
	/**
	 * Computes arrays that indicate, for a given set of starting positions and a given k-mer length, if a sequence
	 * containing this k-mer may yield a score above <code>threshold</code>, choosing the best-scoring option among
	 * all non-specified positions (i.e., those outside the k-mer).
	 * This method is implemented as an upper bound on the scores, i.e., there may be k-mers that are considered
	 * to score above threshold (i.e., that have entry <code>true</code>) although they are not, but there may not be k-mers that are considered to be below
	 * threshold (i.e., that have entry <code>false</code>), although there exist sequences containing this k-mer that do.
	 * The returned array is indexed by the starting positions (in the same order as provided in <code>starts</code>) in the first dimension, and in the second dimension
	 * it is indexed by an integer representation of the k-mers, assigning the highest priority to the first k-mer position, i.e.,
	 * 
	 * $$\sum_{\ell=0}^{L-1} A^{L-1-\ell} \cdot e(x_\ell) $$
	 * where \( A \) denotes the size of the alphabet, \( L \) the length of the k-mer (starting at 0 in this case), 
	 * and \( e(\ldots) \) denotes the function encoding symbols from the alphabet as integers (see {@link de.jstacs.data.alphabets.DiscreteAlphabet}).
	 * 
	 * @param kmer the k-mer length
	 * @param thresh the threshold
	 * @param start the starting position(s)
	 * @return the indicator arrays
	 */
	public boolean[][] getInfixFilter( int kmer, double thresh, int... start );
	
}
