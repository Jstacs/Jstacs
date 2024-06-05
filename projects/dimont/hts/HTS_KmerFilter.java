package projects.dimont.hts;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map.Entry;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.utils.ComparableElement;

public enum HTS_KmerFilter {
	
	HUDDINGE_DISTANCE,
	MINIMUM_HAMMING_DISTANCE;
	
	static ComparableElement<String, Double>[] apply(int numWanted, ComparableElement<String, Double>[] kmersWithWeights, HTS_KmerFilter type) throws Exception {
		
		ComparableElement<String, Double>[] kmerSequenceStatistic = null;
		Arrays.sort(kmersWithWeights);
		
		// find top kmers
		if(numWanted > kmersWithWeights.length){
			numWanted = kmersWithWeights.length;
		}
		kmerSequenceStatistic = new ComparableElement[numWanted];
		Sequence[] previouslyFoundKmers = new Sequence[numWanted];
		int nextKmerIndex=kmerSequenceStatistic.length-1;
		
		switch (type){
		case MINIMUM_HAMMING_DISTANCE:
			outerloop:
				for(int i=kmersWithWeights.length-1;i>=0;i--){
					String s = kmersWithWeights[i].getElement();
					Sequence potentialTopKmer = Sequence.create( DNAAlphabetContainer.SINGLETON, s);
					
					for(int j=kmerSequenceStatistic.length-1;j>nextKmerIndex;j--){
						Sequence previouslyFoundKmer = previouslyFoundKmers[j];
						if(DimontTool.getMinimumHammingDistance( potentialTopKmer, previouslyFoundKmer ) < 2){
							continue outerloop;
						}
					}
					kmerSequenceStatistic[nextKmerIndex] = kmersWithWeights[i];
					previouslyFoundKmers[nextKmerIndex] = potentialTopKmer;
					nextKmerIndex--;
					if(nextKmerIndex < 0){
						break;
					}
				}
			break;
			
		case HUDDINGE_DISTANCE:
			outerloop:
				for(int i=kmersWithWeights.length-1;i>=0;i--){
					String s = kmersWithWeights[i].getElement();
					Sequence potentialTopKmer = Sequence.create( DNAAlphabetContainer.SINGLETON, s);
					
					for(int j=kmerSequenceStatistic.length-1;j>nextKmerIndex;j--){
						Sequence previouslyFoundKmer = previouslyFoundKmers[j];
						if(getHuddingeDistance( potentialTopKmer, previouslyFoundKmer ) < 2){ // !
							continue outerloop;
						}
					}
					kmerSequenceStatistic[nextKmerIndex] = kmersWithWeights[i];
					previouslyFoundKmers[nextKmerIndex] = potentialTopKmer;
					nextKmerIndex--;
					if(nextKmerIndex < 0){
						break;
					}
				}
			break;
		}
						
		return kmerSequenceStatistic;
	}
	
	private static int getHuddingeDistance(Sequence seq1, Sequence seq2) throws OperationNotSupportedException {	
		// set seq1 to the longer sequence of both
		if (seq2.getLength()>seq1.getLength()) {
			Sequence helper = seq1;
			seq1  = seq2;
			seq2 = helper;
		}
		
		Sequence revComplement1 = seq1.reverseComplement();
		int minHammingDist = seq2.getLength();
		int lengthDiff = seq1.getLength() - seq2.getLength();
		
		// find minimum Hamming distance
		// part I: align seq2 without shifts 
		for (int start1=0; start1<=lengthDiff && minHammingDist>0; start1++) {
			minHammingDist = helper(seq1, seq2, start1, 0, minHammingDist);
			minHammingDist = helper(revComplement1, seq2, start1, 0, minHammingDist);
		}
		// part II: shift and align seq2
		for (int shift=1; shift<minHammingDist; shift++) {
			// left shift of seq2
			minHammingDist = shift + helper(seq1, seq2, 0, shift, minHammingDist-shift);
			minHammingDist = shift + helper(revComplement1, seq2, 0, shift, minHammingDist-shift);
			// right shift of seq2
			minHammingDist = shift + helper(seq1, seq2, lengthDiff + shift, 0, minHammingDist-shift);
			minHammingDist = shift + helper(revComplement1, seq2, lengthDiff + shift, 0, minHammingDist-shift);			
		}
				
		return minHammingDist + lengthDiff;
	}
	
	private static int helper(Sequence seq1, Sequence seq2, int start1, int start2, int currentDist) {
		int overlapLength = Math.min( seq1.getLength()-start1, seq2.getLength()-start2 );
		int distance = 0;
		for (int i=0; i<overlapLength && distance<currentDist; i++) {
			distance += (seq1.continuousVal( start1+i ) == seq2.continuousVal( start2+i ) ? 0 : 1);
		}
		return distance;
	}

}
