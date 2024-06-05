package projects.dimont.hts;

import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map.Entry;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModelFactory;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.ToolBox;

public enum HTS_KmerRanking {
	
	Z_SCORE,
	DIMONT_SCORE,
	WEIGHTED_FREQUENCY;
	
	static ComparableElement<String, Double>[] compute(int k, DataSet data, double[] weights, HTS_KmerRanking type)
			throws Exception {
		
		ComparableElement<String, Double>[] kmersWithWeights = null;

		if (type == Z_SCORE) {
			//create homogeneous Markov model of order 2
			TrainableStatisticalModel hmm = TrainableStatisticalModelFactory.createHomogeneousMarkovModel( data.getAlphabetContainer(), 4, (byte)2 );
			//train it on the input data
			hmm.train( data );
			
			// count ALL kmer occurrences
			Hashtable<String, Double> weightedKmerCounts = new Hashtable<String, Double>();
					
			//run over all sequences
			int lastStart, startPosition, n;
			Sequence seq;
			String[] s = new String[2];
			for( n = 0; n < weights.length; n++ ) {
				if (weights[n]==0)
					continue;
				seq = data.getElementAt( n );
				s[0] = seq.toString();
				s[1] = seq.reverseComplement().toString();
				lastStart = seq.getLength()-k;
				
				//run over all k-mers
				for( startPosition = 0; startPosition <= lastStart; startPosition++ ) {
					String kmer = s[0].substring( startPosition, startPosition+k );
					String revComplKmer = s[1].substring( s[0].length()-k-startPosition, s[0].length()-startPosition );
					kmer = kmer.compareTo(revComplKmer) < 0 ? kmer : revComplKmer;
					if( !weightedKmerCounts.containsKey( kmer ) ) {
						weightedKmerCounts.put( kmer, weights[n] );
					}
					else {
						weightedKmerCounts.put( kmer, weightedKmerCounts.get(kmer)+weights[n]);
					}
					
				}
			}
			kmersWithWeights = new ComparableElement[weightedKmerCounts.size()];
			Iterator<Entry<String, Double>> it = weightedKmerCounts.entrySet().iterator(); 
			Entry<String, Double> weightedKmerCount;
			double numberOfPatterns = ToolBox.sum(weights) * (data.getElementAt(0).getLength()-k+1) ;
			
			String key;
			double value, expectedValue_bg;
			for( int i = 0; i < kmersWithWeights.length; i++ ) {
				weightedKmerCount = it.next();
				key = weightedKmerCount.getKey();
				seq = Sequence.create(data.getAlphabetContainer(), key);
				value = weightedKmerCount.getValue();
				expectedValue_bg = numberOfPatterns*(Math.exp(hmm.getLogProbFor(seq))+Math.exp(hmm.getLogProbFor(seq.reverseComplement())));
				kmersWithWeights[i] = new ComparableElement<String, Double>( key , (value-expectedValue_bg)/Math.sqrt(expectedValue_bg) );
			}
				
		}
		else {
			// count kmers
			Hashtable<String, double[]> weightedKmerCounts = getWeightedKmerCounts( k, data, weights );
			// sort kmers by default ( log(n_fg)*n_fg/n_bg ) or enrichment
			kmersWithWeights = new ComparableElement[weightedKmerCounts.size()];
			Iterator<Entry<String, double[]>> it = weightedKmerCounts.entrySet().iterator();
			Entry<String, double[]> weightedKmerCount;
			
			
			if (type == WEIGHTED_FREQUENCY) {
				for( int i = 0; i < kmersWithWeights.length; i++ ) {
					weightedKmerCount = it.next();
					kmersWithWeights[i] = new ComparableElement<String, Double>( weightedKmerCount.getKey(), weightedKmerCount.getValue()[0] );
				}
			}
			else if (type == DIMONT_SCORE) { // DEFAULT
				for( int i = 0; i < kmersWithWeights.length; i++ ) {
					weightedKmerCount = it.next();
					double[] val = weightedKmerCount.getValue();
					kmersWithWeights[i] = new ComparableElement<String, Double>( weightedKmerCount.getKey(), Math.log(val[0]+1)*(val[0]+1)/(val[1]+1) );
					
				}
			}
		}
		
		
		return(kmersWithWeights);
	}
	
	private static Hashtable<String, double[]> getWeightedKmerCounts(int k, DataSet data, double[] weights) throws WrongAlphabetException, OperationNotSupportedException{
		AlphabetContainer con = data.getAlphabetContainer();
		if( !con.isSimple() || !con.isDiscrete() ) {
			throw new WrongAlphabetException();
		}
		Hashtable<String, double[]> weightedKmerCounts = new Hashtable<String, double[]>();
		HashSet<String> kmersOfSingleSequence = new HashSet<String>();
				
		//run over all sequences
		int lastStart, startPosition, n;
		Sequence seq;
		String[] s = new String[2];
		for( n = 0; n < weights.length; n++ ) {
			seq = data.getElementAt( n );
			s[0] = seq.toString();
			s[1] = seq.reverseComplement().toString();
			lastStart = seq.getLength()-k;
			
			//run over all k-mers
			kmersOfSingleSequence.clear();
			for( startPosition = 0; startPosition <= lastStart; startPosition++ ) {
				String kmer = s[0].substring( startPosition, startPosition+k );
				String revComplKmer = s[1].substring( s[0].length()-k-startPosition, s[0].length()-startPosition );
				kmer = kmer.compareTo(revComplKmer) < 0 ? kmer : revComplKmer;
				if( !kmersOfSingleSequence.contains( kmer ) ) {
					kmersOfSingleSequence.add( kmer );
				}
				
			}
			
			Iterator<String> it = kmersOfSingleSequence.iterator();
			double[] kmerCount;
			while( it.hasNext() ) {
				s[0] = it.next();
				if( weightedKmerCounts.containsKey( s[0] ) ) {
					kmerCount = weightedKmerCounts.get( s[0] );
					kmerCount[0] += weights[n];
					kmerCount[1] += 1d-weights[n];
				} else {
					weightedKmerCounts.put( s[0], new double[]{weights[n],1d-weights[n]} );
				}			
			}
		}
		return weightedKmerCounts;
	}

}
