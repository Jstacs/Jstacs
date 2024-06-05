package projects.dimont.hts;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map.Entry;

import javax.naming.OperationNotSupportedException;

import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.ToolBox.TiedRanks;
import projects.dimont.Interpolation;

public enum HTS_SequenceWeights {
	
	SEQUENCE_SPECIFIC_RMAX,
	NOT_HTSELEX,
	//SIMPLE_APPROACH,
	CYCLE_SPECIFIC,
	SEQUENCE_SPECIFIC_SUM,
	SEQUENCE_SPECIFIC_MAX,
	SEQUENCE_SPECIFIC_RSUM,
	SEQUENCE_SPECIFIC_RPROD,
	SEQUENCE_SPECIFIC_LOGSUM,
	CYCLE_SPECIFIC_SMOOTHED;
	
	static double[] compute(DataSet seqs, double[] signals, double wf, HTS_SequenceWeights type, String[] protocolString) throws Exception{
		double[] weights = new double[signals.length];
		double[] weightsPerCycle;
		int minCycle, maxCycle;
			      
		switch(type) {
		case NOT_HTSELEX:
			weights = Interpolation.getWeight( seqs, signals, wf, Interpolation.RANK_LOG );
			break;

		case CYCLE_SPECIFIC:
			weightsPerCycle = getAllFractionsOfNonspecificSequences(seqs, signals);
			minCycle = (int) ToolBox.min(signals);
			
			weightsPerCycle = getCycleSpecificWeightsFromFractions(weightsPerCycle, minCycle, protocolString);
			
			// assign cycle weights
			for (int n=0;n<signals.length;n++) {
				weights[n] = weightsPerCycle[(int) signals[n] - minCycle];
			}
			break;
			
		case CYCLE_SPECIFIC_SMOOTHED:
			weightsPerCycle = getAllFractionsOfNonspecificSequences(seqs, signals);
			minCycle = (int) ToolBox.min(signals);
			
			weightsPerCycle = getSmoothedCycleSpecificWeightsFromFractions(weightsPerCycle, minCycle, protocolString);
			
			// assign cycle weights
			for (int n=0;n<signals.length;n++) {
				weights[n] = weightsPerCycle[(int) signals[n] - minCycle];
			}
			break;
			
		case SEQUENCE_SPECIFIC_SUM:
			weights = getAllSequenceSpecificWeights(seqs, signals, SEQUENCE_SPECIFIC_SUM, protocolString);
			break;
			
		case SEQUENCE_SPECIFIC_RPROD:
			weights = getAllSequenceSpecificWeights(seqs, signals, SEQUENCE_SPECIFIC_RPROD, protocolString);
			break;
			
		case SEQUENCE_SPECIFIC_MAX:
			weights = getAllSequenceSpecificWeights(seqs, signals, SEQUENCE_SPECIFIC_MAX, protocolString);
			break;
			
		case SEQUENCE_SPECIFIC_RMAX:
			weights = getAllSequenceSpecificWeights(seqs, signals, SEQUENCE_SPECIFIC_RMAX, protocolString);
			break;
			
		case SEQUENCE_SPECIFIC_RSUM:
			weights = getAllSequenceSpecificWeights(seqs, signals, SEQUENCE_SPECIFIC_RSUM, protocolString);
			break;
			
		case SEQUENCE_SPECIFIC_LOGSUM:
			weights = getAllSequenceSpecificWeights(seqs, signals, SEQUENCE_SPECIFIC_LOGSUM, protocolString);
			break;
	}	
		
		return weights;
	}
	
	private static double[] getCycleSpecificWeightsFromFractions(double[] fractions, int minCycle, String[] protocolString) {
		double[] weightsPerCycle = new double[fractions.length];
		// correction of "fractions"
		// -> later cycles are known to have a lower fraction of nonspecific sequences
		weightsPerCycle[0] = fractions[0];
		for (int cycle=0; cycle<weightsPerCycle.length-1; cycle++) {
			weightsPerCycle[cycle+1] = fractions[cycle+1];
			// if fraction_k < fraction_k+1
			if (weightsPerCycle[cycle] < weightsPerCycle[cycle+1])
				weightsPerCycle[cycle+1] = weightsPerCycle[cycle];
		}
		
		// compute cycle weights
		for (int i=weightsPerCycle.length-1;i>=0;i--) {
			weightsPerCycle[i] = 1.0 - weightsPerCycle[i]/weightsPerCycle[0];
		}
		
		protocolString[0] = protocolString[0] + "Cycle weights\n";
		if (ToolBox.sum(weightsPerCycle)==0 ) { // enrichment experiment failed
			weightsPerCycle[0] = 0.0;
			protocolString[0] = protocolString[0] +("cycle "+(minCycle)+": "+weightsPerCycle[0])+"\n";
			for (int i=1;i<weightsPerCycle.length;i++) {
				weightsPerCycle[i] = 1.0;
				protocolString[0] = protocolString[0] +("cycle "+(minCycle+i)+": "+weightsPerCycle[i])+"\n";
			}
		}
		else {
			for (int i=0;i<weightsPerCycle.length;i++) {
				weightsPerCycle[i] /= weightsPerCycle[weightsPerCycle.length-1];
				protocolString[0] = protocolString[0] +("cycle "+(minCycle+i)+": "+weightsPerCycle[i])+"\n";
			}
		}
		return weightsPerCycle;
	}
	
	private static double[] getSmoothedCycleSpecificWeightsFromFractions(double[] fractions, int minCycle, String[] protocolString) {
		double[] weightsPerCycle = new double[fractions.length];
		double maxWeight = ToolBox.max(fractions);
		// pre-compute cycle weights
		for (int i=0;i<weightsPerCycle.length;i++) {
			weightsPerCycle[i] = 1.0 - fractions[i]/maxWeight;
		}
		maxWeight = ToolBox.max(weightsPerCycle);
		for (int i=0;i<weightsPerCycle.length;i++) {
			weightsPerCycle[i] /= maxWeight ;
		}
		
		double enrichmentConstant = approximateEnrichment(weightsPerCycle, minCycle, Math.pow(10, -9));
		protocolString[0] = protocolString[0] +"Enrichment constant: "+enrichmentConstant+"\n";
		
		protocolString[0] = protocolString[0] + "Cycle weights\n";
		int maxCycle = minCycle + weightsPerCycle.length -1;
		for (int i=minCycle;i<=maxCycle;i++) {
			weightsPerCycle[i-minCycle] = Math.pow(enrichmentConstant, i-maxCycle);
			protocolString[0] = protocolString[0] +("cycle "+i+": "+weightsPerCycle[i-minCycle])+"\n";
		}
		return weightsPerCycle;
	}
	
	static double[] getAllFractionsOfNonspecificSequences(DataSet seqs, double[] signals) throws WrongAlphabetException, OperationNotSupportedException {
		Pair<Double,Integer>[] informationPerCycle = getNonspecificSequenceInformation( seqs, signals );
		double[] fractions = new double[informationPerCycle.length];
		for (int i=0;i<informationPerCycle.length;i++) {
			fractions[i] = informationPerCycle[i].getFirstElement() / informationPerCycle[i].getSecondElement();
		}
		return fractions;
	}
	
	private static double[] getAllFractionsOfNonspecificSequences(Pair<Double,Integer>[] informationPerCycle) throws WrongAlphabetException, OperationNotSupportedException {
		double[] fractions = new double[informationPerCycle.length];
		for (int i=0;i<informationPerCycle.length;i++) {
			fractions[i] = informationPerCycle[i].getFirstElement() / informationPerCycle[i].getSecondElement();
		}
		return fractions;
	}

	private static Pair<Double,Integer>[] getNonspecificSequenceInformation(DataSet seqs, double[] signals) throws WrongAlphabetException, OperationNotSupportedException {
		AlphabetContainer con = seqs.getAlphabetContainer();
		if( !con.isSimple() || !con.isDiscrete() ) {
			throw new WrongAlphabetException();
		}
		
		// sort sequences by cycle number
		ComparableElement<Sequence, Integer>[] data = new ComparableElement[signals.length];
		for (int n=0;n<signals.length;n++) {
			data[n] = new ComparableElement<Sequence, Integer>( seqs.getElementAt( n ), (int) signals[n] );
		}
		Arrays.sort(data);
		int minCycle = data[0].getWeight();
		int maxCycle = data[signals.length-1].getWeight();
		
		// compute "fraction" of any set of nonspecific sequences out of total sequences for each cycle
		Pair<Double,Integer>[] informationPerCycle = new Pair[maxCycle-minCycle+1];
		int cycle = minCycle;
		ArrayList<Sequence> cycleSequences = new ArrayList<Sequence>();
		for (int n=0;n<signals.length;n++) {
			if (data[n].getWeight()>cycle) { // next cycle
				Sequence[] cycleData = cycleSequences.toArray( new Sequence[cycleSequences.size()] );
				cycleSequences.clear();
				informationPerCycle[cycle-minCycle] = getSequenceInformationOfCycle( cycleData ).getSecondElement();
				cycle++;
			}
			else { // current cycle
				cycleSequences.add(data[n].getElement());
			}
		}
		Sequence[] cycleData = cycleSequences.toArray( new Sequence[cycleSequences.size()] );
		informationPerCycle[maxCycle-minCycle] = getSequenceInformationOfCycle( cycleData ).getSecondElement();
		
		return informationPerCycle;
	}
	
	private static double[] getAllSequenceSpecificWeights(DataSet seqs, double[] signals, HTS_SequenceWeights type, String[] protocolString) throws WrongAlphabetException, OperationNotSupportedException {
		AlphabetContainer con = seqs.getAlphabetContainer();
		if( !con.isSimple() || !con.isDiscrete() ) {
			throw new WrongAlphabetException();
		}
		
		// sort sequences by cycle number
		int[] dataSetIndices = ToolBox.order(signals, false);
		ComparableElement<Sequence, Integer>[] data = new ComparableElement[signals.length];
		for (int n=0;n<signals.length;n++) {
			data[dataSetIndices[n]] = new ComparableElement<Sequence, Integer>( seqs.getElementAt( n ), (int) signals[n] );
		}
		int minCycle = data[0].getWeight();
		int maxCycle = data[signals.length-1].getWeight();
		ArrayList<Sequence[]> dataPerCycle = new ArrayList<Sequence[]>(maxCycle-minCycle+1);
		
		// group data by cycle number and compute non specific sequence information
		Pair<ComparableElement<String,Double>[],Pair<Double,Integer>>[] informationPerCycle = new Pair[maxCycle-minCycle+1];
		int cycle = minCycle;
		ArrayList<Sequence> cycleSequences = new ArrayList<Sequence>();
		for (int n=0;n<signals.length;n++) {
			if (data[n].getWeight()>cycle) { // next cycle
				Sequence[] cycleData = cycleSequences.toArray( new Sequence[cycleSequences.size()] );
				dataPerCycle.add(cycleData);
				cycleSequences.clear();
				informationPerCycle[cycle-minCycle] = getSequenceInformationOfCycle( cycleData );
				cycle++;
			}
			else { // current cycle
				cycleSequences.add(data[n].getElement());
			}
		}
		Sequence[] cycleData = cycleSequences.toArray( new Sequence[cycleSequences.size()] );
		dataPerCycle.add(cycleData);
		informationPerCycle[maxCycle-minCycle] = getSequenceInformationOfCycle( cycleData );
		
		// compute "fraction" of any set of nonspecific sequences out of total sequences for each cycle
		Pair<Double,Integer>[] fractionInfo = new Pair[informationPerCycle.length];
		ArrayList<ComparableElement<String,Double>[]> specific8mers = new ArrayList<ComparableElement<String,Double>[]>(informationPerCycle.length);
		for (int i=0;i<informationPerCycle.length;i++) {
			fractionInfo[i] = informationPerCycle[i].getSecondElement();
			specific8mers.add(informationPerCycle[i].getFirstElement());
		}
		double[] fractions = getAllFractionsOfNonspecificSequences(fractionInfo);
		// compute cycle specific weights
		double[] weightsPerCycle = getCycleSpecificWeightsFromFractions(fractions, minCycle, protocolString);
		
		// compute sequence specific weights
		ArrayList<double[]> sequenceSpecificWeightsPerCycle = getSequenceSpecificWeightsPerCycle(dataPerCycle, specific8mers, weightsPerCycle, type, protocolString);
		double[] sequenceSpecificWeights = new double[signals.length];
		int firstCycleIndex = 0;
		for (int i=0;i<sequenceSpecificWeightsPerCycle.size();i++) {
			for (int j=0; j<sequenceSpecificWeightsPerCycle.get(i).length; j++) {
				sequenceSpecificWeights[firstCycleIndex+j] = sequenceSpecificWeightsPerCycle.get(i)[j];
			}
			firstCycleIndex += sequenceSpecificWeightsPerCycle.get(i).length;
		}
		double[] weights = new double[signals.length];
		for (int n=0;n<signals.length;n++) {
			weights[n] = sequenceSpecificWeights[dataSetIndices[n]];
		}
		
		return weights;
	}
	
	private static Pair<ComparableElement<String,Double>[],Pair<Double,Integer>> getSequenceInformationOfCycle(Sequence[] cycleData) throws WrongAlphabetException, OperationNotSupportedException{
		// count all 8mers
		Hashtable<String, Double> kmerCounts = new Hashtable<String, Double>();
		HashSet<String> kmersOfSingleSequence = new HashSet<String>();
		int k = 8;
				
		// run over all sequences
		int lastStart, startPosition, n;
		Sequence seq;
		String[] s = new String[2];
		for( n = 0; n < cycleData.length ; n++ ) {
			seq = cycleData[n];
			s[0] = seq.toString();
			s[1] = seq.reverseComplement().toString();
			lastStart = seq.getLength()-k;
			
			// run over all k-mers
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
			double kmerCount;
			while( it.hasNext() ) {
				s[0] = it.next();
				kmerCount = kmerCounts.containsKey( s[0] )? kmerCounts.get( s[0] )+1 : 1.0;
				kmerCounts.put( s[0], kmerCount );	
			}			
		}
		
		// sum occurrence counts of kmers that rank between 25% and 75% in relative abundance
		// and remember all interesting kmers that rank 75% or more in relative abundance
		int first_index = (int) Math.ceil( 0.25*kmerCounts.size() );
		int last_index = (int) Math.floor( 0.75*kmerCounts.size() );
		
		ComparableElement<String,Double>[] comparableKmerCounts = new ComparableElement[kmerCounts.size()];
		Iterator<Entry<String, Double>> it = kmerCounts.entrySet().iterator();
		Entry<String, Double> entry;
		for( int a = 0; a < comparableKmerCounts.length; a++ ) {
			entry = it.next();
			comparableKmerCounts[a] = new ComparableElement(entry.getKey(),entry.getValue());
		}
		Arrays.sort(comparableKmerCounts);
		
		double nonspecificOccurrences = 0.0;
		ComparableElement<String,Double>[] interestingKmers = new ComparableElement[kmerCounts.size()-last_index-1];
		for (int i=first_index; i<=last_index; i++) {
			nonspecificOccurrences += comparableKmerCounts[i].getWeight();
		}
		for (int i=last_index+1; i<kmerCounts.size(); i++) {
			interestingKmers[i-last_index-1] = comparableKmerCounts[i];
		}
		// return specific sequences and nonspecific occurrence number and number of sequences
		return new Pair<ComparableElement<String,Double>[],Pair<Double,Integer>>(interestingKmers,
				new Pair<Double,Integer>(nonspecificOccurrences,cycleData.length)
				);
	}

	private static ArrayList<double[]> getSequenceSpecificWeightsPerCycle(ArrayList<Sequence[]> cycleData, 
			ArrayList<ComparableElement<String,Double>[]> topKmers, double[] cycleWeights, HTS_SequenceWeights type, String[] protocolString) throws OperationNotSupportedException {
		
		ArrayList<double[]> sequenceSpecificWeightsPerCycle = new ArrayList<double[]>(cycleWeights.length);
		double[] cycleWeights2 = cycleWeights.clone();
		
		for (int cycleIndex=cycleWeights.length-1; cycleIndex>=0; cycleIndex--) {
			
			// count specific kmers per sequence
			double[] sequenceWeights = new double[cycleData.get(cycleIndex).length];
			if (cycleWeights2[cycleIndex] == 0.0) {
				sequenceSpecificWeightsPerCycle.add(sequenceWeights);
				continue;
			}
			double[] counter = new double[cycleData.get(cycleIndex).length];
			double[] sums = new double[cycleData.get(cycleIndex).length];
			double[] rSums = new double[cycleData.get(cycleIndex).length];
			double[] logSums = new double[cycleData.get(cycleIndex).length];
			double[] rProds = new double[cycleData.get(cycleIndex).length];
			double[] maxs = new double[cycleData.get(cycleIndex).length];
			double[] rMaxs = new double[cycleData.get(cycleIndex).length];
			double curMax;
			double curRMin;
			
			double[] counts = new double[topKmers.get(cycleIndex).length];
			for (int i=0;i<topKmers.get(cycleIndex).length;i++) {
				counts[i] = topKmers.get(cycleIndex)[i].getWeight();
			}
			int[] ranks = ToolBox.rank( counts, TiedRanks.CONTIGUOUS ); // greatest value obtains the lowest rank (0)
			int maxRank = 0;

			Hashtable<String, Double> topKmerSet = new Hashtable(topKmers.get(cycleIndex).length, (float) 1);
			Hashtable<String, Double> topKmerSetWithRanks = new Hashtable(topKmers.get(cycleIndex).length, (float) 1);
			for (int i=0;i<topKmers.get(cycleIndex).length;i++) {
				topKmerSet.put(topKmers.get(cycleIndex)[i].getElement(),topKmers.get(cycleIndex)[i].getWeight());
				topKmerSetWithRanks.put(topKmers.get(cycleIndex)[i].getElement(), (double) ranks[i]);
				if (maxRank < ranks[i])
					maxRank = ranks[i];
			}
			int k = 8;
					
			// run over all sequences
			int lastStart, startPosition, n;
			double rank, count;
			Sequence seq;
			String[] s = new String[2];
			for( n = 0; n < cycleData.get(cycleIndex).length ; n++ ) {
				seq = cycleData.get(cycleIndex)[n];
				s[0] = seq.toString();
				s[1] = seq.reverseComplement().toString();
				lastStart = seq.getLength()-k;
				curMax = 0;
				curRMin = Double.MAX_VALUE;
				rProds[n] = 1.0;
				
				// run over all k-mers
				for( startPosition = 0; startPosition <= lastStart; startPosition++ ) {
					String kmer = s[0].substring( startPosition, startPosition+k );
					String revComplKmer = s[1].substring( s[0].length()-k-startPosition, s[0].length()-startPosition );
					kmer = kmer.compareTo(revComplKmer) < 0 ? kmer : revComplKmer;
					if( topKmerSet.containsKey( kmer ) ) {
						counter[n] += 1;
						count = topKmerSet.get(kmer);
						rank = topKmerSetWithRanks.get(kmer);
						sums[n] += count;
						logSums[n] += Math.log(count);
						rSums[n] += rank;
						rProds[n] *= (maxRank - rank +1);
						if (count > curMax) {
							curMax = count;
						}
						if (rank < curRMin) {
							curRMin = rank;
						}
					}
				}
				rProds[n] = (curRMin == Double.MAX_VALUE)? 0 : Math.pow(rProds[n], 1.0/counter[n]);
				maxs[n] = curMax;
				rMaxs[n] = (curRMin == Double.MAX_VALUE)? 0 : maxRank - curRMin +1;
			}

			
			double cycleDifference, prevWeight;
			if (cycleIndex == 0) {
				cycleDifference = cycleWeights2[0];
				prevWeight = 0;
			}
			else {
				cycleDifference = cycleWeights2[cycleIndex]-cycleWeights2[cycleIndex-1];
				cycleDifference = Math.min(0.5, cycleDifference);
				prevWeight = cycleWeights2[cycleIndex] - cycleDifference;
				if (cycleDifference==0.0) {
					cycleDifference = Math.min(0.01, cycleWeights2[cycleIndex-1]/2);
					cycleWeights2[cycleIndex-1] -= cycleDifference;
					// correction of cycleWeights2
					// -> later cycles are known to have a higher (or equal) weight
					for (int cycle=cycleIndex-1; cycle>0; cycle--) {
						if (cycleWeights2[cycle-1] > cycleWeights2[cycle])
							cycleWeights2[cycle-1] = cycleWeights2[cycle];
					}
					prevWeight = cycleWeights2[cycleIndex-1];
				}
			}
			
			if (type == SEQUENCE_SPECIFIC_SUM) {
				double maxSum = ToolBox.max(sums);
				for( n = 0; n < cycleData.get(cycleIndex).length ; n++ ) {
					sequenceWeights[n] = prevWeight + cycleDifference*sums[n]/maxSum ;
				}
			}
			else if (type == SEQUENCE_SPECIFIC_MAX) {
				double maxMax = ToolBox.max(maxs);
				for( n = 0; n < cycleData.get(cycleIndex).length ; n++ ) {
					sequenceWeights[n] = prevWeight + cycleDifference*maxs[n]/maxMax ;
				}
			}
			else if (type == SEQUENCE_SPECIFIC_RPROD) {
				double maxRProd = ToolBox.max(rProds);
				for( n = 0; n < cycleData.get(cycleIndex).length ; n++ ) {
					sequenceWeights[n] = prevWeight + cycleDifference*rProds[n]/maxRProd ;
				}
			}
			else if (type == SEQUENCE_SPECIFIC_RSUM) {
				double maxRSum = ToolBox.max(rSums);
				for( n = 0; n < cycleData.get(cycleIndex).length ; n++ ) {
					sequenceWeights[n] = prevWeight + cycleDifference*rSums[n]/maxRSum ;
				}
			}
			else if (type == SEQUENCE_SPECIFIC_RMAX) {
				double maxRMax = ToolBox.max(rMaxs);
				for( n = 0; n < cycleData.get(cycleIndex).length ; n++ ) {
					sequenceWeights[n] = prevWeight + cycleDifference*rMaxs[n]/maxRMax ;
				}
			}
			else if (type == SEQUENCE_SPECIFIC_LOGSUM) {
				double maxLogSum = ToolBox.max(logSums);
				for( n = 0; n < cycleData.get(cycleIndex).length ; n++ ) {
					sequenceWeights[n] = prevWeight + cycleDifference*logSums[n]/maxLogSum ;
				}
			}
			
			
			sequenceSpecificWeightsPerCycle.add(sequenceWeights);
		}
		
		protocolString[0] = protocolString[0] + "-----------------------------------------\n" + "Effective cycle weights:\n";
		for (int i=0;i<cycleWeights2.length;i++) {
			protocolString[0] = protocolString[0] +cycleWeights2[i]+"\n";
		}
		
		Collections.reverse(sequenceSpecificWeightsPerCycle);
		return sequenceSpecificWeightsPerCycle;
	}
	
	private static double approximateEnrichment(double[] weightsPerCycle, int minCycle, double precision) {
		// argmin_e  sum_cycle ( e^(cycle-maxCyle)-weightsPerCycle[cycle] )^2
		double old_e = 1.0;
		double new_e;
		double direction = getDirection(old_e, weightsPerCycle, minCycle);
		double stepsize = 1.0;
		
		while (true) {
			new_e = old_e + direction*stepsize;
			direction = getDirection(new_e, weightsPerCycle, minCycle);
			if (direction == 0) break;
			if ((old_e<new_e && direction<0) || (old_e>new_e && direction>0)) { //change in direction
				if (stepsize>precision) stepsize *= 0.1;
				else break;
			}
			if (old_e == new_e || new_e >= 50.0 ) break;
			old_e = new_e;
		}
		return new_e;
	}
	
	private static double getDirection(double e, double[] weightsPerCycle, int minCycle) {
		int maxCycle = minCycle + weightsPerCycle.length -1;
		double derivation = 0.0;
		for (int i=minCycle; i<=maxCycle; i++) {
			derivation += (2 *Math.pow(e, i-maxCycle) *(i-maxCycle) *(Math.pow(e, i-maxCycle)-weightsPerCycle[i-minCycle]) ) /e;
		}
		return - Math.signum(derivation);
	}
	
	static double findPwmCorrectionFactor(int cycle1, int cycle2, DataSet seqs, double[] signals) throws OperationNotSupportedException, WrongAlphabetException {
		if (cycle2<cycle1) {
			int helper = cycle1;
			cycle1  = cycle2;
			cycle2 = helper;
		}
		
		int minCycle = (int) ToolBox.min(signals);
		
		Pair<Double,Integer>[] information = getNonspecificSequenceInformation(seqs, signals);
		double[] fractions = new double[information.length];
		
		// compute "fraction" of any set of nonspecific sequences out of total sequences for each cycle
		for (int i=0;i<fractions.length;i++) {
			fractions[i] = information[i].getFirstElement() / information[i].getSecondElement();
		}
		
		// compute pwm correction factor
		double factor = fractions[cycle2-minCycle] / fractions[cycle1-minCycle];
		if (factor>1) factor = 1.0;
		factor *= ((double) information[cycle2-minCycle].getSecondElement()  / information[cycle1-minCycle].getSecondElement() );
		
		return factor;
	}

}
