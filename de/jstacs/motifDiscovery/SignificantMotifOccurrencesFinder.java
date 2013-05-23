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

package de.jstacs.motifDiscovery;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.classifiers.utils.PValueComputation;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.ComplementableDiscreteAlphabet;
import de.jstacs.data.sequences.PermutedSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.MotifAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile;
import de.jstacs.results.NumericalResult;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;


/**
 * This class enables the user to predict motif occurrences given a specific significance level. 
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class SignificantMotifOccurrencesFinder {

	
	/**
	 * 
	 * @author Jan Grau, Jens Keilwagen
	 */
	public enum RandomSeqType{

		/**
		 * A enum constant that indicates to use the a background set to determine the significance level.
		 */
		BACKGROUND(-2),
		/**
		 * A enum constant that indicates to use permuted instances of the sequence to determine the significance level.
		 */
		PERMUTED(-1),
		/**
		 * A enum constant that indicates to use sequences drawn from a homogeneous Markov model of order 0 to determine the significance level.
		 */
		hMM0(0),
		/**
		 * A enum constant that indicates to use sequences drawn from a homogeneous Markov model of order 1 to determine the significance level.
		 */
		hMM1(1),
		/**
		 * A enum constant that indicates to use sequences drawn from a homogeneous Markov model of order 3 to determine the significance level.
		 */
		hMM2(2),
		/**
		 * A enum constant that indicates to use sequences drawn from a homogeneous Markov model of order 3 to determine the significance level.
		 */
		hMM3(3),
		/**
		 * A enum constant that indicates to use sequences drawn from a homogeneous Markov model of order 4 to determine the significance level.
		 */
		hMM4(4),
		/**
		 * A enum constant that indicates to use sequences drawn from a homogeneous Markov model of order 5 to determine the significance level.
		 */
		hMM5(5);
		
		private final int order;
		
		RandomSeqType(int order){
			this.order = order;
		}
		
		/**
		 * This method returns the Markov order.
		 * 
		 * @return the Markov order
		 */
		public int getOrder(){
			return order;
		}
		
	};
	
	private RandomSeqType type;
	private boolean oneHistogram;
	private DataSet bg;
	private double[] weights;
	private MotifDiscoverer disc;
	private int numSequences;
	private double sign;
	private double[][] sortedScores;
	private int[][][] globalToLocalIndexes;
	private JoinMethod joinMethod;
	
	/**
	 * This constructor creates an instance of {@link SignificantMotifOccurrencesFinder} that uses the given {@link RandomSeqType} to determine the siginificance level.
	 * 
	 * @param disc the {@link MotifDiscoverer} for the prediction
	 * @param type the type that determines how the significance level is determined 
	 * @param oneHistogram a switch to decide whether to use one background distribution histogram for all sequence or sequence specific background distribution histograms
	 * @param numSequences the number of sampled sequence instances used to determine the significance level
	 * @param sign the significance level
	 */
	public SignificantMotifOccurrencesFinder(MotifDiscoverer disc, RandomSeqType type, boolean oneHistogram, int numSequences, double sign){
		this(disc,type,new SumOfProbabilities(),oneHistogram,numSequences,sign);
	}
	/**
	 * This constructor creates an instance of {@link SignificantMotifOccurrencesFinder} that uses the given {@link RandomSeqType} to determine the siginificance level.
	 * 
	 * @param disc the {@link MotifDiscoverer} for the prediction
	 * @param type the type that determines how the significance level is determined 
	 * @param joiner the {@link JoinMethod} that defines how the profiles of the same motif in different components shall be joined
	 * @param oneHistogram a switch to decide whether to use one background distribution histogram for all sequence or sequence specific background distribution histograms
	 * @param numSequences the number of sampled sequence instances used to determine the significance level
	 * @param sign the significance level
	 */
	public SignificantMotifOccurrencesFinder(MotifDiscoverer disc, RandomSeqType type, JoinMethod joiner, boolean oneHistogram, int numSequences, double sign){
		this.disc = disc;
		if( type == RandomSeqType.BACKGROUND ) {
			throw new IllegalArgumentException( "This type can not be used in this constructor." );
		}
		this.type = type;
		this.joinMethod = joiner;
		this.oneHistogram = oneHistogram;
		this.numSequences = numSequences;
		this.sign = sign;
		prepareIndices();
	}

	/**
	 * This constructor creates an instance of {@link SignificantMotifOccurrencesFinder} that uses a {@link DataSet} to determine the siginificance level.
	 * 
	 * @param disc the {@link MotifDiscoverer} for the prediction
	 * @param bg the background data set
	 * @param weights the weights of the background data set, can be <code>null</code>
	 * @param sign the significance level
	 */
	public SignificantMotifOccurrencesFinder(MotifDiscoverer disc, DataSet bg, double[] weights, double sign){
		this(disc,new SumOfProbabilities(),bg,weights,sign);
	}
	
	/**
	 * This constructor creates an instance of {@link SignificantMotifOccurrencesFinder} that uses a {@link DataSet} to determine the siginificance level.
	 * 
	 * @param disc the {@link MotifDiscoverer} for the prediction
	 * @param joiner the {@link JoinMethod} that defines how the profiles of the same motif in different components shall be joined
	 * @param bg the background data set
	 * @param weights the weights of the background data set, can be <code>null</code>
	 * @param sign the significance level
	 */
	public SignificantMotifOccurrencesFinder(MotifDiscoverer disc, JoinMethod joiner, DataSet bg, double[] weights, double sign){
		this.disc = disc;
		this.type = RandomSeqType.BACKGROUND;
		this.joinMethod = joiner;
		this.oneHistogram = true;
		this.numSequences = bg.getNumberOfElements();
		this.bg = bg;
		this.weights = weights==null ? null : weights.clone();
		this.sign = sign;
		prepareIndices();
	}

	
	private void createBgDataSet( DataSet s ) throws Exception {
		switch( type ) {
			case BACKGROUND:
				//already existing
				break;
			case PERMUTED:
				Sequence[] seqs = new Sequence[s.getNumberOfElements()*numSequences];
				Sequence current;
				for( int n = 0, j, i = 0; i < s.getNumberOfElements(); i++ ) {
					current = s.getElementAt( i );
					for( j = 0; j < numSequences; j++, n++ ) {
						seqs[n] = new PermutedSequence( current );
					}
				}
				bg = new DataSet( "permuted " + s.getAnnotation(), seqs );
				break;
			case hMM0:
			case hMM1:
			case hMM2:
			case hMM3:
			case hMM4:
			case hMM5:
				int order = type.getOrder();
				HomogeneousMMDiffSM hmm = new HomogeneousMMDiffSM(s.getAlphabetContainer(),order,0,new double[order+1],true,true,1);
				hmm.initializeFunction( 0, false, new DataSet[]{s}, null );
				
				if( order > 0 ) {
					double[][][] condProbs = hmm.getAllConditionalStationaryDistributions();
					DoubleList list = new DoubleList( (int) (1.5*Math.pow( s.getAlphabetContainer().getAlphabetLengthAt( 0 ), condProbs.length )) );//TODO
					for( int k, j, i = 0; i < condProbs.length; i++ ) {
						for( j = 0; j < condProbs[i].length; j++ ) {
							for( k = 0; k < condProbs[i][j].length; k++ ) {
								list.add( Math.log( condProbs[i][j][k] ) );
							}
						}
					}
					hmm.setParameters( list.toArray(), 0 );
				}
				
				bg = hmm.emit( numSequences*s.getNumberOfElements(), s.getElementLength() );
				break;
			default:
				//XXX
		}
	}
	
	private void createBgDataSet( Sequence seq ) throws Exception {
		createBgDataSet( new DataSet( "", seq ) );
	}
	
	private double[][] getAllProfilesOfScoresFor(int motif, Sequence seq, int start) throws Exception{
		double[][] profiles = new double[globalToLocalIndexes[motif][0].length][];
		for(int i=0;i<globalToLocalIndexes[motif][0].length;i++){
			profiles[i] = disc.getProfileOfScoresFor( globalToLocalIndexes[motif][0][i], globalToLocalIndexes[motif][1][i], seq, start, KindOfProfile.UNNORMALIZED_JOINT );
		}
		return profiles;
	}
	
	private double[] getJointProfileOfScoresFor(int motif, Sequence seq, int start) throws Exception{
		double[][] profiles = getAllProfilesOfScoresFor( motif, seq, start );
		return joinMethod.joinProfiles( profiles );
	}
/*	
public void test() throws Exception {
	fillSortedScoresArray(0, 0);
	System.out.println("0: " + Arrays.toString(sortedScores[0]));
	System.out.println("1: " + Arrays.toString(sortedScores[1]));
	for( int i = 0; i < sortedScores[0].length; i++ ) {
		System.out.println( i 
				+ "\t" + (sortedScores[0][i]*0.99) + "\t" + getPValue(sortedScores[0], sortedScores[0][i]*0.99, sortedScores[1])
				+ "\t" + (sortedScores[0][i]) + "\t" + getPValue(sortedScores[0], sortedScores[0][i], sortedScores[1])
				+ "\t" + (sortedScores[0][i]*1.01) + "\t" + getPValue(sortedScores[0], sortedScores[0][i]*1.01, sortedScores[1])
		);
	}
}
*/	
	private void fillSortedScoresArray( int motif, int start ) throws Exception {
//XXX von hier
		int num = 0, i;
		Sequence bgSeq = null;
		double[] temp = null;
		double sum = 0;
		
		int h = start+disc.getMotifLength( motif )-1, t;
		double w, thresh;
		for( i = 0; i<numSequences; i++ ){
			w = weights==null?1:weights[i];
			bgSeq = bg.getElementAt( i );
			t = bgSeq.getLength() - h;
			num +=  t;
			sum += t*w;
		}
		
		thresh = sum * sign;
		num = (int)Math.ceil( sign * num );
		sortedScores = new double[2][num];

		sortedScores[0][0] = Double.NEGATIVE_INFINITY;
		Arrays.fill( sortedScores[1], 0 );
		
		int stop = 1;
		double partSum = 0;
		for( i = 0 ; i<numSequences; i++ ){
			bgSeq = bg.getElementAt( i );
			temp = getJointProfileOfScoresFor( motif, bgSeq, start );
			Arrays.sort(temp);
			w = weights==null?1:weights[i];
			
//XXX gain performance for first insertions?
			
			//add in sorted list
			for( t = temp.length-1; t >= 0; t-- ) {
				//find position
				int idx = Arrays.binarySearch( sortedScores[0], 0, stop, temp[t] );
				boolean match = idx >= 0; 
				//compute insertion point => sortedScores[0][index-1] < temp[t] <= sortedScores[0][index]
				int index = !match ? (-idx - 1) : idx;
				//System.out.println(stop + "\t" + temp[t] + "\t" + sortedScores[0][0] + "\t" + sortedScores[0][stop-1] + "\t" + sum + "\t" + all );
				if( index == 0 ) {
					//below the significance threshold
					break;
				} else if( index == 1 && partSum > thresh && !match ) {
					//new the significance threshold
					
					//sortedScores[0][0] < temp[t] < sortedScores[0][1]
					sortedScores[0][0] = temp[t];
					break;
				} else {
					partSum+=w;
					//determine how elements that can be integrated into sortedScores[][0]
					int j = 1;
					double del = 0;
					while( j < index && partSum - (del + sortedScores[1][j]) > thresh ) {
						del += sortedScores[1][j];
						j++;
					}
					if( j > 1 ) {
						sortedScores[0][0] = sortedScores[0][j-1];
						sortedScores[1][0] += del;
						partSum -= del;
					}
					
					//keep everything from j to stop
					
					//expand if: stop-(j-1)+(match?0:1) > sortedScores[0].length
					if( stop == sortedScores[0].length && !match && j==1 ) {
						double[][] help = new double[2][2*sortedScores[0].length];
						System.arraycopy( sortedScores[0], 0, help[0], 0, stop );
						System.arraycopy( sortedScores[1], 0, help[1], 0, stop );
						sortedScores = help;
					}
					
					//shift lower part
					if( j != 1 ) {
						System.arraycopy( sortedScores[0], j, sortedScores[0], 1, index-1-j+1 );
						System.arraycopy( sortedScores[1], j, sortedScores[1], 1, index-1-j+1 );
					}
					
					//shift upper part
					int s= index-j+1 + (match?0:1);
					if( s != index ) {
						System.arraycopy( sortedScores[0], index, sortedScores[0], s, stop-index );
						System.arraycopy( sortedScores[1], index, sortedScores[1], s, stop-index );
					}
					
					//set/insert
					if( match ) {
						sortedScores[1][s] += w;	
					} else {
						s--;
						sortedScores[0][s] = temp[t];
						sortedScores[1][s] = w;
					}
					stop = stop-(j-1)+(match?0:1);
				}
			}
			if( t >= 0 ) {
				//above the significance threshold
				sortedScores[1][0] += w*(t+1);
			}
		}
		
		//cut down;
		if( stop < sortedScores[0].length ); {
			double[][] help = new double[2][stop];
			System.arraycopy( sortedScores[0], 0, help[0], 0, stop );
			System.arraycopy( sortedScores[1], 0, help[1], 0, stop );
			sortedScores = help;
		}
		
		
		sum = 0;
		for( i = sortedScores[0].length-1; i>=0; i-- ){
			sum += sortedScores[1][i];
			sortedScores[1][i] = sum;			
		}
		
		for( i = sortedScores[0].length-1; i>0; i-- ){
			sortedScores[1][i] = sortedScores[1][i]/sum;		
		}
		sortedScores[1][0] = Double.NaN;
//XXX bis hier
	}
	
	private void findSignificantMotifOccurrences(int motif, Sequence seq, int start, AbstractList<MotifAnnotation> annotation, int addMax, AbstractList<Sequence> sites, int addLeftSymbols, int addRightSymbols ) throws Exception {
		if( !oneHistogram ) {
			createBgDataSet( seq );
			fillSortedScoresArray( motif, start );
		}
		 
		//System.out.println( scores[0] + " .. " + scores[scores.length-1] );
		
		double[][] allProfs = getAllProfilesOfScoresFor( motif, seq, start );
		double[] joined = joinMethod.joinProfiles( allProfs );

		double thresh = sortedScores[0][1], pVal;
		//XXX Problem
		/*
		System.out.println( sign + "\t" + signIndex + "\t" + thresh + "\t" + sortedScores[1][signIndex]);
		System.out.println( Arrays.toString(joined) );
		*/
		int length = disc.getMotifLength( motif );
		int annotIndex = annotation != null ? annotation.size() : 0;
		int siteIndex = sites != null ? sites.size() : 0;
		DoubleList pValues = new DoubleList();
		double[] probs = new double[2];
		double[] tempStrandProbs = null;
		double[] tempCompProbs = new double[allProfs.length];
		for(int j=0;j<joined.length;j++){
			if(joined[j] > thresh){
				try{
					
					for(int c=0;c<tempCompProbs.length;c++) {
						tempCompProbs[c] = Math.exp(allProfs[c][j] - joined[j]);
					}
					probs[0] = probs[1] = 0;
					for(int i=0;i<tempCompProbs.length;i++){
						tempStrandProbs = disc.getStrandProbabilitiesFor( globalToLocalIndexes[motif][0][i], globalToLocalIndexes[motif][1][i], seq, start+j );
						for(int s=0;s<probs.length;s++){
							probs[s] += tempStrandProbs[s]*tempCompProbs[i];
						}
					}
					
					if( sites != null ) {
						Sequence site;
						
						if( probs[1] > probs[0] ) {
							site = seq.getSubSequence( j+start-addRightSymbols, length + addLeftSymbols + addRightSymbols ).reverseComplement();
						} else {
							site =  seq.getSubSequence( j+start-addLeftSymbols, length + addLeftSymbols + addRightSymbols );
						}
						sites.add( site );
					}
					pVal = getPValue(sortedScores[0],joined[j],sortedScores[1]);
					pValues.add( pVal );
					if( annotation != null ) {
						annotation.add( new MotifAnnotation( "motif* " + motif, j+start, length,
								probs[1] > probs[0] ? Strand.REVERSE : Strand.FORWARD,
								new NumericalResult( "component", "the component of the model where this motif was found with the highest probability", globalToLocalIndexes[motif][0][ getIndexOfMax( tempCompProbs ) ] ),
								new NumericalResult( "p-value", "", pVal ),
								new NumericalResult( "score", "", joined[j] ),
								new NumericalResult( "forward probability", "probability of the forward strand", probs[0] )) );
					}
				} catch( Exception ex ) {
					//TODO
					ex.printStackTrace();
				}
			}

		}
		
		if( pValues.length() > addMax ) {
			//reduce prediction
			double[] array = pValues.toArray();
			Arrays.sort( array );
			
			//System.out.println( Arrays.toString( array ) );
			//System.out.println( array[addMax]);
			
			int i = 0;
			while( i < pValues.length() ) {
				if( pValues.get( i ) >= array[addMax] ) {
					if( annotation != null ) {
						annotation.remove( annotIndex );
					}
					if( sites != null ) {
						sites.remove( siteIndex );
					}
				} else {
					annotIndex++;
					siteIndex++;
				}
				i++;
			}
		}
	}
	
	private static double getPValue( double[] sortedScores, double myScore, double[] cumulative ) {
		int idx = PValueComputation.getIndex( sortedScores, myScore, 0 );
		return idx >= cumulative.length ? 0 : cumulative[idx];
	}
	
	private int getLocalIndexOfMotifInComponent(int component, int motif){
		for(int i=0;i<disc.getNumberOfMotifsInComponent( component );i++){
			if(disc.getGlobalIndexOfMotifInComponent( component, i ) == motif){
				return i;
			}
		}
		return -1;
	}
	
	private void prepareIndices() {
		int n = disc.getNumberOfMotifs();
		this.globalToLocalIndexes = new int[n][][];
		for(int i=0;i<n;i++){
			this.globalToLocalIndexes[i] = computeIndices( i );
		}
	}
	
	private int[][] computeIndices( int motif )
	{
		int num = 0;
		for(int i=0;i<disc.getNumberOfComponents();i++){
			int loc = getLocalIndexOfMotifInComponent( i, motif );
			if(loc > -1){
				num++;
			}
		}
		int[][] idxs = new int[2][num];
		num = 0;
		for(int i=0;i<disc.getNumberOfComponents();i++){
			int loc = getLocalIndexOfMotifInComponent( i, motif );
			if(loc > -1){
				idxs[0][num] = i;
				idxs[1][num] = loc;
				num++;
			}
		}
		return idxs;
	}
	
	/**
	 * This method finds the significant motif occurrences in the sequence.
	 * 
	 * @param motif the motif index
	 * @param seq the sequence
	 * @param start the start position
	 * 
	 * @return an array of {@link MotifAnnotation} for the sequence
	 * 
	 * @throws Exception if the background sample could not be created, or some of the scores could not be computed
	 */
	public MotifAnnotation[] findSignificantMotifOccurrences( int motif, Sequence seq, int start ) throws Exception {
		return findSignificantMotifOccurrences( motif, seq, Integer.MAX_VALUE, start );
	}
	
	/**
	 * This method finds the significant motif occurrences in the sequence.
	 * 
	 * @param motif the motif index
	 * @param seq the sequence
	 * @param addMax the number of motif occurrences that can at most be annotated
	 * @param start the start position
	 * 
	 * @return an array of {@link MotifAnnotation} for the sequence
	 * 
	 * @throws Exception if the background sample could not be created, or some of the scores could not be computed
	 */
	public MotifAnnotation[] findSignificantMotifOccurrences( int motif, Sequence seq, int addMax, int start ) throws Exception {
		LinkedList<MotifAnnotation> list = new LinkedList<MotifAnnotation>();
		if(oneHistogram){
			fillSortedScoresArray( motif, start );
		}
		findSignificantMotifOccurrences( motif, seq, start,list, addMax, null, 0, 0 );
		return list.toArray( new MotifAnnotation[0] );
	}	
	
	public double[][] getPWM( int motif, DataSet data, double[] weights, int addLeft, int addRight ) throws Exception {
		ArrayList<MotifAnnotation> list = new ArrayList<MotifAnnotation>();
		if(oneHistogram){
			fillSortedScoresArray( motif, 0 );
		}
		double w = 1;
		ComplementableDiscreteAlphabet abc = null;
		try {
			abc = (ComplementableDiscreteAlphabet) data.getAlphabetContainer().getAlphabetAt(0);
		} catch( Exception e ) {}
		double[][] pwm = new double[addLeft+addRight+disc.getMotifLength(motif)][(int)data.getAlphabetContainer().getAlphabetLengthAt(0)];
		for( int i = 0; i < data.getNumberOfElements(); i++ ) {
			if( weights != null ) {
				w = weights[i];
			}
			Sequence seq = data.getElementAt(i);
			findSignificantMotifOccurrences( motif, seq, 0, list, Integer.MAX_VALUE, null, 0, 0 );
			for( int l, j = 0; j < list.size(); j++ ) {
				MotifAnnotation ma = list.get(j);
				int start = ma.getPosition() - addLeft;
				int end = ma.getEnd() + addRight;
				Strand strand = ma.getStrandedness();
				switch( strand ) { 
					case FORWARD: 
						l = 0;
						while( start < end && start < seq.getLength() ) {
							if( start >= 0 ) {
								pwm[l][seq.discreteVal(start)] += w;
							}
							l++;
							start++;
						}
						break;
					case REVERSE:
						l = pwm.length-1;
						while( start < end && start < seq.getLength() ) {
							if( start >= 0 ) {
								pwm[l][abc.getComplementaryCode( seq.discreteVal(start) )] += w;
							}
							l--;
							start++;
						}
						break;
					default: throw new RuntimeException();
				}
			}
			list.clear();
		}
		for( int i = 0; i < pwm.length; i++ ) {
			Normalisation.sumNormalisation(pwm[i]);
		}
		return pwm;
	}
	
	
	public Pair<double[][], int[][]> getPWMAndPositions( int motif, DataSet data, double[] weights, int addLeft, int addRight ) throws Exception {
		ArrayList<MotifAnnotation> list = new ArrayList<MotifAnnotation>();
		if(oneHistogram){
			fillSortedScoresArray( motif, 0 );
		}
		int[][] positions = new int[data.getNumberOfElements()][];
		double w = 1;
		ComplementableDiscreteAlphabet abc = null;
		try {
			abc = (ComplementableDiscreteAlphabet) data.getAlphabetContainer().getAlphabetAt(0);
		} catch( Exception e ) {}
		int nbs = 0;
		double[][] pwm = new double[addLeft+addRight+disc.getMotifLength(motif)][(int)data.getAlphabetContainer().getAlphabetLengthAt(0)];
		for( int i = 0; i < data.getNumberOfElements(); i++ ) {
			if( weights != null ) {
				w = weights[i];
			}
			Sequence seq = data.getElementAt(i);
			findSignificantMotifOccurrences( motif, seq, 0, list, Integer.MAX_VALUE, null, 0, 0 );
			positions[i] = new int[list.size()];
			for( int l, j = 0; j < list.size(); j++ ) {
				nbs++;
				MotifAnnotation ma = list.get(j);
				positions[i][j] = ma.getPosition();
				int start = ma.getPosition() - addLeft;
				int end = ma.getEnd() + addRight;
				Strand strand = ma.getStrandedness();
				switch( strand ) { 
					case FORWARD: 
						l = 0;
						while( start < end && start < seq.getLength() ) {
							if( start >= 0 ) {
								pwm[l][seq.discreteVal(start)] += w;
							}
							l++;
							start++;
						}
						break;
					case REVERSE:
						l = pwm.length-1;
						while( start < end && start < seq.getLength() ) {
							if( start >= 0 ) {
								pwm[l][abc.getComplementaryCode( seq.discreteVal(start) )] += w;
							}
							l--;
							start++;
						}
						break;
					default: throw new RuntimeException();
				}
			}
			list.clear();
		}
		System.out.println("nbs: "+nbs);
		for( int i = 0; i < pwm.length; i++ ) {
			Normalisation.sumNormalisation(pwm[i]);
		}
		return new Pair<double[][],int[][]>(pwm,positions);
	}
	
	public Pair<double[][], double[]> getPWMAndPosDist( int motif, DataSet data, double[] weights, double[] mean, int addLeft, int addRight ) throws Exception {
		ArrayList<MotifAnnotation> list = new ArrayList<MotifAnnotation>();
		if(oneHistogram){
			fillSortedScoresArray( motif, 0 );
		}
		double sd = 0, n = 0;
		int nbs = 0;
		double w = 1;
		ComplementableDiscreteAlphabet abc = null;
		try {
			abc = (ComplementableDiscreteAlphabet) data.getAlphabetContainer().getAlphabetAt(0);
		} catch( Exception e ) {}
		double[][] pwm = new double[addLeft+addRight+disc.getMotifLength(motif)][(int)data.getAlphabetContainer().getAlphabetLengthAt(0)];
		for( int i = 0; i < data.getNumberOfElements(); i++ ) {
			if( weights != null ) {
				w = weights[i];
			}
			Sequence seq = data.getElementAt(i);
			findSignificantMotifOccurrences( motif, seq, 0, list, Integer.MAX_VALUE, null, 0, 0 );
			w /= list.size();
			for( int l, j = 0; j < list.size(); j++ ) {
				MotifAnnotation ma = list.get(j);
				int start = ma.getPosition() - addLeft;
				sd += w*(start-mean[i])*(start-mean[i]);
				n += w;
				nbs++;
				int end = ma.getEnd() + addRight;
				Strand strand = ma.getStrandedness();
				switch( strand ) { 
					case FORWARD: 
						l = 0;
						while( start < end && start < seq.getLength() ) {
							if( start >= 0 ) {
								pwm[l][seq.discreteVal(start)] += w;
							}
							l++;
							start++;
						}
						break;
					case REVERSE:
						l = pwm.length-1;
						while( start < end && start < seq.getLength() ) {
							if( start >= 0 ) {
								pwm[l][abc.getComplementaryCode( seq.discreteVal(start) )] += w;
							}
							l--;
							start++;
						}
						break;
					default: throw new RuntimeException();
				}
			}
			list.clear();
		}
		System.out.println("nbs: "+nbs);
		for( int i = 0; i < pwm.length; i++ ) {
			Normalisation.sumNormalisation(pwm[i]);
		}
		return new Pair<double[][],double[]>(pwm,new double[]{Math.sqrt( sd/n )});
	}
	
	/**
	 * This method annotates a {@link DataSet}.
	 * 
	 * @param data the {@link DataSet}
	 * @param motifIndex the index of the motif
	 * 
	 * @return an annotated {@link DataSet}
	 * 
	 * @throws Exception if something went wrong
	 * 
	 * @see SignificantMotifOccurrencesFinder#annotateMotif(int, DataSet, int)
	 */
	public DataSet annotateMotif( DataSet data, int motifIndex ) throws Exception
	{
		return annotateMotif( 0, data, motifIndex );
	}
	
	/**
	 * This method annotates a {@link DataSet} starting in each sequence at <code>startPos</code>.
	 * 
	 * @param startPos the start position used for all sequences
	 * @param data the {@link DataSet}
	 * @param motifIndex the index of the motif
	 * 
	 * @return an annotated {@link DataSet}
	 * 
	 * @throws Exception if something went wrong
	 * 
	 * @see SignificantMotifOccurrencesFinder#annotateMotif(int, DataSet, int)
	 */
	public DataSet annotateMotif( int startPos, DataSet data, int motifIndex ) throws Exception
	{
		return annotateMotif( startPos, data, motifIndex, Integer.MAX_VALUE, false );
	}
	
	/**
	 * This method annotates a {@link DataSet}.
	 * At most, <code>addMax</code> motif occurrences of the motif instance will be annotated.
	 * 
	 * @param data the {@link DataSet}
	 * @param motifIndex the index of the motif
	 * @param addMax the number of motif occurrences that can at most be annotated for each motif instance
	 * 
	 * @return an annotated {@link DataSet}
	 * 
	 * @throws Exception if something went wrong
	 * 
	 * @see SignificantMotifOccurrencesFinder#annotateMotif(int, DataSet, int)
	 */
	public DataSet annotateMotif( DataSet data, int motifIndex, int addMax ) throws Exception
	{
		return annotateMotif( 0, data, motifIndex, addMax, false );
	}
	
	/**
	 * This method annotates a {@link DataSet} starting in each sequence at <code>startPos</code>.
	 * At most, <code>addMax</code> motif occurrences of the motif instance will be annotated.
	 * 
	 * @param startPos the start position used for all sequences
	 * @param data the {@link DataSet}
	 * @param motifIndex the index of the motif
	 * @param addMax the number of motif occurrences that can at most be annotated for each motif instance
	 * @param addAnnotation a switch whether to add or replace the current annotation
	 * 
	 * @return an annotated {@link DataSet}
	 * 
	 * @throws Exception if something went wrong
	 * 
	 * @see SignificantMotifOccurrencesFinder#annotateMotif(int, DataSet, int)
	 */
	public DataSet annotateMotif( int startPos, DataSet data, int motifIndex, int addMax, boolean addAnnotation ) throws Exception {
		return (DataSet) predictBS( startPos, data, null, motifIndex, addMax, 0, 0, addAnnotation ).get( 0 );
	}
	
	/**
	 * This method returns a {@link DataSet} containing the predicted binding sites.
	 * 
	 * @param data the {@link DataSet}
	 * @param motifIndex the index of the motif
	 * 
	 * @return a {@link DataSet} containing the predicted binding sites
	 * 
	 * @throws Exception if something went wrong
	 */
	public DataSet getBindingSites( DataSet data, int motifIndex ) throws Exception
	{
		return getBindingSites( 0, data, motifIndex, Integer.MAX_VALUE, 0, 0 );
	}
	
	/**
	 * This method returns a {@link DataSet} containing the predicted binding sites.
	 * 
	 * @param startPos the start position used for all sequences
	 * @param data the {@link DataSet}
	 * @param motifIndex the index of the motif
	 * @param addMax the number of motif occurrences that can at most be annotated for each motif instance
	 * @param addLeft number of positions added to the left of the predicted motif occurrence
	 * @param addRight number of positions added to the right of the predicted motif occurrence
	 * @return a {@link DataSet} containing the predicted binding sites
	 * 
	 * @throws Exception if something went wrong
	 */
	public DataSet getBindingSites( int startPos, DataSet data, int motifIndex, int addMax, int addLeft, int addRight ) throws Exception
	{
		return (DataSet) predictBS( startPos, data, null, motifIndex, addMax, addLeft, addRight, false ).get( 1 );
	}
	
	/**
	 * This method returns a list of start positions of binding sites.
	 * 
	 * @param startPos the start position used for all sequences
	 * @param data the {@link DataSet}
	 * @param motifIndex the index of the motif
	 * @param addMax the number of motif occurrences that can at most be annotated for each motif instance
	 * 
	 * @return a list of start positions
	 * 
	 * @throws Exception if something went wrong
	 */
	public IntList getStartPositions( int startPos, DataSet data, int motifIndex, int addMax ) throws Exception
	{
		return (IntList) predictBS( startPos, data, null, motifIndex, addMax, 0, 0, false ).get( 3 );
	}
	
	/**
	 * Returns the number of sequences in <code>data</code> that are predicted to be bound at least once by motif no. <code>motif</code>.
	 * 
	 * @param data the data
	 * @param weights the weights of the data
	 * @param motifIndex the index of the motif
	 * 
	 * @return the number of sequences in <code>data</code> bound by motif <code>motif</code>
	 * 
	 * @throws Exception if the background sample for the prediction could not be created or some of the scores could not be computed
	 */
	public double getNumberOfBoundSequences( DataSet data, double[] weights, int motifIndex ) throws Exception
	{
		return (Double) predictBS( 0, data, weights, motifIndex, Integer.MAX_VALUE, 0, 0, false ).get(2);
	}
	
	/**
	 * This method returns an offset that must be added to scores for computing PR curves. If this {@link SignificantMotifOccurrencesFinder} 
	 * was instantiated using <code>oneHistogram=true</code>, the {@link SignificantMotifOccurrencesFinder#getValuesForEachNucleotide(DataSet, int, boolean)} returns scores and no offset is needed. Otherwise,
	 * it returns p-values and, hence, 1-(p-value) must be used for the PR curve and the offset is 1.
	 * 
	 * @return the offset
	 * 
	 * @see MotifDiscoveryAssessment#getSortedValuesForMotifAndFlanking(DataSet, double[][], double, double, String)
	 */
	public double getOffsetForAucPR() {
		return oneHistogram ? 0 : 1;
	}

	/**
	 * This method returns a factor that must be multiplied to scores for computing PR curves. If this {@link SignificantMotifOccurrencesFinder} 
	 * was instantiated using <code>oneHistogram=true</code>, the {@link SignificantMotifOccurrencesFinder#getValuesForEachNucleotide(DataSet, int, boolean)} returns scores and a factor of 1 is appropriate. Otherwise,
	 * it returns p-values and, hence, 1-(p-value) must be used for the PR curve and the factor is -1.
	 * 
	 * @return the factor
	 * 
	 * @see MotifDiscoveryAssessment#getSortedValuesForMotifAndFlanking(DataSet, double[][], double, double, String)
	 */
	public double getFactorForAucPR() {
		return oneHistogram ? 1 : -1;
	}
	
	/**
	 * This method determines a score for each possible starting position in each of the sequences in <code>data</code> 
	 * that this position is covered by at least one motif occurrence of the
	 * motif with index <code>index</code>. If the {@link SignificantMotifOccurrencesFinder}
	 * was constructed using <code>oneHistogram=true</code> the returned values are arbitrary scores, and p-values otherwise.
	 * 
	 * @param data the {@link DataSet}
	 * @param motif the motif index
	 * @param addOnlyBest a switch whether to add only the best
	 * 
	 * @return an array containing for each sequence an array with the scores for each starting position in the sequence
	 * 
	 * @throws Exception if something went wrong during the computation of the scores of the {@link MotifDiscoverer} 
	 * 
	 * @see MotifDiscoveryAssessment#getSortedValuesForMotifAndFlanking(DataSet, double[][], double, double, String)
	 * @see #getOffsetForAucPR()
	 * @see #getFactorForAucPR()
	 */
	public double[][] getValuesForEachNucleotide( DataSet data, int motif, boolean addOnlyBest ) throws Exception {
		double[][] res = new double[data.getNumberOfElements()][];
		if( oneHistogram ) {
			fillSortedScoresArray( motif, 0 );
		}
		for( int i = 0; i < res.length; i++ ) {
			res[i] = getValueForNucleotides(data.getElementAt(i), 0, motif, addOnlyBest);
		}
		return res;
	}
	
	
	private static int getIndexOfMax( double... values ) {
		int idx = 0, i = 1;
		for( ; i < values.length; i++ ) {
			if( values[i] > values[idx] ) {
				idx = i;
			}
		}
		return idx;
	}
	
	private double[] getValueForNucleotides(Sequence seq, int start, int motif, boolean addOnlyBest ) throws Exception{
		if( !oneHistogram ) {
			createBgDataSet( seq );
			fillSortedScoresArray( motif, start );
		}
		
		double[] temp = getJointProfileOfScoresFor( motif, seq, start );
	
		//naive approach
		double[] res;
		int i, length = disc.getMotifLength( motif );
		
		if( addOnlyBest ) {
			res = new double[seq.getLength()-start];
			Arrays.fill( res, oneHistogram? Double.NEGATIVE_INFINITY : 1 );
			int idx = getIndexOfMax( temp );
			double best = oneHistogram ? temp[idx] : getPValue(sortedScores[0],temp[idx],sortedScores[1]);
			for( i = 0; i < length; i++ ){
				res[idx+i] = best;
			}
		} else {
			res = smooth( temp, seq.getLength()-start, length );
			if( !oneHistogram ) {
				for(i=0;i<temp.length;i++){
					res[i] = getPValue(sortedScores[0],temp[i],sortedScores[1]);
				}
			}
		}
		return res;
	}
	
	private static double[] smooth( double[] temp, int seqLength, int motifLength ) {
		double[] res = new double[seqLength];
		Arrays.fill( res, temp[temp.length-1] );
		System.arraycopy( temp, 0, res, 0, temp.length );
		
		for(int k,j,i=res.length-1; i >= 0; i--) {
			for(k=i-1,j=1; j < motifLength && k >= 0; j++, k--){
				if( res[i] < res[k] ) {
					res[i] = res[k];
				}
			}
		}
		return res;
	}
	
	private ArrayList predictBS( int startPos, DataSet data, double[] weights, int motif, int addMax, int addLeft, int addRight, boolean addAnnotation ) throws Exception
	{
		int i, n = data.getNumberOfElements();
		Sequence[] seqs = new Sequence[n];
		IntList posList = new IntList();
		LinkedList<MotifAnnotation> seqAn = new LinkedList<MotifAnnotation>();
		LinkedList<Sequence> bsList = new LinkedList<Sequence>();
		LinkedList<Sequence> currentList = new LinkedList<Sequence>();
		SequenceAnnotation[] empty = new SequenceAnnotation[0];
		MotifAnnotation current;

		if( oneHistogram ) {
			createBgDataSet( data );
			fillSortedScoresArray( motif, startPos );
		}
		
		double weight = 1,  bound = 0;
		for( i = 0; i < n; i++ )
		{
			seqs[i] = data.getElementAt(i);
			if( weights != null ) {
				weight = weights[i];
			}
			//collect annotation in seqAn
			seqAn.clear();
			findSignificantMotifOccurrences( motif, seqs[i], startPos, seqAn, addMax, currentList, addLeft, addRight );
			if( currentList.size() > 0 ) {
				bound += weight;
				bsList.addAll( currentList );
				currentList.clear();
			}
			
			//replace annotation with those currently computed
			seqs[i] = seqs[i].annotate( addAnnotation, seqAn.toArray( empty ) );
			for( int k = 0; k < seqAn.size(); k++ )
			{
				current = seqAn.get( k );
				posList.add( current.getPosition() );
			}
		}
		ArrayList res = new ArrayList(4);
		res.add( new DataSet( "annotated sample", seqs ) );
		DataSet bs;
		try {
			bs = new DataSet( "annotated binding sites", bsList.toArray( new Sequence[0] ) );
		} catch( Exception e ) {
			bs = null;
		}
		res.add( bs );
		res.add( bound );
		res.add( posList );
		return res;
	}
	
	/**
	 * This method returns a clone of the internally used {@link MotifDiscoverer}.
	 * 
	 * @return clone of the internally used {@link MotifDiscoverer}
	 * @throws CloneNotSupportedException if the {@link MotifDiscoverer} can not be cloned correctly
	 */
	public MotifDiscoverer getMotifDiscoverer() throws CloneNotSupportedException {
		return disc.clone();
	}
	
	/**
	 * Interface for methods that combine several profiles over the same sequence
	 * into one common profile
	 * @author Jan Grau
	 *
	 */
	public static interface JoinMethod {
		
		/**
		 * Joins the profiles in <code>profiles</code> into a common
		 * profile, where each line of <code>profiles</code> corresponds to one profile
		 * to be joined
		 * @param profiles the profiles to be joined
		 * @return the joined profiles
		 * @throws Exception if the profiles could not be joined, e.g. because they are not of equal length
		 */
		public double[] joinProfiles(double[][] profiles) throws Exception;
		
	}
	
	/**
	 * Joins several profiles containing log-probabilities into one profile containing
	 * the logarithm of the sum of the probabilities of the single profiles.
	 * 
	 * @author Jan Grau
	 *
	 */
	public static class SumOfProbabilities implements JoinMethod {
		
		@Override
		public double[] joinProfiles(double[][] profiles){
			double[] res = new double[profiles[0].length];
			for(int i=0;i<profiles.length;i++){
				if(profiles[i].length != res.length){
					throw new IllegalArgumentException( "Profiles must be of same length, but profile "+i+" has length "+profiles[i].length+" instead of "+res.length );
				}
			}
			double[] temp = new double[profiles.length];
			for(int i=0;i<res.length;i++){
				for(int j=0;j<temp.length;j++){
					temp[j] = profiles[j][i];
				}
				res[i] = Normalisation.getLogSum( temp );
			}

			return res;
		}
	}
}
