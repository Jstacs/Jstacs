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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.Map.Entry;

import javax.naming.OperationNotSupportedException;

import de.jstacs.WrongAlphabetException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.DataSet.WeightedDataSetFactory;
import de.jstacs.data.DataSet.WeightedDataSetFactory.SortOperation;
import de.jstacs.io.ArrayHandler;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Pair;

/**
 * This class enables the user to get some statistics of a {@link DataSet} in an easy way.
 * 
 * @author Jens Keilwagen
 */
public final class KMereStatistic {
	
	private AlphabetContainer abc;
	private int length, k, n;
	private int[] powers;
	private int[][] counts;
	
	/**
	 * This constructor creates an internal statistic counting all <code>k</code>-mers in the <code>data</code>.
	 * 
	 * @param data the data
	 * @param k the number of symbols in each counted word
	 */
	public KMereStatistic( DataSet data, int k ) {
		abc = data.getAlphabetContainer();
		length = data.getElementLength();
		if( !abc.isSimple() || length == 0 ) {
			throw new IllegalArgumentException( "Can not compute the statistic: check the sample" );
		}
		if( k < 1 || k >= length ) {
			throw new IllegalArgumentException( "Can not compute the statistic: check the order" );
		}
		this.k = k;
		powers = new int[Math.max( k + 1, 2 )];
		powers[0] = 1;
		powers[1] = (int) abc.getAlphabetLengthAt( 0 );
		int i = 2, l, idx, s;
		for( ; i < powers.length; i++ )
		{
			powers[i] = powers[1] * powers[i - 1];
		}
		
		counts = new int[length-k+1][powers[k]];
		Sequence seq;
		for( i = 0; i< data.getNumberOfElements(); i++ ) {
			seq = data.getElementAt( i );
			idx = 0;
			for( l = 0; l < k-1; l++ )
			{
				idx = idx * powers[1] + seq.discreteVal( l );
			}
			s = 0;
			for( ; l < length; l++, s++ )
			{
				idx = ( (idx * powers[1])  + seq.discreteVal( l ) ) % powers[k];
				counts[s][idx]++;
			}
		}
		n = data.getNumberOfElements();
	}
	
	/**
	 * This method returns an array of smoothed profiles. For each k-mere it returns one profile.
	 * The order of the profile is the same as the order of the k-meres.
	 * 
	 * @param window the window length, for no smoothing use 1 
	 * @param kmere the k-mere
	 * 
	 * @return an array of smoothed profiles
	 * 
	 * @see #getSmoothedProfile(int, Sequence...)
	 * @see Sequence#create(AlphabetContainer, String)
	 */
	public double[][] getSmoothedProfile( int window, String... kmere )
	{
		Sequence[] seq = null;
		try {
			seq = new Sequence[kmere.length];
			for( int i = 0; i < seq.length; i++ ) {
				seq[i] = Sequence.create( abc, kmere[i] );
			}
		} catch( Exception e ) {
			seq = null;
		}
		if( seq == null ) {
			throw new IllegalArgumentException();
		}
		return getSmoothedProfile( window, seq );
	}
	
	/**
	 * This method returns an array of smoothed profiles. For each k-mere it returns one profile.
	 * The order of the profile is the same as the order of the k-meres.
	 * 
	 * @param window the window length, for no smoothing use 1 
	 * @param seq the {@link Sequence} instances containing the k-meres
	 * 
	 * @return an array of smoothed profiles
	 */
	public double[][] getSmoothedProfile( int window, Sequence... seq )
	{
		if( window < 1 )
		{
			throw new IllegalArgumentException( "The window has to have at least length 1." );
		}
		
		if( seq == null ) {
			throw new IllegalArgumentException( "check the subsequences" );
		} else {
			int i = 0, l, idx, s, m = n*window;
			double[][] res = new double[seq.length][counts.length-window+1];
			for( i = 0; i < seq.length; i++ ) {
				//check
				if( seq[i].getLength() != k ) {
					throw new IllegalArgumentException();
				}
				
				//get idx;
				idx = 0;
				for( l = 0; l < k; l++ )
				{
					idx = idx * powers[1] + seq[i].discreteVal( l );
				}
				
				//sum;
				for( l = 0; l < window; l++ )
				{
					res[i][0] += counts[l][idx];
				}
				s = 1;
				for( ; s < res[i].length; l++, s++ )
				{
					res[i][s] = res[i][s-1] - counts[l-window][idx] + counts[l][idx];
					res[i][s-1] /= m;
				}
				res[i][s-1] /= m;
			}
			return res;
		}		
	}

	/**
	 * This method returns an array of strings of length
	 * <code>motifLength</code> so that each String is contained in all
	 * sequences of the sample respectively in the sample and the reverse
	 * complementary sample.
	 * 
	 * @param data
	 *            the sample of sequences
	 * @param motifLength
	 *            the motif length
	 * @param bothStrands
	 *            the switch for using both strand <code>true</code> or only
	 *            forward strand <code>false</code>
	 * 
	 * @return an array of Strings of length <code>motifLength</code> so that
	 *         each String is contained in <code>data</code> respectively on
	 *         one strand of the <code>data</code>
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static Sequence[] getCommonString(DataSet data, int motifLength, boolean bothStrands) throws Exception {
		int f = bothStrands ? 2 : 1;
		LinkedList<Sequence> candidates = new LinkedList<Sequence>();
		HashSet<Sequence> current = new HashSet<Sequence>(f * data.getMaximalElementLength());
		Sequence s = data.getElementAt(0);
		int i = 0, end = s.getLength() - motifLength;
		Sequence help;
		while( i <= end ){
			help = s.getSubSequence( i, motifLength );
			if( !current.contains( help ) ) {
				current.add( help );
				candidates.add( help );
			}
			i++;
		}
		for( int j = 1; candidates.size() > 0 && j < data.getNumberOfElements(); j++ ) {
			current.clear();
			add(current, data.getElementAt(j), motifLength);
			intersection(candidates, current, bothStrands);
		}
		Sequence[] empty = new Sequence[0];
		empty = candidates.toArray(empty);
		Arrays.sort(empty);
		return empty;
	}

	private static void add(HashSet<Sequence> hash, Sequence sequence, int motifLength) {
		int i = 0, end = sequence.getLength() - motifLength;
		Sequence help;
		while (i <= end) {
			help = sequence.getSubSequence( i, motifLength );
			if( !hash.contains( help ) ) {
				hash.add(help);
			}
			i++;
		}
	}

	private static void intersection(LinkedList<Sequence> candidates, HashSet<Sequence> current, boolean bothStrands ) throws OperationNotSupportedException {
		int i = 0;
		while( i < candidates.size() ) {
			Sequence seq = candidates.get( i );
			if ( !( current.contains( seq ) || current.contains( seq.reverseComplement() ) ) ) {
				candidates.remove( i );
			} else {
				i++;
			}
		}
	}

	/**
	 * This method enables the user to get a statistic over all <code>k</code>-mers
	 * in the <code>data</code>. That is it counts the outcome of each
	 * <code>k</code>-mere in the complete <code>data</code>.
	 * 
	 * @param data
	 *            the sample of sequences
	 * @param k
	 *            the motif length
	 * @param bothStrands
	 *            the switch for using both strand <code>true</code> or only
	 *            forward strand <code>false</code>. If <code>true</code>
	 *            for each <code>k</code>-mer only this <code>k</code>-mere
	 *            or its reverse complement is contained in the returned
	 *            WeightedSampleFactory.
	 * 
	 * @return a WeightedSampleFactory containing all <code>k</code>-mers and
	 *         their absolute frequencies in <code>data</code> respectively on
	 *         one strand of the <code>data</code>
	 * 
	 * @throws Exception
	 *             if something went wrong
	 *             
	 * @see KMereStatistic#getAbsoluteKMereFrequencies(DataSet, int, boolean, DataSet.WeightedDataSetFactory.SortOperation)
	 * @see SortOperation#NO_SORT
	 */
	public static WeightedDataSetFactory getAbsoluteKMereFrequencies( DataSet data, int k, boolean bothStrands) throws Exception {
		return getAbsoluteKMereFrequencies( data, k, bothStrands, SortOperation.NO_SORT );
	}
	
	/**
	 * This method enables the user to get a statistic over all <code>k</code>-mers
	 * in the <code>data</code>. That is it counts the outcome of each
	 * <code>k</code>-mere in the complete <code>data</code>.
	 * 
	 * @param data
	 *            the sample of sequences
	 * @param k
	 *            the motif length
	 * @param bothStrands
	 *            the switch for using both strand <code>true</code> or only
	 *            forward strand <code>false</code>. If <code>true</code>
	 *            for each <code>k</code>-mer only this <code>k</code>-mere
	 *            or its reverse complement is contained in the returned
	 *            WeightedSampleFactory.
	 * @param sortOp
	 * 			  the way how the result should be sorted
	 * 
	 * @return a WeightedSampleFactory containing all <code>k</code>-mers and
	 *         their absolute frequencies in <code>data</code> respectively on
	 *         one strand of the <code>data</code>
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static WeightedDataSetFactory getAbsoluteKMereFrequencies( DataSet data, int k, boolean bothStrands, SortOperation sortOp ) throws Exception {
		DataSet myData = data;
		if (bothStrands) {
			Sequence[] seqs = new Sequence[2 * data.getNumberOfElements()];
			for (int j = 0; j < data.getNumberOfElements(); j++) {
				seqs[2 * j] = data.getElementAt(j);
				seqs[2 * j + 1] = seqs[2 * j].reverseComplement();
			}
			myData = new DataSet("both strands of " + data.getAnnotation(), seqs);
		}
		WeightedDataSetFactory wsf = new WeightedDataSetFactory( sortOp, myData, null, k);
		if (bothStrands) {
			wsf = removeReverseComplements(wsf, 2, sortOp );
		}
		return wsf;
	}
	
	/**
	 * This method enables the user to get a statistic over all <code>k</code>-mers
	 * in the sequences. That is, it creates for each occurring <code>k</code>-mer an array
	 * of {@link BitSet}s indicating for each data set and each sequence whether it contains
	 * the <code>k</code>-mer (or its reverse complement) or not.
	 * 
	 * @param data
	 *            the {@link DataSet}s of {@link Sequence}s
	 * @param k
	 *            the motif length
	 * @param bothStrands
	 *            the switch for using both strand <code>true</code> or only
	 *            forward strand <code>false</code>. If <code>true</code>
	 *            for each <code>k</code>-mer only this <code>k</code>-mere
	 *            or its reverse complement is contained in the returned
	 *            WeightedSampleFactory.
	 * @param addIndex
	 *            the maximal index for inserting new k-meres 
	 * 
	 * @return a {@link Hashtable} on {@link Sequence}s and arrays of {@link BitSet}s; each
	 * 		   entry encodes a <code>k</code>-mer and the occurrence of this <code>k</code>-mer
	 *         in each data set and sequence; if a <code>k</code>-mer occurs in data set
	 *         <code>d</code> in sequence <code>n</code> the <code>n</code>-th bit of the
	 *         <code>d</code>-th {@link BitSet} is true.
	 *          
	 * @throws WrongAlphabetException if the {@link AlphabetContainer}s of the {@link DataSet}s do not match or if they are not simple and discrete 
	 * @throws OperationNotSupportedException if the <code>bothStrands==true</code> but the reverse complement could not be computed
	 * 
	 * @see Hashtable
	 * @see KMereStatistic#merge(Hashtable, int, boolean)
	 */
	public static Hashtable<Sequence, BitSet[]> getKmereSequenceStatistic( int k, boolean bothStrands, int addIndex, DataSet... data ) throws WrongAlphabetException, OperationNotSupportedException {
		AlphabetContainer con = data[0].getAlphabetContainer();
		if( !con.isSimple() || !con.isDiscrete() ) {
			throw new WrongAlphabetException();
		}
		int[] anz = new int[data.length];
		for( int d = 0; d < data.length; d++ ) {
			if( !con.checkConsistency( data[d].getAlphabetContainer() ) ) {
				throw new WrongAlphabetException();
			}
			anz[d] = data[d].getNumberOfElements();
		}
		Hashtable<Sequence, BitSet[]> res = new Hashtable<Sequence, BitSet[]>();
		Sequence seq, current;
		BitSet[] b;
		boolean add = false;
		//run over all sequences
		for( int h, m, l, n, d = 0; d < data.length; d++ ) {
			for( n = 0; n < anz[d]; n++ ) {
				seq = data[d].getElementAt( n );
				m = seq.getLength()-k+1;
				//run over all k-mers
				for( l = 0; l < m; l++ ) {
					current = seq.getSubSequence( con, l, k );
					b = res.get( current );
					if( b == null && bothStrands ) {
						b = res.get( current.reverseComplement() );
					}
					if( b != null || d <= addIndex ) {
						if( b == null ) {
							b = createBitSets( anz );
							add = true;
						}
						b[d].set( n );
						if( add ) {
							res.put( current, b );
							add = false;
						}
					}
				}
			}
		}
		return res;
	}
	
	private static BitSet[] createBitSets( int[] anz ) {
		BitSet[] b = new BitSet[anz.length];
		for( int h = 0; h < anz.length; h++ ) {
			b[h] = new BitSet( anz[h] );
		}
		return b;
	}
	
	/**
	 * This method enables the user to get a statistic for a set of <code>k</code>-mers.
	 * That is, it creates for each <code>k</code>-mer from <code>filter</code> an array
	 * of {@link BitSet}s indicating for each data set and each sequence whether it contains
	 * the <code>k</code>-mer (or its reverse complement) or not.
	 * 
	 * @param bothStrands
	 *            the switch for using both strand <code>true</code> or only
	 *            forward strand <code>false</code>. If <code>true</code>
	 *            for each <code>k</code>-mer only this <code>k</code>-mere
	 *            or its reverse complement is contained in the returned
	 *            WeightedSampleFactory.
	 * @param maxMismatch
	 *            the maximal number of mismatches
	 * @param filter
	 *            a filter containing all interesting <code>k</code>-mers
	 * @param data
	 *            the {@link DataSet}s of {@link Sequence}s
	 * 
	 * @return a {@link Hashtable} on {@link Sequence}s and arrays of {@link BitSet}s; each
	 * 		   entry encodes a <code>k</code>-mer and the occurrence of this <code>k</code>-mer
	 *         in each data set and sequence; if a <code>k</code>-mer occurs in data set
	 *         <code>d</code> in sequence <code>n</code> the <code>n</code>-th bit of the
	 *         <code>d</code>-th {@link BitSet} is true.
	 *         
	 * @throws WrongAlphabetException if the {@link AlphabetContainer}s of the {@link DataSet}s do not match or if they are not simple and discrete 
	 * @throws OperationNotSupportedException if the <code>bothStrands==true</code> but the reverse complement could not be computed
	 * 
	 * @see Hashtable
	 * @see KMereStatistic#merge(Hashtable, int, boolean)
	 */
	public static Pair<Sequence, BitSet[]>[] getKmereSequenceStatistic( boolean bothStrands, int maxMismatch, HashSet<Sequence> filter, DataSet... data ) throws WrongAlphabetException, OperationNotSupportedException {
		AlphabetContainer con = data[0].getAlphabetContainer();
		if( !con.isSimple() || !con.isDiscrete() ) {
			throw new WrongAlphabetException();
		}
		int[] anz = new int[data.length];
		for( int d = 0; d < data.length; d++ ) {
			if( !con.checkConsistency( data[d].getAlphabetContainer() ) ) {
				throw new WrongAlphabetException();
			}
			anz[d] = data[d].getNumberOfElements();
		}
		Pair<Sequence, BitSet[]>[] res = new Pair[filter.size()];
		Iterator<Sequence> it = filter.iterator();
		int i = 0;
		while( it.hasNext() ) {
			res[i++] = new Pair( it.next(), createBitSets(anz) );
		}
		Sequence seq, current;
		BitSet[] b;
		//run over all sequences
		Iterator<Pair<Sequence, BitSet[]>> it2;
		for( int n, d = 0; d < data.length; d++ ) {
			for( n = 0; n < anz[d]; n++ ) {
				//System.out.println(d + "\t" + n);
				seq = data[d].getElementAt( n );
				for( i = 0; i < res.length; i++ ) {
					current = res[i].getFirstElement();
					b = res[i].getSecondElement();
					if( contains(seq, current, maxMismatch, bothStrands) ) {
						b[d].set( n );
					}
				}
			}
		}
		return res;
	}
	
	private static boolean contains( Sequence seq, Sequence kmere, int maxMismatch, boolean bothStrands ) throws WrongAlphabetException, OperationNotSupportedException {
		boolean res = seq.matches(maxMismatch, kmere);
		if( !res && bothStrands ) {
			res = seq.matches(maxMismatch, kmere.reverseComplement());
		}
		return res;
	}
	
	/**
	 * This method allows to merge the statistics of k-mers by allowing mismatches.
	 * 
	 * @param statistic a statistic as obtained from {@link KMereStatistic#getKmereSequenceStatistic(int, boolean, int, DataSet...)}
	 * @param maximalMissmatch the maximal number of allowed mismatches 
	 * @param bothStrands the switch for using both strand <code>true</code> or only forward strand <code>false</code>.
	 * 
	 * @return a merged statistic
	 * 
	 * @throws OperationNotSupportedException if the <code>bothStrands==true</code> but the reverse complement could not be computed
	 * @throws CloneNotSupportedException if an array of {@link BitSet} can not be cloned
	 * @throws WrongAlphabetException see {@link Sequence#getHammingDistance(Sequence)}
	 * @throws WrongLengthException see {@link Sequence#getHammingDistance(Sequence)}
	 * 
	 * @see Sequence#getHammingDistance(Sequence)
	 * @see KMereStatistic#getKmereSequenceStatistic(int, boolean, int, DataSet...)
	 */
	public static Hashtable<Sequence, BitSet[]> merge( Hashtable<Sequence, BitSet[]> statistic, int maximalMissmatch, boolean bothStrands ) throws OperationNotSupportedException, CloneNotSupportedException, WrongLengthException, WrongAlphabetException {
		Hashtable<Sequence, BitSet[]> res = new Hashtable<Sequence, BitSet[]>();
		Set<Entry<Sequence, BitSet[]>> set = statistic.entrySet();
		Sequence s1, s2, rc = null;
		BitSet[] b1, b2, o1, o2;
		int idx, d;
		
		Entry<Sequence,BitSet[]>[] array = (Entry<Sequence,BitSet[]>[]) ArrayHandler.cast( set.toArray() );
		for( int i = 0; i < array.length; i++ ) {
			res.put( array[i].getKey(), ArrayHandler.clone( array[i].getValue() ) );
		}
		for( int arrayIndex2, arrayIndex1 = 0; arrayIndex1 < array.length; arrayIndex1++ ) {
			s1 = array[arrayIndex1].getKey();
			o1 = array[arrayIndex1].getValue();
			b1 = res.get( s1 );
			if( bothStrands ) {
				rc = s1.reverseComplement();
			}
			
			for( arrayIndex2 = arrayIndex1+1; arrayIndex2 < array.length; arrayIndex2++ ) {
				s2 = array[arrayIndex2].getKey();
				d = s1.getHammingDistance( s2 );
				if( bothStrands ) {
					d = Math.min( d, rc.getHammingDistance( s2 ) );
				}
				if( d >= 0 && d <= maximalMissmatch ) {
					o2 = array[arrayIndex2].getValue();
					b2 = res.get( s2 );
					for( idx = 0; idx < b1.length; idx++ ) {
						b1[idx].or( o2[idx] );
						b2[idx].or( o1[idx] );
					}
				}
			}
		}
		return res;
	}

	/**
	 * This method returns a list of {@link Sequence}s. Each entry corresponds to a sequence
	 * or a set of sequences (depending on the input of the <code>statistic</code>) that occurs
	 * in more than <code>threshold</code> {@link Sequence}s of the data set.
	 * 
	 * @param statistic a statistic as obtained from {@link KMereStatistic#getKmereSequenceStatistic(int, boolean, int, DataSet...)} or {@link KMereStatistic#merge(Hashtable, int, boolean)}
	 * @param dataSetIndex the index of the {@link BitSet} to be used
	 * @param threshold a threshold that has to be exceeded by {@link BitSet#cardinality()} to be declared as a conserved pattern
	 * 
	 * @return a list of conserved patterns
	 * 
	 * @see KMereStatistic#getKmereSequenceStatistic(int, boolean, int, DataSet...)
	 * @see KMereStatistic#merge(Hashtable, int, boolean)
	 */
	public static LinkedList<Sequence> getConservedPatterns( Hashtable<Sequence, BitSet[]> statistic, int dataSetIndex, int threshold ) {
		Iterator<Entry<Sequence, BitSet[]>> it = statistic.entrySet().iterator();
		Entry<Sequence, BitSet[]> e;
		LinkedList<Sequence> list = new LinkedList<Sequence>();
		while( it.hasNext() ) {
			e = it.next();
			if( e.getValue()[dataSetIndex].cardinality() >= threshold ) {
				list.add( e.getKey() );
			}
		}
		return list;
	}
	
	/**
	 * This method allows to remove those entries from the statistic that have a lower weighted foreground cardinality than the weighted background cardinality.
	 * 
	 * @param statistic a statistic as obtained from {@link KMereStatistic#getKmereSequenceStatistic(int, boolean, int, DataSet...)} or {@link KMereStatistic#merge(Hashtable, int, boolean)}
	 * @param fgIndex the foreground index of the {@link BitSet} to be used
	 * @param bgIndex the background index of the {@link BitSet} to be used
	 * @param fgWeight the weight used to weight the foreground cardinality
	 * @param bgWeight the weight used to weight the background cardinality
	 * 
	 * @return a {@link Hashtable} containing only the positive entries
	 */
	public static Hashtable<Sequence, BitSet[]> removeBackground( Hashtable<Sequence, BitSet[]> statistic, int fgIndex, int bgIndex, double fgWeight, double bgWeight ) {
		Hashtable<Sequence, BitSet[]> res = new Hashtable<Sequence, BitSet[]>();
		Iterator<Entry<Sequence, BitSet[]>> it = statistic.entrySet().iterator();
		Entry<Sequence, BitSet[]> e;
		BitSet[] b;
		while( it.hasNext() ) {
			e = it.next();
			b = e.getValue();
			if( b[fgIndex].cardinality() * fgWeight > b[bgIndex].cardinality() * bgWeight ) {
				res.put( e.getKey(), e.getValue() );
			}
		}
		return res;
	}	

	private static WeightedDataSetFactory removeReverseComplements( WeightedDataSetFactory wsf, int div, SortOperation so ) throws Exception {
		ArrayList<Sequence> seqs = new ArrayList<Sequence>();
		DoubleList weight = new DoubleList();
		Sequence seq;
		for (int j = 0; j < wsf.getNumberOfElements(); j++) {
			seq = wsf.getElementAt(j);
			if (seq.equals(seq.reverseComplement())) {
				seqs.add(seq);
				weight.add(wsf.getWeight(j) / div);
			} else {
				if (seq.compareTo(seq.reverseComplement()) < 0) {
					seqs.add(seq);
					weight.add(wsf.getWeight(j));
				}
			}
		}
		return new WeightedDataSetFactory(so, new DataSet(null, seqs.toArray(new Sequence[0])), weight.toArray());
	}
}