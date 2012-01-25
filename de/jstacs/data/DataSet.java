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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.DataSet.WeightedDataSetFactory.SortOperation;
import de.jstacs.data.sequences.ArbitrarySequence;
import de.jstacs.data.sequences.ByteSequence;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.ShortSequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.NullSequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotationParser;
import de.jstacs.io.AbstractStringExtractor;
import de.jstacs.io.SymbolExtractor;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Pair;

/**
 * This is the class for any data set of {@link Sequence}s. All {@link Sequence}s
 * in a {@link DataSet} have to have the same {@link AlphabetContainer}. The
 * {@link Sequence}s may have different lengths.
 * 
 * <br>
 * 
 * For the internal representation the class {@link Sequence} is used, where the
 * external alphabet is converted to integral numerical values. The class
 * {@link DataSet} knows about this coding via instances of class
 * {@link AlphabetContainer} and accordingly {@link de.jstacs.data.alphabets.Alphabet}.
 * 
 * <br>
 * <br>
 * 
 * <a name="access"> There are different ways to access the elements of a
 * {@link DataSet}. If one needs random access there is the method
 * {@link #getElementAt(int)}. For fast sequential access it is recommended to
 * use an {@link ElementEnumerator}. </a>
 * 
 * <br>
 * <br>
 * 
 * {@link DataSet} is immutable.
 * 
 * @author Jens Keilwagen, Andre Gohr, Jan Grau
 * 
 * @see AlphabetContainer
 * @see de.jstacs.data.alphabets.Alphabet
 * @see Sequence
 */
public class DataSet implements Iterable<Sequence>{

	/**
	 * This <code>enum</code> defines different partition methods for a
	 * {@link DataSet}.
	 * 
	 * @author Jens Keilwagen
	 * 
	 * @see DataSet#partition(PartitionMethod, double...)
	 * @see DataSet#partition(int, PartitionMethod)
	 * @see DataSet#partition(double, PartitionMethod, int)
	 * @see DataSet#partition(double[], PartitionMethod, double...)
	 * @see DataSet#partition(double[], int, PartitionMethod)
	 */
	public static enum PartitionMethod {
		/**
		 * This value indicates that the {@link DataSet} will be split by the
		 * number of elements.
		 */
		PARTITION_BY_NUMBER_OF_ELEMENTS,

		/**
		 * This value indicates that the {@link DataSet} will be split by the
		 * number of &quot;symbols&quot;.
		 */
		PARTITION_BY_NUMBER_OF_SYMBOLS,
		
		/**
		 * This value indicates that the {@link DataSet} will be split by weights of the sequences.
		 * If this is not possible (e.g. <code>weights == null</code>), it works like {@link #PARTITION_BY_NUMBER_OF_ELEMENTS}.
		 * 
		 * @see DataSet#partition(double[], PartitionMethod, double...)
		 * @see DataSet#partition(double[], int, PartitionMethod)
		 */
		PARTITION_BY_WEIGHTS;
	}

	/**
	 * Returns the annotation for an array of {@link DataSet}s.
	 * 
	 * @param s
	 *            an array of {@link DataSet}s
	 * 
	 * @return the annotation
	 * 
	 * @see DataSet#getAnnotation()
	 */
	public static final String getAnnotation( DataSet... s ) {
		if( s == null || s.length == 0 ) {
			return "[]";
		} else {
			StringBuffer sb = new StringBuffer( s.length * 100 );
			sb.append( s[0].getAnnotation() );
			for( int i = 1; i < s.length; i++ ) {
				sb.append( ", " );
				sb.append( s[i].getAnnotation() );
			}
			return "[" + sb.toString() + "]";
		}
	}
	
	/**
	 * This method computes the difference between the {@link DataSet} <code>data</code> and
	 * the {@link DataSet}s <code>samples</code>.
	 * 
	 * @param data
	 * 			  the minuend
	 * @param samples
	 *            the subtrahends
	 * 
	 * @return the difference
	 * 
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer}s do not match, i.e., if the {@link DataSet}s are from different domains
	 * @throws EmptyDataSetException
	 *             if the difference is empty
	 */
	public static final DataSet diff( DataSet data, DataSet... samples ) throws EmptyDataSetException, WrongAlphabetException {
		Hashtable<Sequence, int[]> hash = new Hashtable<Sequence, int[]>( data.getNumberOfElements() * 2 );
		AlphabetContainer abc = data.getAlphabetContainer();
		int n = 0, i = 0, anz = data.getNumberOfElements();
		Sequence seq;
		int[] occurrence;
		// put all sequence in the hash table
		for( ; n < anz; n++ ) {
			seq = data.getElementAt( n );
			occurrence = hash.get( seq );
			if( occurrence != null ) {
				occurrence[0]++;
			}else {
				hash.put( seq, new int[]{1} );
			}
		}
		// remove 
		for( ; i < samples.length && hash.size() > 0; i++ ) {
			if( !abc.checkConsistency( samples[i].getAlphabetContainer() ) ) {
				throw new WrongAlphabetException( "The data sets do not have the same AlphabetContainer." );
			}
			for( n = 0; n < samples[i].getNumberOfElements(); n++ ) {
				seq = samples[i].getElementAt( n );
				occurrence = hash.get( seq );
				if( occurrence != null ) {
					if( occurrence[0] == 1 ) {
						hash.remove( seq );
					} else {
						occurrence[0]--;
					}
					anz--;
				}
			}
		}
		// create new sample
		Sequence[] seqs = new Sequence[anz];
		Iterator<Entry<Sequence,int[]>> it = hash.entrySet().iterator();
		Entry<Sequence,int[]> current;
		anz = 0;
		while( it.hasNext() ) {
			current = it.next();
			seq = current.getKey();
			occurrence = current.getValue();
			for( i = 0; i < occurrence[0]; i++ ) {
				seqs[anz++] = seq;
			}
		}
		return new DataSet( "diff of " + data.getAnnotation() + " and " + getAnnotation( samples ), seqs );
	}

	/**
	 * This method computes the intersection between all elements/{@link DataSet}
	 * s of the array, i.e. it returns a {@link DataSet} containing only
	 * {@link Sequence}s that are contained in all {@link DataSet}s of the array.
	 * 
	 * @param samples
	 *            the array of {@link DataSet}s
	 * 
	 * @return the intersection of the elements/{@link DataSet}s in the array
	 * 
	 * @throws IllegalArgumentException
	 *             if the elements of the array are from different domains
	 * @throws EmptyDataSetException
	 *             if the intersection is empty
	 */
	public static final DataSet intersection( DataSet... samples ) throws IllegalArgumentException, EmptyDataSetException {
		WeightedDataSetFactory[] wsf = new WeightedDataSetFactory[samples.length];
		int[] index = new int[samples.length];
		int i = 0, len = -1;
		double anz;
		AlphabetContainer abc = samples[i].getAlphabetContainer();
		while( i < wsf.length ) {
			if( !abc.checkConsistency( samples[i].getAlphabetContainer() ) ) {
				throw new IllegalArgumentException( "The data sets do not have the same AlphabetContainer." );
			}
			try {
				wsf[i] = new WeightedDataSetFactory( SortOperation.SORT_BY_SEQUENCE, samples[i++] );
			} catch ( WrongAlphabetException doesNotHappen ) {
				RuntimeException r = new RuntimeException( doesNotHappen.getMessage() );
				r.setStackTrace( doesNotHappen.getStackTrace() );
				throw r;
			} catch ( WrongLengthException doesNotHappen ) {
				RuntimeException r = new RuntimeException( doesNotHappen.getMessage() );
				r.setStackTrace( doesNotHappen.getStackTrace() );
				throw r;
			}
		}
		boolean goOn = true, same;
		ArrayList<Sequence> list = new ArrayList<Sequence>( 100 );
		String current, help;
		do {
			// find candidate
			current = wsf[0].getElementAt( index[0] ).toString();
			for( i = 1; i < samples.length; i++ ) {
				help = wsf[i].getElementAt( index[i] ).toString();
				if( current.compareTo( help ) < 0 ) {
					current = help;
				}
			}

			// System.out.print( current + " - " );

			same = true;
			for( i = 0; i < samples.length; i++ ) {
				help = "";
				while( index[i] < wsf[i].getNumberOfElements() && current.compareTo( ( help = wsf[i].getElementAt( index[i] ).toString() ) ) > 0 ) {
					index[i]++;
				}
				same &= current.equals( help );
				// System.out.print( help + "\t" );
			}
			// System.out.println();
			if( same ) {
				anz = wsf[0].getWeight( index[0] );
				for( i = 1; i < samples.length; i++ ) {
					if( anz > wsf[i].getWeight( index[i] ) ) {
						anz = wsf[i].getWeight( index[i] );
					}
					index[i]++;
				}

				if( list.size() == 0 ) {
					len = wsf[0].getElementAt( index[0] ).getLength();
				} else if( len != wsf[0].getElementAt( index[0] ).getLength() ) {
					len = 0;
				}
				for( i = 0; i < anz; i++ ) {
					list.add( wsf[0].getElementAt( index[0] ) );
				}
				index[0]++;
			}

			for( i = 0; i < samples.length; i++ ) {
				if( index[i] == wsf[i].getNumberOfElements() ) {
					goOn = false;
				}
			}
		} while( goOn );
		return new DataSet( abc, list.toArray( new Sequence[0] ), len, "intersection of " + getAnnotation( samples ) );
	}

	/**
	 * This method unites all {@link DataSet}s of the array <code>s</code>
	 * regarding the array <code>in</code>.
	 * 
	 * @param s
	 *            the array of {@link DataSet}s
	 * @param in
	 *            an array indicating which {@link DataSet} is used in the union,
	 *            if <code>in[i]==true</code> the {@link DataSet}
	 *            <code>s[i]</code> is used
	 * 
	 * @return the united {@link DataSet}
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>s.length != in.length</code> or the {@link de.jstacs.data.alphabets.Alphabet}
	 *             s do not match
	 * @throws EmptyDataSetException
	 *             if the union is empty
	 * 
	 * @see DataSet#union(DataSet[], boolean[], int)
	 */
	public static final DataSet union( DataSet[] s, boolean[] in ) throws IllegalArgumentException, EmptyDataSetException {
		try {
			return union( s, in, 0 );
		} catch ( WrongLengthException doesNotHappen ) {
			IllegalArgumentException i = new IllegalArgumentException( doesNotHappen.getMessage() );
			i.setStackTrace( doesNotHappen.getStackTrace() );
			throw i;
		}
	}

	/**
	 * Unites all {@link DataSet}s of the array <code>s</code>.
	 * 
	 * @param s
	 *            the array of {@link DataSet}s
	 * 
	 * @return the united {@link DataSet}
	 * 
	 * @throws IllegalArgumentException
	 *             if the {@link de.jstacs.data.alphabets.Alphabet}s do not match
	 * 
	 * @see DataSet#union(DataSet[], boolean[])
	 */
	public static final DataSet union( DataSet... s ) throws IllegalArgumentException {
		if( s == null || s.length == 0 ) {
			return null;
		} else {
			boolean[] in = new boolean[s.length];
			Arrays.fill( in, true );
			try {
				return union( s, in );
			} catch ( EmptyDataSetException doesNotHappen ) {
				// since each given sample is not empty, the union can't be empty
				return null;
			}
		}
	}

	/**
	 * This method unites all {@link DataSet}s of the array <code>s</code>
	 * regarding the array <code>in</code> and sets the element length in the
	 * united {@link DataSet} to <code>subsequenceLength</code>.
	 * 
	 * @param s
	 *            the array of {@link DataSet}s
	 * @param in
	 *            an array indicating which {@link DataSet} is used in the union,
	 *            if <code>in[i]==true</code> the {@link DataSet}
	 *            <code>s[i]</code> is used
	 * @param subsequenceLength
	 *            the length of the elements in the united {@link DataSet}
	 * 
	 * @return the united {@link DataSet}
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>s.length != in.length</code> or the {@link de.jstacs.data.alphabets.Alphabet}
	 *             s do not match
	 * @throws EmptyDataSetException
	 *             if the union is empty
	 * @throws WrongLengthException
	 *             if the united {@link DataSet} does not support this
	 *             <code>subsequenceLength</code>
	 */
	public static final DataSet union( DataSet[] s, boolean[] in, int subsequenceLength ) throws IllegalArgumentException,
			EmptyDataSetException,
			WrongLengthException {
		if( s == null || s.length == 0 ) {
			return null;
		} else {
			if( in.length != s.length ) {
				throw new IllegalArgumentException( "The arrays have to have the same dimension." );
			}
			int i = 0, l = s.length, len, anz, start;
			while( i < l && !in[i] ) {
				i++;
			}
			if( i == l ) {
				return null;
			}
			start = i;
			len = s[i].getElementLength();
			anz = s[i].getNumberOfElements();
			String annot = "the union of [" + s[i].getAnnotation();
			boolean seq = s[i++].indexOfFirstSubseq == null, subseq = !seq;
			while( i < l && ( !in[i] || s[start].alphabetContainer.checkConsistency( s[i].alphabetContainer ) ) ) {
				if( in[i] ) {
					if( s[i].getMinimalElementLength() < subsequenceLength ) {
						// would be asked later anyway, but saves time to do now
						throw new WrongLengthException( subsequenceLength );
					}
					anz += s[i].getNumberOfElements();
					if( len != 0 && len != s[i].getElementLength() ) {
						len = 0;
					}
					seq &= s[i].indexOfFirstSubseq == null;
					subseq &= s[i].indexOfFirstSubseq != null;
					annot += ", " + s[i].getAnnotation();
				}
				i++;
			}
			if( i < l ) {
				throw new IllegalArgumentException( "The alphabets of the data sets do not match." );
			}

			Sequence[] seqs = new Sequence[anz];
			ElementEnumerator ei;
			anz = 0;
			for( i = 0; i < l; i++ ) {
				if( in[i] ) {
					ei = new ElementEnumerator( s[i] );
					while( ei.hasMoreElements() ) {
						seqs[anz++] = ei.nextElement();
					}
				}
			}

			DataSet res = new DataSet( s[start].alphabetContainer, seqs, len, annot + "]" );
			res.setSubsequenceLength( subsequenceLength );

			return res;
		}
	}

	/**
	 * This method unites all {@link DataSet}s of the array <code>s</code> and
	 * sets the element length in the united sample to
	 * <code>subsequenceLength</code>.
	 * 
	 * @param s
	 *            the array of {@link DataSet}s
	 * @param subsequenceLength
	 *            the length of the elements in the united {@link DataSet}
	 * 
	 * @return the united {@link DataSet}
	 * 
	 * @throws IllegalArgumentException
	 *             if the {@link de.jstacs.data.alphabets.Alphabet}s do not match
	 * @throws WrongLengthException
	 *             if the united {@link DataSet} does not support this
	 *             <code>subsequenceLength</code>
	 * 
	 * @see DataSet#union(DataSet[], boolean[], int)
	 */
	public static final DataSet union( DataSet[] s, int subsequenceLength ) throws IllegalArgumentException, WrongLengthException {
		if( s == null || s.length == 0 ) {
			return null;
		} else {
			boolean[] in = new boolean[s.length];
			Arrays.fill( in, true );
			try {
				return union( s, in, subsequenceLength );
			} catch ( EmptyDataSetException doesNotHappen ) {
				// since each given sample is not empty, the union can't be empty
				return null;
			}
		}
	}

	/**
	 * Some annotation for the {@link DataSet}.
	 */
	private String annotation;

	/**
	 * The {@link AlphabetContainer} of the {@link DataSet}.
	 */
	private AlphabetContainer alphabetContainer;

	/**
	 * All sequences of the {@link DataSet}.
	 */
	private Sequence[] seqs;

	/**
	 * The length of the elements.
	 */
	private int length;

	/**
	 * The index of the first subsequence for each {@link Sequence}. Since for
	 * {@link Sequence} 0 the index is always 0 the indices are shifted one
	 * position. That is why the last entry of the array is empty and is used
	 * for the number of elements.
	 */
	private int[] indexOfFirstSubseq;

	/**
	 * Creates a new {@link DataSet} from an {@link AlphabetContainer}, an array
	 * of {@link Sequence}s, a given length for the elements of the
	 * {@link DataSet} and an annotation for the {@link DataSet}.<br>
	 * This constructor is for the <code>partition</code>- and
	 * <code>union</code>-methods. You can decide whether to copy the
	 * {@link Sequence}s in a new array or not.
	 * 
	 * <br>
	 * <br>
	 * 
	 * <b>You have to ensure that all {@link Sequence}s are defined over the
	 * {@link de.jstacs.data.alphabets.Alphabet}(s) in the {@link AlphabetContainer} <code>abc</code>
	 * since it is not checked internally! Furthermore you have to ensure that
	 * the length is correct!</b>
	 * 
	 * @param abc
	 *            the {@link AlphabetContainer} with the {@link de.jstacs.data.alphabets.Alphabet}s
	 * @param seqs
	 *            the array of {@link Sequence}s
	 * @param length
	 *            the length of the {@link Sequence}s
	 * 
	 * @throws EmptyDataSetException
	 *             if the array <code>seqs</code> is <code>null</code> or its
	 *             length is 0
	 */
	private DataSet( AlphabetContainer abc, Sequence[] seqs, int length, String annotation ) throws EmptyDataSetException {
		if( seqs == null || seqs.length == 0 ) {
			throw new EmptyDataSetException();
		}
		this.alphabetContainer = abc;
		this.seqs = seqs;
		this.length = length;
		this.annotation = annotation;
	}

	/**
	 * Creates a new {@link DataSet} from a {@link de.jstacs.io.StringExtractor}
	 * using the given {@link AlphabetContainer}.
	 * 
	 * @param abc
	 *            the {@link AlphabetContainer}
	 * @param se
	 *            the {@link de.jstacs.io.StringExtractor}
	 * 
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} is not suitable
	 * @throws EmptyDataSetException
	 *             if the {@link DataSet} would be empty
	 * @throws WrongLengthException
	 *             never happens (forwarded from
	 *             {@link DataSet#DataSet(AlphabetContainer, AbstractStringExtractor, String, int)}
	 *             )
	 * 
	 * @see DataSet#DataSet(AlphabetContainer, AbstractStringExtractor, String, int)
	 */
	public DataSet( AlphabetContainer abc, AbstractStringExtractor se ) throws WrongAlphabetException, EmptyDataSetException,
																		WrongLengthException {
		this( abc, se, abc.getDelim(), 0 );
	}

	/**
	 * Creates a new {@link DataSet} from a {@link de.jstacs.io.StringExtractor}
	 * using the given {@link AlphabetContainer} and all overlapping windows of
	 * length <code>subsequenceLength</code>.
	 * 
	 * @param abc
	 *            the {@link AlphabetContainer}
	 * @param se
	 *            the {@link de.jstacs.io.StringExtractor}
	 * @param subsequenceLength
	 *            the length of the window sliding on the {@link String} of
	 *            <code>se</code>, if <code>len</code> is 0 (zero) then the
	 *            {@link Sequence}s are used as given from the
	 *            {@link de.jstacs.io.StringExtractor}
	 * 
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} is not suitable
	 * @throws WrongLengthException
	 *             if the subsequence length is not supported
	 * @throws EmptyDataSetException
	 *             if the {@link DataSet} would be empty
	 * 
	 * @see DataSet#DataSet(AlphabetContainer, AbstractStringExtractor, String, int)
	 */
	public DataSet( AlphabetContainer abc, AbstractStringExtractor se, int subsequenceLength ) throws WrongAlphabetException,
																								WrongLengthException, EmptyDataSetException {
		this( abc, se, abc.getDelim(), subsequenceLength );
	}

	/**
	 * Creates a new {@link DataSet} from a {@link de.jstacs.io.StringExtractor}
	 * using the given {@link AlphabetContainer} and a delimiter
	 * <code>delim</code>.
	 * 
	 * @param abc
	 *            the {@link AlphabetContainer}
	 * @param se
	 *            the {@link de.jstacs.io.StringExtractor}
	 * @param delim
	 *            the delimiter for parsing the {@link String}s
	 * 
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} is not suitable
	 * @throws EmptyDataSetException
	 *             if the {@link DataSet} would be empty
	 * @throws WrongLengthException
	 *             never happens (forwarded from
	 *             {@link DataSet#DataSet(AlphabetContainer, AbstractStringExtractor, String, int)}
	 *             )
	 * 
	 * @see DataSet#DataSet(AlphabetContainer, AbstractStringExtractor, String,
	 *      int)
	 */
	public DataSet( AlphabetContainer abc, AbstractStringExtractor se, String delim ) throws WrongAlphabetException, EmptyDataSetException,
																					WrongLengthException {
		this( abc, se, delim, 0 );
	}

	/**
	 * Creates a new {@link DataSet} from a {@link de.jstacs.io.StringExtractor}
	 * using the given {@link AlphabetContainer}, the given delimiter
	 * <code>delim</code> and all overlapping windows of length
	 * <code>subsequenceLength</code>.
	 * 
	 * @param abc
	 *            the {@link AlphabetContainer}
	 * @param se
	 *            the {@link de.jstacs.io.StringExtractor}
	 * @param delim
	 *            the delimiter for parsing the {@link String}s
	 * @param subsequenceLength
	 *            the length of the window sliding on the {@link String} of
	 *            <code>se</code>, if <code>len</code> is 0 (zero) then the
	 *            {@link Sequence}s are used as given from the
	 *            {@link de.jstacs.io.StringExtractor}
	 * 
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} is not suitable
	 * @throws EmptyDataSetException
	 *             if the {@link DataSet} would be empty
	 * @throws WrongLengthException
	 *             if the subsequence length is not supported
	 */
	public DataSet( AlphabetContainer abc, AbstractStringExtractor se, String delim, int subsequenceLength ) throws EmptyDataSetException,
																											WrongAlphabetException,
																											WrongLengthException {
		alphabetContainer = abc;
		LinkedList<Sequence> newSeqs = new LinkedList<Sequence>();
		SymbolExtractor temp = new SymbolExtractor( delim );
		length = -1;
		SequenceAnnotation[] annot;
		
		try {
			if( alphabetContainer.isDiscrete() ) {
				// create pure discrete sample
				int l = (int)alphabetContainer.getMaximalAlphabetLength();
				if( l <= Byte.MAX_VALUE ) {
					while( se.hasMoreElements() ) {
						annot = se.getCurrentSequenceAnnotations();
						temp.setStringToBeParsed( se.nextElement() );
						if( length < 0 ) {
							length = temp.countElements();
						} else if( length > 0 && temp.countElements() != length ) {
							length = 0;
						}
						newSeqs.add( new ByteSequence( alphabetContainer, annot, temp ) );
					}
				} else if( l <= Short.MAX_VALUE ) {
					while( se.hasMoreElements() ) {
						annot = se.getCurrentSequenceAnnotations();
						temp.setStringToBeParsed( se.nextElement() );
						if( length < 0 ) {
							length = temp.countElements();
						} else if( length > 0 && temp.countElements() != length ) {
							length = 0;
						}
						newSeqs.add( new ShortSequence( alphabetContainer, annot, temp ) );
					}
				} else if( l <= Integer.MAX_VALUE ) {
					while( se.hasMoreElements() ) {
						annot = se.getCurrentSequenceAnnotations();
						temp.setStringToBeParsed( se.nextElement() );
						if( length < 0 ) {
							length = temp.countElements();
						} else if( length > 0 && temp.countElements() != length ) {
							length = 0;
						}
						newSeqs.add( new IntSequence( alphabetContainer, annot, temp ) );
					}
				} else {
					throw new WrongAlphabetException( "Could not encode. Too many symbols." );
				}
			} else {
				// create hybrid or pure continuous sample
				if( delim.length() == 0 ) {
					throw new IllegalArgumentException( "delim has to be not empty" );
				}
				while( se.hasMoreElements() ) {
					annot = se.getCurrentSequenceAnnotations();
					temp.setStringToBeParsed( se.nextElement() );
					if( length < 0 ) {
						length = temp.countElements();
					} else if( length > 0 && temp.countElements() != length ) {
						length = 0;
					}
					newSeqs.add( new ArbitrarySequence( alphabetContainer, annot, temp ) );
				}
			}
		} catch ( WrongSequenceTypeException e ) {
			RuntimeException doesNotHappen = new RuntimeException( e.getMessage() );
			doesNotHappen.setStackTrace( e.getStackTrace() );
			throw doesNotHappen;
		}
		seqs = new Sequence[newSeqs.size()];
		if( seqs.length == 0 ) {
			throw new EmptyDataSetException();
		}
		newSeqs.toArray( seqs );

		setSubsequenceLength( subsequenceLength );
		if( subsequenceLength > 0 ) {
			this.annotation = "all subsequences of length " + subsequenceLength + " from " + se.getAnnotation();
		} else {
			this.annotation = se.getAnnotation();
		}
	}

	/**
	 * Creates a new {@link DataSet} from a given {@link DataSet} and a given
	 * length <code>subsequenceLength</code>.<br>
	 * This constructor enables you to use subsequences of the elements of a
	 * {@link DataSet}.
	 * 
	 * <br>
	 * <br>
	 * 
	 * It can also be used to ensure that all sequences that can be accessed by
	 * {@link #getElementAt(int)} are real objects and do not have to be created
	 * at the invocation of the method. (The same holds for the
	 * {@link ElementEnumerator}. In those cases both ways to access the
	 * {@link Sequence} are approximately equally fast.)
	 * 
	 * @param s
	 *            the given {@link DataSet}
	 * @param subsequenceLength
	 *            the new element length
	 * 
	 * @throws WrongLengthException
	 *             if something is wrong with <code>subsequenceLength</code>
	 */
	public DataSet( DataSet s, int subsequenceLength ) throws WrongLengthException {
		this( s, subsequenceLength, false );
	}

	/**
	 * Creates a new {@link DataSet} from a given {@link DataSet} and a given
	 * length <code>subsequenceLength</code>. The given value for
	 * <code>copy</code> determines if the subsequences shall be copied.<br>
	 * This constructor enables you to use subsequences of the elements of a
	 * {@link DataSet}.
	 * 
	 * <br>
	 * <br>
	 * 
	 * It can also be used to ensure that all {@link Sequence}s that can be
	 * accessed by {@link #getElementAt(int)} are real objects and do not have
	 * to be created at the invocation of the method. (The same holds for the
	 * {@link ElementEnumerator}. In those cases both ways to access the
	 * {@link Sequence} are approximately equally fast.)
	 * 
	 * If <code>copy</code> is <code>true</code> all subsequences are copied to
	 * form a new {@link DataSet}.
	 * 
	 * @param s
	 *            the {@link DataSet}
	 * @param subsequenceLength
	 *            the new element length
	 * @param copy
	 *            indicates if the subsequences shall be copied
	 * 
	 * @throws WrongLengthException
	 *             if something is wrong with <code>subsequenceLength</code>
	 */
	private DataSet( DataSet s, int subsequenceLength, boolean copy ) throws WrongLengthException {
		if( copy ) {
			this.alphabetContainer = s.alphabetContainer;
			this.seqs = s.getAllElements();
			setSubsequenceLength( subsequenceLength );
			this.seqs = this.getAllElements();
			this.indexOfFirstSubseq = null;
			this.length = subsequenceLength;
			this.annotation = "all subsequences of length " + subsequenceLength + " from " + s.annotation;
		} else {
			this.alphabetContainer = s.alphabetContainer;
			if( s.indexOfFirstSubseq == null ) {
				this.seqs = s.seqs;
			} else {
				this.seqs = s.getAllElements();
			}
			this.length = s.length;
			setSubsequenceLength( subsequenceLength );
			this.annotation = "all subsequences of length " + subsequenceLength + " from " + s.annotation;
		}

	}

	/**
	 * Creates a new {@link DataSet} from an array of {@link Sequence}s and a
	 * given annotation.<br>
	 * This constructor is specially designed for the method
	 * {@link de.jstacs.sequenceScores.statisticalModels.StatisticalModel#emitDataSet(int, int...)}
	 * 
	 * @param annotation
	 *            the annotation of the {@link DataSet}
	 * @param seqs
	 *            the {@link Sequence}(s)
	 * 
	 * @throws EmptyDataSetException
	 *             if the array <code>seqs</code> is <code>null</code> or the
	 *             length is 0
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer}s do not match
	 */
	public DataSet( String annotation, Sequence... seqs ) throws EmptyDataSetException, WrongAlphabetException {
		if( seqs == null || seqs.length == 0 ) {
			throw new EmptyDataSetException();
		}
		this.alphabetContainer = seqs[0].getAlphabetContainer();
		this.seqs = new Sequence[seqs.length];
		int i = 1;
		length = seqs[0].getLength();
		this.seqs[0] = seqs[0];
		while( i < seqs.length ) {
			this.seqs[i] = seqs[i];
			if( length != seqs[i].getLength() ) {
				length = 0;
			}
			if( !alphabetContainer.checkConsistency( seqs[i++].getAlphabetContainer() ) ) {
				throw new WrongAlphabetException( "The sequences are not defined over the same AlphabetContainer." );
			}
		}
		indexOfFirstSubseq = null;
		this.annotation = annotation;
	}

	/**
	 * Returns an array of {@link Sequence}s containing all elements of this
	 * {@link DataSet}.
	 * 
	 * @return all elements ({@link Sequence}s) of this {@link DataSet}
	 * 
	 * @see ElementEnumerator
	 */
	public Sequence[] getAllElements() {
		Sequence[] res = new Sequence[getNumberOfElements()];
		ElementEnumerator ei = new ElementEnumerator( this );
		for( int i = 0; i < res.length; i++ ) {
			res[i] = ei.nextElement();
		}
		return res;
	}

	/**
	 * Returns the {@link AlphabetContainer} of this {@link DataSet}.
	 * 
	 * @return the {@link AlphabetContainer} of this {@link DataSet}
	 */
	public final AlphabetContainer getAlphabetContainer() {
		return alphabetContainer;
	}

	/**
	 * Returns some annotation of the {@link DataSet}.
	 * 
	 * @return some annotation of the {@link DataSet}
	 */
	public final String getAnnotation() {
		return annotation;
	}

	/**
	 * This method enables you to use only composite {@link Sequence}s of all
	 * elements in the current {@link DataSet}. Each composite {@link Sequence}
	 * will be build from one corresponding {@link Sequence} in this
	 * {@link DataSet} and all composite {@link de.jstacs.data.sequences.Sequence}s
	 * will be returned in a new {@link DataSet}.
	 * 
	 * @param starts
	 *            the start positions of the chunks
	 * @param lengths
	 *            the lengths of the chunks
	 * 
	 * @return a composite {@link DataSet}
	 * 
	 * @throws IllegalArgumentException
	 *             if either <code>starts</code> or <code>lengths</code> or both
	 *             in combination are not suitable
	 * 
	 * @see Sequence#getCompositeSequence(AlphabetContainer, int[], int[])
	 */
	public final DataSet getCompositeDataSet( int[] starts, int[] lengths ) throws IllegalArgumentException {
		AlphabetContainer abc = alphabetContainer.getCompositeContainer( starts, lengths );
		Sequence[] n = new Sequence[getNumberOfElements()];
		ElementEnumerator ei = new ElementEnumerator( this );
		int i = 0, length = 0;
		while( i < n.length ) {
			n[i++] = ei.nextElement().getCompositeSequence( abc, starts, lengths );
		}
		for( i = 0; i < lengths.length; i++ ) {
			length += lengths[i];
		}
		try {
			return new DataSet( abc, n, length, "composite data set (starts=" + Arrays.toString( starts )
												+ ", lengths="
												+ Arrays.toString( lengths )
												+ ") of "
												+ annotation );
		} catch ( EmptyDataSetException doesNotHappen ) {
			// since the current sample is not empty, a sample of infixes can't be empty
			return null;
		}
	}

	/**
	 * This method returns the element, i.e. the {@link Sequence}, with index
	 * <code>i</code>. <a name="getElementAt"> See also <a href="#access">this
	 * comment</a>.
	 * 
	 * @param i
	 *            the index of the element, i.e. the {@link Sequence}
	 * 
	 * @return the element, i.e. the {@link Sequence}, with index <code>i</code>
	 */
	public Sequence getElementAt( int i ) {
		if( indexOfFirstSubseq == null ) {
			return seqs[i];
		} else {
			int seqInd = getIndexOfSeq( i ), startPos = i - ( ( seqInd == 0 ) ? 0 : indexOfFirstSubseq[seqInd - 1] );
			if( length == 0 ) {
				return seqs[seqInd].getSubSequence( startPos );
			} else {
				return seqs[seqInd].getSubSequence( startPos, length );
			}
		}
	}

	/**
	 * Returns the length of the elements, i.e. the {@link Sequence}s, in this
	 * {@link DataSet}.
	 * 
	 * @return the length of the elements, i.e. the {@link Sequence}s, in this
	 *         {@link DataSet}
	 */
	public int getElementLength() {
		return length;
	}
	
	/**
	 * Returns the average length of all {@link Sequence}s in this {@link DataSet}.
	 * @return the average length
	 */
	public double getAverageElementLength(){
		if(length != 0){
			return length;
		}else{
			double meanLength = 0;
			for(int i=0;i<seqs.length;i++){
				meanLength += seqs[i].getLength();
			}
			meanLength /= seqs.length;
			return meanLength;
		}
	}

	/**
	 * This method enables you to use only an infix of all elements, i.e. the
	 * {@link Sequence}s, in the current {@link DataSet}. The subsequences will
	 * be returned in an new {@link DataSet}.
	 * 
	 * <br>
	 * <br>
	 * 
	 * This method can also be used to create a {@link DataSet} of prefixes if
	 * the element length is not zero.
	 * 
	 * @param start
	 *            the start position of the infix
	 * @param length
	 *            the length of the infix, has to be positive
	 * 
	 * @return a {@link DataSet} of the specified infixes
	 * 
	 * @throws IllegalArgumentException
	 *             if either <code>start</code> or <code>length</code> or both
	 *             in combination are not suitable
	 */
	public final DataSet getInfixDataSet( int start, int length ) throws IllegalArgumentException {
		if( length <= 0 ) {
			throw new IllegalArgumentException( "The length has to be positive." );
		}
		int i = this.length == 0 ? getMinimalElementLength() : this.length;
		if( i >= start + length ) {
			if( start == 0 && length == this.length ) {
				return this;
			} else {
				AlphabetContainer abc = alphabetContainer.getSubContainer( start, length );
				Sequence[] n = new Sequence[getNumberOfElements()];
				ElementEnumerator ei = new ElementEnumerator( this );
				for( i = 0; i < n.length; i++ ) {
					n[i] = ei.nextElement().getSubSequence( abc, start, length );
				}
				try {
					return new DataSet( abc, n, length, "infix data set (start=" + start + ", length=" + length + ") of " + annotation );
				} catch ( EmptyDataSetException doesNotHappen ) {
					// since the current sample is not empty, a sample of infixes can't be empty
					return null;
				}
			}
		} else {
			throw new IllegalArgumentException( "The values for start and length or not suitable." );
		}
	}

	/**
	 * Returns a {@link DataSet} that contains the reverse complement of all {@link Sequence}s in
	 * this {@link DataSet}.
	 * @return the reverse complements
	 * @throws OperationNotSupportedException if the {@link AlphabetContainer} of any of the {@link Sequence}s in this {@link DataSet}
	 * 						is not complementable
	 */
	public DataSet getReverseComplementaryDataSet() throws OperationNotSupportedException{
		Sequence[] rc = new Sequence[seqs.length];
		for(int i=0;i<seqs.length;i++){
			rc[i] = seqs[i].reverseComplement();
		}
		try{
			return new DataSet(annotation == null ? null : "reverse complement of "+annotation,rc);
		}catch (EmptyDataSetException e) {
			//cannot happen since this sample is not empty
			return null;
		}catch(WrongAlphabetException ex){
			//cannot happen since alphabet constructed
			return null;
		}
	}
	
	/**
	 * Returns the minimal length of an element, i.e. a {@link Sequence}, in
	 * this {@link DataSet}.
	 * 
	 * @return the minimal length of an element, i.e. a {@link Sequence}, in
	 *         this {@link DataSet}
	 */
	public int getMinimalElementLength() {
		if( length != 0 ) {
			return length;
		} else if( indexOfFirstSubseq != null ) {
			return 0;
		} else {
			int l, min = Integer.MAX_VALUE;
			ElementEnumerator ei = new ElementEnumerator( this );
			while( ei.hasMoreElements() && min != 0 ) {
				l = ei.nextElement().getLength();
				if( l < min ) {
					min = l;
				}
			}
			return min;
		}
	}

	/**
	 * Returns the maximal length of an element, i.e. a {@link Sequence}, in
	 * this {@link DataSet}.
	 * 
	 * @return the maximal length of an element, i.e. a {@link Sequence}, in
	 *         this {@link DataSet}
	 */
	public int getMaximalElementLength() {
		if( length != 0 ) {
			return length;
		} else {
			int l, max = Integer.MIN_VALUE;
			ElementEnumerator ei = new ElementEnumerator( this );
			while( ei.hasMoreElements() ) {
				l = ei.nextElement().getLength();
				if( l > max ) {
					max = l;
				}
			}
			return max;
		}
	}

	/**
	 * Returns the number of elements, i.e. the {@link Sequence}s, in this
	 * {@link DataSet}.
	 * 
	 * @return the number of elements, i.e. the {@link Sequence}s, in this
	 *         {@link DataSet}
	 */
	public int getNumberOfElements() {
		if( indexOfFirstSubseq == null ) {
			return seqs.length;
		} else {
			return indexOfFirstSubseq[seqs.length - 1];
		}
	}
	
	public Iterator<Sequence> iterator() {
		return new ElementEnumerator(this);
	}

	/**
	 * Returns the number of overlapping elements that can be extracted.
	 * 
	 * @param len
	 *            the length of the elements
	 * 
	 * @return the number of elements with the specified length
	 * 
	 * @throws WrongLengthException
	 *             if the given length is bigger than the minimal element length
	 *             
	 * @see #getNumberOfElementsWithLength(int, double[])
	 */
	public int getNumberOfElementsWithLength( int len ) throws WrongLengthException {
		return (int) getNumberOfElementsWithLength( len, null );
	}
	
	/**
	 * Returns the weighted number of overlapping elements that can be extracted.
	 * 
	 * @param len
	 *            the length of the elements
	 * @param weights
	 *            the weights of each element of the sample (see {@link #getElementAt(int)}), can be <code>null</code>
	 * 
	 * @return the weighted number of elements with the specified length
	 * 
	 * @throws WrongLengthException
	 *             if the given length is bigger than the minimal element length
	 * @throws IllegalArgumentException
	 *             if the weights have a wrong dimension
	 */
	public double getNumberOfElementsWithLength( int len, double[] weights ) throws WrongLengthException, IllegalArgumentException {
		double w = 1, all = 0;
		if( weights != null && weights.length != getNumberOfElements() ) {
			throw new IllegalArgumentException( "The weights array has the wrong dimension" );
		}
		if( length == 0 || weights != null ) {
			for( int l, i = 0; i < seqs.length; i++ ) {
				l = seqs[i].getLength();
				if( weights != null ) {
					w = weights[i];
				}
				if( l < len ) {
					throw new WrongLengthException( len );
				} else {
					all += w*(l - len + 1);
				}
			}
		} else {
			if( length < len ) {
				throw new WrongLengthException( len );
			}
			all = ( length - len + 1 ) * seqs.length;
		}
		return all;
	}

	/**
	 * This method enables you to use only a suffix of all elements, i.e. the
	 * {@link Sequence}, in the current {@link DataSet}. The subsequences will be
	 * returned in an new {@link DataSet}.
	 * 
	 * @param start
	 *            the start position of the suffix
	 * 
	 * @return a {@link DataSet} of specified suffixes
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>start</code> is not suitable
	 */
	public final DataSet getSuffixDataSet( int start ) throws IllegalArgumentException {
		int l = 0;
		if( length != 0 ) {
			l = length - start;
		}
		AlphabetContainer abc;
		if( alphabetContainer.isSimple() ) {
			abc = alphabetContainer;
		}
		{
			abc = alphabetContainer.getSubContainer( start, l );
		}
		Sequence[] n = new Sequence[getNumberOfElements()];
		ElementEnumerator ei = new ElementEnumerator( this );
		for( int i = 0; i < n.length; i++ ) {
			n[i] = ei.nextElement().getSubSequence( abc, start );
		}
		try {
			return new DataSet( abc, n, l, "suffix data set (start=" + start + ") of " + annotation );
		} catch ( EmptyDataSetException doesNotHappen ) {
			// since the current sample is not empty, a sample of suffixes can't be empty
			return null;
		}
	}

	/**
	 * This method indicates whether all random variables are defined over the
	 * same range, i.e. all positions use the same (fixed) alphabet.
	 * 
	 * @return <code>true</code> if the {@link DataSet} is simple,
	 *         <code>false</code> otherwise
	 * 
	 * @see AlphabetContainer#isSimple()
	 */
	public final boolean isSimpleDataSet() {
		return alphabetContainer.isSimple();
	}

	/**
	 * This method indicates if all positions use discrete values.
	 * 
	 * @return <code>true</code> if the {@link DataSet} is discrete,
	 *         <code>false</code> otherwise
	 * 
	 * @see AlphabetContainer#isDiscrete()
	 */
	public final boolean isDiscreteDataSet() {
		return alphabetContainer.isDiscrete();
	}

	private static double getSumOfWeights( double[] weights ) {
		double res = 0;
		for( int i = 0; i < weights.length; i++ ) {
			res += weights[i];
		}
		return res;
	}
	
	/**
	 * This method partitions the elements, i.e. the {@link Sequence}s, of the
	 * {@link DataSet} in two distinct parts. The second part (test sample) holds
	 * the percentage of <code>p</code>, the first the rest (train sample). The
	 * first part has element length as the current {@link DataSet}, the second
	 * has element length <code>subsequenceLength</code>, which might be
	 * necessary for testing.
	 * 
	 * @param p
	 *            the percentage for the second part, the second part holds at
	 *            least this percentage of the full {@link DataSet}
	 * @param method
	 *            the method how to partition the sample (partitioning
	 *            criterion)
	 * @param subsequenceLength
	 *            the element length of the second part, if 0 (zero) then the
	 *            sequences are used as given in this {@link DataSet}
	 * 
	 * @return the array of partitioned {@link DataSet}s
	 * 
	 * @throws WrongLengthException
	 *             if something is wrong with <code>subsequenceLength</code>
	 * @throws UnsupportedOperationException
	 *             if the {@link DataSet} is not simple
	 * @throws EmptyDataSetException
	 *             if at least one of the created partitions is empty
	 * 
	 * @see DataSet.PartitionMethod
	 * @see DataSet.PartitionMethod#PARTITION_BY_NUMBER_OF_ELEMENTS
	 * @see DataSet.PartitionMethod#PARTITION_BY_NUMBER_OF_SYMBOLS
	 * @see DataSet#partition(PartitionMethod, double...)
	 */
	public DataSet[] partition( double p, PartitionMethod method, int subsequenceLength ) throws WrongLengthException,
			UnsupportedOperationException,
			EmptyDataSetException {
		if( !isSimpleDataSet() && length != subsequenceLength ) {
			throw new UnsupportedOperationException( "This method can only be used for simple data sets." );
		}
		DataSet[] parts = partition( method, 1d - p, p );
		parts[1].setSubsequenceLength( subsequenceLength );
		return parts;
	}

	/**
	 * This method partitions the elements, i.e. the {@link Sequence}s, of the
	 * {@link DataSet} in distinct parts where each part holds the corresponding
	 * percentage given in the array <code>percentage</code>.
	 * 
	 * @param method
	 *            the method how to partition the {@link DataSet} (partitioning
	 *            criterion)
	 * @param percentage
	 *            the array of percentages for each &quot;subsample&quot;
	 * 
	 * @return the array of partitioned {@link DataSet}s
	 * 
	 * @throws IllegalArgumentException
	 *             if something with the percentages is not correct (
	 *             <code>sum != 1</code> or one value is not in
	 *             <code>[0,1]</code>)
	 * @throws EmptyDataSetException
	 *             if at least one of the created partitions is empty
	 * 
	 * @see DataSet.PartitionMethod
	 * @see DataSet.PartitionMethod#PARTITION_BY_NUMBER_OF_ELEMENTS
	 * @see DataSet.PartitionMethod#PARTITION_BY_NUMBER_OF_SYMBOLS
	 */
	public DataSet[] partition( PartitionMethod method, double... percentage ) throws IllegalArgumentException, EmptyDataSetException {
		return partition( null, method, percentage ).getFirstElement();
	}
	
	/**
	 * This method partitions the elements, i.e. the {@link Sequence}s, of the
	 * {@link DataSet} and the corresponding weights in distinct parts where each part holds the corresponding
	 * percentage given in the array <code>percentage</code>.
	 * 
	 * @param sequenceWeights
	 * 			  the weights for the sequences (might be <code>null</code>)
	 * @param method
	 *            the method how to partition the {@link DataSet} (partitioning
	 *            criterion)
	 * @param percentage
	 *            the array of percentages for each &quot;subsample&quot;
	 * 
	 * @return a {@link Pair} containing an array of partitioned {@link DataSet}s and an array of partitioned sequence weights
	 * 
	 * @throws IllegalArgumentException
	 *             if something with the percentages is not correct (
	 *             <code>sum != 1</code> or one value is not in
	 *             <code>[0,1]</code>)
	 * @throws EmptyDataSetException
	 *             if at least one of the created partitions is empty
	 * 
	 * @see DataSet.PartitionMethod
	 * @see DataSet.PartitionMethod#PARTITION_BY_NUMBER_OF_ELEMENTS
	 * @see DataSet.PartitionMethod#PARTITION_BY_NUMBER_OF_SYMBOLS
	 */
	public Pair<DataSet[],double[][]> partition( double[] sequenceWeights, PartitionMethod method, double... percentage ) throws IllegalArgumentException, EmptyDataSetException {
		if( percentage == null | percentage.length <= 1 ) {
			return new Pair<DataSet[],double[][]>( new DataSet[]{ this }, new double[][]{sequenceWeights} );
		}
		int i = 0;
		double sum = 0;
		for( ; i < percentage.length; i++ ) {
			if( 0 > percentage[i] || 1 < percentage[i] ) {
				throw new IllegalArgumentException( "The value of percentage[" + i + "] is not in [0,1]." );
			}
			sum += percentage[i];
		}
		if( Math.abs( 1d-sum ) > 1E-10 ) {//TODO
			throw new IllegalArgumentException( "The sum of the percentages is not 1. (sum = " + sum + ")" );
		}
		if( sequenceWeights == null && method == PartitionMethod.PARTITION_BY_WEIGHTS ) {
			method = PartitionMethod.PARTITION_BY_NUMBER_OF_ELEMENTS;
		}
		
		double[] anz = new double[percentage.length];
		double l;
		switch( method ) {
			case PARTITION_BY_NUMBER_OF_ELEMENTS:
				l = getNumberOfElements();
				break;
			case PARTITION_BY_NUMBER_OF_SYMBOLS:
				// count all nucleotides
				l = 0;
				ElementEnumerator ei = new ElementEnumerator( this );
				while( ei.hasMoreElements() ) {
					l += ei.nextElement().getLength();
				}
				break;
			case PARTITION_BY_WEIGHTS:
				l = getSumOfWeights(sequenceWeights);
				break;
			default:
				throw new IllegalArgumentException( "The partitioning criterion is unknown." );
		}

		double sumAnz = 0;
		for( i = 0; i < anz.length; i++ ) {
			anz[i] = Math.ceil( l * percentage[i] );
			sumAnz += anz[i];
		}
		i = anz.length - 1;
		double d = Math.ceil( (l-sumAnz)/percentage.length );
		while( i >= 0 && l-sumAnz > d ) {
			anz[i] += d;
			sumAnz += d;
			i--;
		}
		if( i >= 0 ) {
			anz[i] += (l-sumAnz);
		}
		return partitionDataSetAndWeights( anz, method, sequenceWeights );
	}

	/**
	 * This method partitions the elements, i.e. the {@link Sequence}s, of the
	 * {@link DataSet} in <code>k</code> distinct parts.
	 * 
	 * @param k
	 *            the number of distinct parts
	 * @param method
	 *            the method how to partition the {@link DataSet} (partitioning
	 *            criterion)
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>k</code> is not correct
	 * @throws EmptyDataSetException
	 *             if at least one of the created partitions is empty
	 * 
	 * @return the array of partitioned {@link DataSet}s
	 * 
	 * @see DataSet.PartitionMethod
	 * @see DataSet.PartitionMethod#PARTITION_BY_NUMBER_OF_ELEMENTS
	 * @see DataSet.PartitionMethod#PARTITION_BY_NUMBER_OF_SYMBOLS
	 */
	public DataSet[] partition( int k, PartitionMethod method ) throws IllegalArgumentException, EmptyDataSetException {
		return partition(null, k, method).getFirstElement();
	}
	
	/**
	 * This method partitions the elements, i.e. the {@link Sequence}s, of the
	 * {@link DataSet} and the corresponding weights in <code>k</code> distinct parts.
	 *
	 * @param sequenceWeights
	 * 			  the weights for the sequences (might be <code>null</code>)
	 * @param k
	 *            the number of distinct parts
	 * @param method
	 *            the method how to partition the {@link DataSet} (partitioning
	 *            criterion)
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>k</code> is not correct
	 * @throws EmptyDataSetException
	 *             if at least one of the created partitions is empty
	 * 
	 * @return a {@link Pair} containing an array of partitioned {@link DataSet}s and an array of partitioned sequence weights
	 * 
	 * @see DataSet.PartitionMethod
	 * @see DataSet.PartitionMethod#PARTITION_BY_NUMBER_OF_ELEMENTS
	 * @see DataSet.PartitionMethod#PARTITION_BY_NUMBER_OF_SYMBOLS
	 */
	public Pair<DataSet[],double[][]> partition( double[] sequenceWeights, int k, PartitionMethod method ) throws IllegalArgumentException, EmptyDataSetException {
		if( k < 1 ) {
			throw new IllegalArgumentException( "Can't partition in " + k + " parts." );
		}
		if( k == 1 ) {
			return new Pair<DataSet[],double[][]>( new DataSet[]{ this }, new double[][]{sequenceWeights} );
		}
		double[] percentage = new double[k];
		Arrays.fill(percentage,1d/k);
		return partition(sequenceWeights, method, percentage);
	}

	private Pair<DataSet[], double[][]> partitionDataSetAndWeights( double[] anz, PartitionMethod method, double[] seqWeights ) throws EmptyDataSetException {
		int[] pos = new int[getNumberOfElements()], ends = new int[anz.length];
		int last = pos.length, drawn, help, i = 0;
		for( i = 0; i < last; i++ ) {
			pos[i] = i;
		}
		Random r = new Random();

		double current;
		if( seqWeights == null && method == PartitionMethod.PARTITION_BY_WEIGHTS ) {
			method = PartitionMethod.PARTITION_BY_NUMBER_OF_ELEMENTS;
		}
		// for i = 0 it is not necessary to draw
		for( i = anz.length - 1; i > 0; i-- ) {
			ends[i] = last;
			current = 0;
			while( current < anz[i] ) {
				drawn = r.nextInt( last );
				last--;
				help = pos[drawn];
				pos[drawn] = pos[last];
				pos[last] = help;
				switch( method ) {
					case PARTITION_BY_NUMBER_OF_ELEMENTS:
						current++;
						break;
					case PARTITION_BY_NUMBER_OF_SYMBOLS:
						current += getElementAt( help ).getLength();
						break;
					case PARTITION_BY_WEIGHTS:
						current += seqWeights[help];
						break;
					default:
						throw new IllegalArgumentException( "The partitioning criterion is unknown." );
				}
			}
		}
		ends[i] = last;
		return getPartitionsOfElements( pos, ends, seqWeights );
	}

	/**
	 * Randomly samples elements, i.e. {@link Sequence}s, from the set of all
	 * elements, i.e. the {@link Sequence}s, contained in this {@link DataSet}. <br>
	 * Depending on whether this {@link DataSet} is chosen to contain overlapping
	 * elements (windows of length <code>subsequenceLength</code>) or not, those
	 * elements (overlapping windows, whole sequences) are subsampled.
	 * 
	 * @param number
	 *            the number of {@link Sequence}s that should be drawn from the
	 *            contained set of {@link Sequence}s (with replacement)
	 * 
	 * @return a new {@link DataSet} containing the drawn {@link Sequence}s
	 * 
	 * @throws EmptyDataSetException
	 *             if <code>number</code> is not positive
	 */
	public DataSet subSampling( int number ) throws EmptyDataSetException {
		if( number <= 0 ) {
			throw new EmptyDataSetException();
		}
		Random r = new Random();
		Sequence subsampled_seqs[] = new Sequence[number];

		int numOfElements = this.getNumberOfElements();

		for( int i = 0; i < subsampled_seqs.length; i++ ) {
			subsampled_seqs[i] = this.getElementAt( r.nextInt( numOfElements ) );
		}

		return new DataSet( this.alphabetContainer, subsampled_seqs, this.length, "subsample of " + annotation );
	}

	/**
	 * This method writes the {@link DataSet} to a file <code>f</code>.
	 * 
	 * @param f
	 *            the {@link File}
	 * 
	 * @throws IOException
	 *             if something went wrong with the file
	 *             
	 * @see DataSet#save(OutputStream, char, SequenceAnnotationParser)
	 */
	public final void save( File f ) throws IOException {
		save( new FileOutputStream( f ), AbstractStringExtractor.FASTA, null );
	}

	/**
	 * This method allows to write all {@link Sequence}s including their
	 * {@link SequenceAnnotation}s into a {@link OutputStream}. The
	 * {@link SequenceAnnotation}s are parsed using the
	 * {@link SequenceAnnotationParser}.
	 * 
	 * @param stream
	 *            the stream which is used to write the {@link DataSet}
	 * @param commentChar
	 *            the character that marks comment lines
	 * @param p
	 *            the parser for the {@link SequenceAnnotation}s of the
	 *            {@link Sequence}s
	 * 
	 * @throws IOException
	 *             if something went wrong while writing into the stream.
	 * 
	 * @see SequenceAnnotationParser#parseAnnotationToComment(char,
	 *      SequenceAnnotation...)
	 */
	public final void save( OutputStream stream, char commentChar, SequenceAnnotationParser p ) throws IOException {
        BufferedWriter b = new BufferedWriter( new OutputStreamWriter( stream ) );
        ElementEnumerator ei = new ElementEnumerator( this );
        Sequence current;
        if( p == null ) {
            p = NullSequenceAnnotationParser.DEFAULT_INSTANCE;
        }
        while( ei.hasMoreElements() ) {
                current = ei.nextElement();
                b.write( p.parseAnnotationToComment( commentChar, current.getAnnotation() ) );
                b.newLine();
                b.write( current.toString() );
                if( ei.hasMoreElements() ) {
                        b.newLine();
                }
        }
        b.close();
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		ElementEnumerator ei = new ElementEnumerator( this );
		int l = getNumberOfElements();
		StringBuffer erg = new StringBuffer( l * Math.max( getElementLength(), 10 ) );
		erg.append( "annotation       : " + annotation + "\n\n" );
		erg.append( "AlphabetContainer:\n" + alphabetContainer + "\n" );
		erg.append( "element length   : " + getElementLength() + "\n" );
		erg.append( "# of elements    : " + l + "\n\nsequences:\n" );

		Pattern cp = Pattern.compile( "\n" );
		Matcher m = cp.matcher( erg );
		String temp = m.replaceAll( "\n# " );
		erg.delete( 0, erg.length() );
		erg.append( "# " );
		erg.append( temp );
		erg.append( "\n" );

		while( ei.hasMoreElements() ) {
			erg.append( ei.nextElement() + "\n" );
		}
		return erg.toString();
	}

	// finds the index in O(log seqs.length)
	private int getIndexOfSeq( int overAllIndex ) throws IndexOutOfBoundsException {
		if( overAllIndex < 0 || overAllIndex > indexOfFirstSubseq[seqs.length - 1] ) {
			throw new IndexOutOfBoundsException();
		}
		int lower = 0, sep, upper = seqs.length - 1;
		if( overAllIndex < indexOfFirstSubseq[lower] ) {
			return 0;
		} else {
			do {
				sep = ( upper + lower ) / 2;
				if( overAllIndex < indexOfFirstSubseq[sep] ) {
					upper = sep;
				} else {
					lower = sep;
				}
			} while( upper - lower > 1 );
			return lower + 1;
		}
	}

	private Pair<DataSet[],double[][]> getPartitionsOfElements( int[] pos, int[] ends, double[] seqWeights ) throws EmptyDataSetException {
		int i = 0, j, last = 0;
		DataSet[] parts = new DataSet[ends.length];
		double[][] partWeights = new double[ends.length][];
		Sequence[] seqs;
		DoubleList w = new DoubleList();
		for( i = 0; i < ends.length; i++ ) {
			w.clear();
			seqs = new Sequence[ends[i] - last];
			for( j = 0; j < seqs.length; j++ ) {
				seqs[j] = getElementAt( pos[last+j] );
				if( seqWeights != null ) {
					w.add( seqWeights[ pos[last+j] ] );
				}
			}
			last = ends[i];
			partWeights[i] = w.length() == 0 ? null : w.toArray();
			parts[i] = new DataSet( alphabetContainer, seqs, length, "partition of " + annotation );
		}
		return new Pair<DataSet[], double[][]>(parts,partWeights);
	}

	/**
	 * This method computes the indices if one wants to use subsequences of
	 * length <code>len</code>. If <code>len</code> is 0 (zero) then the
	 * {@link Sequence}s are used as given.
	 * 
	 * @param len
	 *            the subsequence length
	 * 
	 * @throws WrongLengthException
	 *             if the subsequence length is not supported
	 */
	private void setSubsequenceLength( int len ) throws WrongLengthException {
		if( len < 0 ) {
			throw new WrongLengthException( len );
		}
		if( length != len ) {
			if( len == 0 ) {
				return;
			}
			if( indexOfFirstSubseq != null ) {
				throw new UnsupportedOperationException( "operation not supported since indexOfFirstSubseq != null" );
			}
			if( isSimpleDataSet() ) {
				indexOfFirstSubseq = new int[seqs.length];
				if( length == 0 ) {
					int l, i = 0, all = 0;
					while( i < seqs.length ) {
						l = seqs[i].getLength();
						if( l < len ) {
							throw new WrongLengthException( len );
						} else {
							all += l - len + 1;
							indexOfFirstSubseq[i++] = all;
						}
					}
				} else {
					if( length < len ) {
						throw new WrongLengthException( len );
					}
					for( int i = 0, offset = length - len + 1, all = offset; i < seqs.length; i++, all += offset ) {
						indexOfFirstSubseq[i] = all;
					}
				}
				length = len;
			} else {
				throw new UnsupportedOperationException( "For this data set it is impossible to have a sliding window, since the AlphabetContainer is not simple." );
			}
		}
	}
	
	/**
	 * This method returns all {@link SequenceAnnotation} types and the corresponding
	 * identifier which occur in this {@link DataSet}.
	 * 
	 * @return a {@link Hashtable} with key = {@link SequenceAnnotation} type and identifier = {@link SequenceAnnotation} identifier
	 * 
	 * @see SequenceAnnotation
	 */
	public Hashtable<String, HashSet<String>> getAnnotationTypesAndIdentifier() {
		Hashtable<String, HashSet<String>> res = new Hashtable<String, HashSet<String>>();
		HashSet<String> help;
		boolean add;
		for( Sequence s : this ) {
			SequenceAnnotation[] seqAn = s.getAnnotation();
			if( seqAn != null ) {
				for( SequenceAnnotation current : seqAn ) {
					help = res.get( current.getType() );
					add = false;
					if( help == null ) {
						add = true;
						help = new HashSet<String>();
					}
					help.add(  current.getIdentifier() );
					if( add ) {
						res.put( current.getType(), help );
					}
				}
			}
		}
		return res;
	}

	/**
	 * This method creates a matrix which contains the index of the {@link Sequence} with specific {@link SequenceAnnotation}
	 * combination or -1 if the {@link DataSet} does not contain any {@link Sequence} with such a combination. The rows and
	 * columns are indexed according to the {@link Hashtable}s.
	 * 
	 * <br><br>
	 * 
	 * Here is a short example, how to interpret the returned matrix:
	 * <code><pre>
	 * int[][] matrix = s.getSequenceAnnotationIndexMatrix( rowType, rowHash, columnType, columnHash )
	 * 
	 * if( matrix[i][j] < 0 ) {
	 * 	System.out.println( "There is no Sequence in the DataSet with this SequenceAnnotation combination");
	 * } else {
	 * 	System.out.println( "This is the Sequence: " + s.getElementAt( matrix[i][j] ) );
	 * }
	 * </pre></code>
	 * 
	 * 
	 * @param rowType the {@link SequenceAnnotation} type for the rows
	 * @param rowHash a {@link Hashtable} of {@link SequenceAnnotation} identifier and indices for the rows 
	 * @param columnType the {@link SequenceAnnotation} type for the columns
	 * @param columnHash a {@link Hashtable} of {@link SequenceAnnotation} identifier and indices for the columns
	 * 
	 * @return a matrix with the indices of the {@link Sequence}s with each specific combination of
	 * 		   {@link SequenceAnnotation} for code>rowType</code> and <code>columnType</code> and -1
	 * 		   if this combination does not exist in the {@link DataSet}
	 * 
	 * @see DataSet#getAnnotationTypesAndIdentifier()
	 * @see de.jstacs.utils.ToolBox#parseHashSet2IndexHashtable(HashSet)
	 */
	public int[][] getSequenceAnnotationIndexMatrix( String rowType, Hashtable<String, Integer> rowHash, String columnType, Hashtable<String, Integer> columnHash ) {
		
		int[][] result = new int[rowHash.size()][columnHash.size()];
		for( int r = 0; r < result.length; r++ ) {
			Arrays.fill( result[r], -1 );
		}
		String help;
		int idxRow, idxColumn, idx = 0;
		for( Sequence s : this ) {
			SequenceAnnotation[] seqAn = s.getAnnotation();
			idxRow = idxColumn = -1;
			if( seqAn != null ) {
				for( SequenceAnnotation current : seqAn ) {
					help = current.getType();
					if( help.equals( rowType ) ) {
						idxRow = rowHash.get( current.getIdentifier() );
					} else if( help.equals( columnType ) ) {
						idxColumn = columnHash.get( current.getIdentifier() );
					}
				}
			}
			if( idxRow >= 0 && idxColumn >= 0 ) {
				result[idxRow][idxColumn] = idx;
			}
			idx++;
		}
		return result;
	}

	/**
	 * This class can be used to have a fast sequential access to a
	 * {@link DataSet}. <a name="ElementEnumerator"> It enumerates all elements
	 * of a {@link DataSet}.
	 * 
	 * <br>
	 * <br>
	 * 
	 * As further functionality the method {@link #reset()} is implemented to
	 * reuse an {@link ElementEnumerator}.
	 * 
	 * @author Jens Keilwagen </a>
	 */
	public static class ElementEnumerator implements RecyclableSequenceEnumerator, Iterator<Sequence> {

		private int seqCounter, startPosCounter;

		private DataSet s;

		/**
		 * Creates a new {@link ElementEnumerator} on the given {@link DataSet}
		 * <code>data</code>.
		 * 
		 * @param data
		 *            the given {@link DataSet}
		 */
		public ElementEnumerator( DataSet data ) {
			s = data;
			reset();
		}

		/* (non-Javadoc)
		 * @see java.util.Enumeration#hasMoreElements()
		 */
		public boolean hasMoreElements() {
			return seqCounter < s.seqs.length;
		}

		/* (non-Javadoc)
		 * @see java.util.Enumeration#nextElement()
		 */
		public Sequence nextElement() {
			if( s.indexOfFirstSubseq == null ) {
				return s.seqs[seqCounter++];
			} else {
				Sequence current;
				if( s.length != 0 ) {
					current = s.seqs[seqCounter].getSubSequence( startPosCounter, s.length );
				} else {
					current = s.seqs[seqCounter].getSubSequence( startPosCounter );
				}
				if( ++startPosCounter + s.length > s.seqs[seqCounter].getLength() ) {
					seqCounter++;
					startPosCounter = 0;
				}
				return current;
			}
		}

		public void reset() {
			seqCounter = startPosCounter = 0;
		}

		public boolean hasNext() {
			return hasMoreElements();
		}

		public Sequence next() {
			return nextElement();
		}

		public void remove() {
			throw new UnsupportedOperationException("DataSets are immutable");
		}
	}

	/**
	 * This class enables you to eliminate {@link Sequence}s that occur more
	 * than once in one or more {@link DataSet}s. The number of occurrences is
	 * given by the weight for a {@link Sequence}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class WeightedDataSetFactory {

		/**
		 * This <code>enum</code> defines the different types of sort operations
		 * that can be performed while creating a {@link WeightedDataSetFactory}.
		 * 
		 * @author Jens Keilwagen
		 */
		public static enum SortOperation {
			/**
			 * Indicates that no sorting shall be done after eliminating
			 * {@link Sequence}s that occur more than once.
			 */
			NO_SORT,

			/**
			 * Indicates that the {@link Sequence}s shall be sorted after
			 * eliminating {@link Sequence}s that occur more than once. Probably
			 * this is the slowest option.
			 */
			SORT_BY_SEQUENCE,

			/**
			 * This value indicates that the {@link Sequence}s shall be sorted
			 * according to their weights after eliminating {@link Sequence}s
			 * that occur more than once.
			 */
			SORT_BY_WEIGHTS;
		}

		private DataSet res;

		private double[] weights;

		/**
		 * Creates a new {@link WeightedDataSetFactory} on the given
		 * {@link DataSet}(s) with {@link SortOperation} <code>sort</code>.
		 * 
		 * @param sort
		 *            the given {@link SortOperation}
		 * @param data
		 *            the given {@link DataSet}(s)
		 * 
		 * @throws WrongAlphabetException
		 *             if the alphabets of the {@link DataSet}s do not match
		 * @throws WrongLengthException
		 *             does not happen (forwarded from
		 *             {@link de.jstacs.data.DataSet.WeightedDataSetFactory#DataSet.WeightedDataSetFactory(DataSet.WeightedDataSetFactory.SortOperation, DataSet[], double[][], int)}
		 *             )
		 * 
		 * @see de.jstacs.data.DataSet.WeightedDataSetFactory#DataSet.WeightedDataSetFactory(DataSet.WeightedDataSetFactory.SortOperation, DataSet[], double[][], int) 
		 */
		public WeightedDataSetFactory( SortOperation sort, DataSet... data ) throws WrongAlphabetException, WrongLengthException {
			this( sort, data, null, 0 );
		}

		/**
		 * Creates a new {@link WeightedDataSetFactory} on the given
		 * {@link DataSet} and an array of <code>weights</code> with
		 * {@link SortOperation} <code>sort</code>.
		 * 
		 * @param sort
		 *            the given {@link SortOperation}
		 * @param data
		 *            the given {@link DataSet}
		 * @param weights
		 *            the weights for each element in the {@link DataSet}
		 * 
		 * @throws WrongAlphabetException
		 *             if the alphabets of the {@link DataSet}s do not match
		 * @throws WrongLengthException
		 *             does not happen (forwarded from
		 *             {@link de.jstacs.data.DataSet.WeightedDataSetFactory#DataSet.WeightedDataSetFactory(DataSet.WeightedDataSetFactory.SortOperation, DataSet[], double[][], int)}
		 *             )
		 * 
		 * @see de.jstacs.data.DataSet.WeightedDataSetFactory#DataSet.WeightedDataSetFactory(DataSet.WeightedDataSetFactory.SortOperation, DataSet[], double[][], int)
		 */
		public WeightedDataSetFactory( SortOperation sort, DataSet data, double[] weights ) throws WrongAlphabetException,
																							WrongLengthException {
			this( sort, new DataSet[]{ data }, new double[][]{ weights }, 0 );
		}

		/**
		 * Creates a new {@link WeightedDataSetFactory} on the given
		 * {@link DataSet} and an array of <code>weights</code> with a given
		 * <code>length</code> and {@link SortOperation} <code>sort</code>.
		 * 
		 * @param sort
		 *            the given {@link SortOperation}
		 * @param data
		 *            the given {@link DataSet}
		 * @param weights
		 *            the weight for each element in the {@link DataSet}
		 * @param length
		 *            the length of the elements in the resulting
		 *            {@link WeightedDataSetFactory}
		 * 
		 * @throws WrongAlphabetException
		 *             if the alphabets of the {@link DataSet}s do not match
		 * @throws WrongLengthException
		 *             if the length is not supported
		 * 
		 * @see de.jstacs.data.DataSet.WeightedDataSetFactory#DataSet.WeightedDataSetFactory(DataSet.WeightedDataSetFactory.SortOperation, DataSet[], double[][], int)
		 */
		public WeightedDataSetFactory( SortOperation sort, DataSet data, double[] weights, int length ) throws WrongAlphabetException,
																										WrongLengthException {
			this( sort, new DataSet[]{ data }, new double[][]{ weights }, length );
		}

		/**
		 * Creates a new {@link WeightedDataSetFactory} on the given array of
		 * {@link DataSet}s and an array of <code>weights</code> with a given
		 * <code>length</code> and {@link SortOperation} <code>sort</code>.
		 * 
		 * @param sort
		 *            the given {@link SortOperation}
		 * @param data
		 *            the given {@link DataSet}
		 * @param weights
		 *            the weights for each element in each {@link DataSet}
		 * @param length
		 *            the length of the elements in the resulting
		 *            {@link WeightedDataSetFactory}
		 * 
		 * @throws WrongAlphabetException
		 *             if the alphabets of the {@link DataSet}s do not match
		 * @throws WrongLengthException
		 *             if the length is not supported
		 */
		public WeightedDataSetFactory( SortOperation sort, DataSet[] data, double[][] weights, int length ) throws WrongAlphabetException,
																											WrongLengthException {
			Hashtable<Sequence, double[]> ht = new Hashtable<Sequence, double[]>( data.length * data[0].getNumberOfElements() );
			for( int i = 0; i < data.length; i++ ) {
				if( data[0].alphabetContainer.checkConsistency( data[i].alphabetContainer ) ) {
					if( weights != null ) {
						add( ht, data[i], weights[i], length );
					} else {
						add( ht, data[i], null, length );
					}
				} else {
					throw new WrongAlphabetException( "The AlphabetContainer for all DataSet has to be consistent." );
				}
			}
			create( "all sequences" + ( length > 0 ? ( " of length " + length ) : "" ) + " that occur in " + DataSet.getAnnotation( data ),
					sort,
					ht );
		}

		private void add( Hashtable<Sequence, double[]> ht, DataSet data, double[] weights, int length ) throws WrongLengthException {
			Sequence s;
			double w = 1;
			int i = 0, anz = data.getNumberOfElements(), j, l;
			for( ; i < anz; i++ ) {
				//System.out.println( i + "\t" + ht.size() );
				s = data.getElementAt( i );
				if( weights != null ) {
					w = weights[i];
				}
				if( length == 0 ) {
					put( ht, s, w );
				} else {
					l = s.getLength() - length + 1;
					if( l > 0 ) {
						for( j = 0; j < l; j++ ) {
							put( ht, s.getSubSequence( s.getAlphabetContainer(), j, length ), w );
						}
					} else {
						throw new WrongLengthException( length );
					}
				}
			}
		}

		private void put( Hashtable<Sequence, double[]> ht, Sequence s, double w ) {
			double[] value = ht.get( s );
			if( value != null ) {
				value[0] += w;
			} else {
				ht.put( s, new double[]{ w } );
			}
		}

		@SuppressWarnings( "unchecked" )
		private void create( String annotation, SortOperation sort, Hashtable<Sequence, double[]> ht ) {
			Entry<Sequence, double[]>[] array = ht.entrySet().toArray( new Entry[0] );
			switch( sort ) {
				case NO_SORT:
					break;
				case SORT_BY_SEQUENCE:
					Arrays.sort( array, SequenceComparator.DEFAULT );
					break;
				case SORT_BY_WEIGHTS:
					Arrays.sort( array, WeightsComparator.DEFAULT );
					break;
				default:
					throw new IllegalArgumentException( "unknown sort operation" );
			}

			Sequence[] seqs = new Sequence[array.length];
			weights = new double[array.length];
			Entry<Sequence, double[]> e;
			for( int i = 0; i < weights.length; i++ ) {
				e = array[i];
				seqs[i] = e.getKey();
				weights[i] = e.getValue()[0];
			}
			try {
				res = new DataSet( annotation, seqs );
			} catch ( Exception doesNotHappen ) {
				RuntimeException r = new RuntimeException( doesNotHappen.getMessage() );
				r.setStackTrace( doesNotHappen.getStackTrace() );
				throw r;
			}
		}

		/**
		 * Returns the {@link Sequence} with index <code>index</code>.
		 * 
		 * @param index
		 *            the index of the {@link Sequence}
		 * 
		 * @return the {@link Sequence} with index <code>index</code>
		 */
		public Sequence getElementAt( int index ) {
			return res.getElementAt( index );
		}

		/**
		 * Returns the number of elements, i.e. {@link Sequence}s, in the
		 * internal {@link DataSet}.
		 * 
		 * @return the number of elements, i.e. {@link Sequence}s, in the
		 *         internal {@link DataSet}
		 */
		public int getNumberOfElements() {
			return res.getNumberOfElements();
		}

		/**
		 * Returns the {@link DataSet}, where each {@link Sequence} occurs only
		 * once.
		 * 
		 * @return the {@link DataSet}, where each {@link Sequence} occurs only
		 *         once
		 */
		public DataSet getDataSet() {
			return res;
		}

		/**
		 * Returns the sum of all weights.
		 * 
		 * @return the sum of all weights
		 */
		public double getSumOfWeights() {
			return DataSet.getSumOfWeights(weights);
		}

		/**
		 * Returns the weight for the {@link Sequence} with index
		 * <code>index</code>.
		 * 
		 * @param index
		 *            the index of the {@link Sequence}
		 * 
		 * @return the weight for the {@link Sequence} with index
		 *         <code>index</code>
		 */
		public double getWeight( int index ) {
			return weights[index];
		}

		/**
		 * Returns a copy of the weights for the {@link DataSet}.
		 * 
		 * @return the weights for the {@link DataSet}
		 */
		public double[] getWeights() {
			return weights.clone();
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#toString()
		 */
		@Override
		public String toString() {
			StringBuffer sb = new StringBuffer( ( 10 + res.getElementLength() ) * weights.length );
			for( int i = 0; i < weights.length; i++ ) {
				sb.append( i + " " + res.getElementAt( i ) + "\t" + weights[i] + "\n" );
			}
			return sb.toString();
		}

		/**
		 * This comparator can be used to sort the elements in a {@link WeightedDataSetFactory}.
		 * It corresponds to {@link SortOperation#SORT_BY_WEIGHTS}
		 * 
		 * @author Jens Keilwagen
		 */
		private static final class WeightsComparator implements Comparator<Entry<Sequence, double[]>> {

			/**
			 * This constant holds the only instance of the class.
			 */
			public static final WeightsComparator DEFAULT = new WeightsComparator();

			private WeightsComparator() {};

			public int compare( Entry<Sequence, double[]> o1, Entry<Sequence, double[]> o2 ) {
				return (int)Math.signum( o2.getValue()[0] - o1.getValue()[0] );
			}
		}

		/**
		 * This comparator can be used to sort the elements in a {@link WeightedDataSetFactory}.
		 * It corresponds to {@link SortOperation#SORT_BY_SEQUENCE}
		 * 
		 * @author Jens Keilwagen
		 */
		private static final class SequenceComparator implements Comparator<Entry<Sequence, double[]>> {

			/**
			 * This constant holds the only instance of the class.
			 */
			public static final SequenceComparator DEFAULT = new SequenceComparator();

			private SequenceComparator() {};

			public int compare( Entry<Sequence, double[]> o1, Entry<Sequence, double[]> o2 ) {
				return o1.getKey().compareTo( o2.getKey() );
			}
		}
	}
}
