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

package de.jstacs.data.sequences;

import java.util.Arrays;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.ComplementableDiscreteAlphabet;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;

/**
 * This is the main class for all sequences. All sequences are immutable.
 * 
 * @param <T> the type of each position
 * 
 * @author Jens Keilwagen
 */
public abstract class Sequence<T> implements Comparable<Sequence<T>> {

	/**
	 * The underlying alphabets.
	 */
	protected AlphabetContainer alphabetCon;

	/**
	 * The pointer to the reverse complement of the {@link de.jstacs.data.sequences.Sequence}.
	 */
	protected Sequence<T> rc;

	/**
	 * The annotation of the {@link de.jstacs.data.sequences.Sequence}.
	 */
	protected SequenceAnnotation[] annotation;
	
	/**
	 * The hash code for the instance.
	 */
	private Integer hashCode;

	/**
	 * Creates a new {@link de.jstacs.data.sequences.Sequence} with the given {@link AlphabetContainer}
	 * and the given annotation, but without the content. The content has to be
	 * set by the constructor of the extending class.
	 * 
	 * @param container
	 *            the {@link AlphabetContainer} of the {@link de.jstacs.data.sequences.Sequence}
	 * @param annotation
	 *            the annotation of the {@link de.jstacs.data.sequences.Sequence}
	 */
	protected Sequence( AlphabetContainer container, SequenceAnnotation[] annotation ) {
		if( container == null ) {
			throw new NullPointerException();
		}
		this.alphabetCon = container;
		if( annotation != null ) {
			this.annotation = annotation.clone();
		}
		hashCode = null;
	}

	/**
	 * Returns the continuous value at position <code>pos</code> of the
	 * {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * @param pos
	 *            the position of the {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @return the continuous value at position <code>pos</code> of the
	 *         {@link de.jstacs.data.sequences.Sequence}
	 */
	public abstract double continuousVal( int pos );

	/**
	 * Returns the discrete value at position <code>pos</code> of the
	 * {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * @param pos
	 *            the position of the {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @return the discrete value at position <code>pos</code> of the
	 *         {@link de.jstacs.data.sequences.Sequence}
	 */
	public abstract int discreteVal( int pos );

	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals( Object o ) {
		if( o == null ) {
			return false;
		}
		if( !( o instanceof Sequence ) ) {
			return false;
		}
		Sequence s = (Sequence)o;
		return compareTo( s ) == 0 && alphabetCon.checkConsistency( s.alphabetCon );
	}

	/**
	 * Return the alphabets, i.e. the {@link AlphabetContainer}, used in this
	 * {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * @return the alphabets, i.e. the {@link AlphabetContainer}, used in this
	 *         {@link de.jstacs.data.sequences.Sequence}
	 */
	public final AlphabetContainer getAlphabetContainer() {
		return alphabetCon;
	}

	/**
	 * Returns the annotation of the {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * @return the annotation of the {@link de.jstacs.data.sequences.Sequence} (can be <code>null</code>)
	 */
	public final SequenceAnnotation[] getAnnotation() {
		if( annotation != null ) {
			return annotation.clone();
		} else {
			return null;
		}
	}
	
	
	/**
	* Returns the {@link SequenceAnnotation} of this {@link de.jstacs.data.sequences.Sequence} that has type <code>type</code> and identifier <code>identifier</code>.
	* 
	* @param type the chosen type of the {@link SequenceAnnotation}
	* @param identifier the chosen identifier of the {@link SequenceAnnotation}
	* @return the first {@link SequenceAnnotation} that meets the criteria
	*/
	public SequenceAnnotation getSequenceAnnotationByTypeAndIdentifier( String type, String identifier ) {
		if( annotation != null ) {
			for(SequenceAnnotation ann : annotation){
				//System.out.println(ann.getType()+" <-> "+type+" ; "+ann.getIdentifier()+" <-> "+identifier);
				if(ann.getType().equals( type ) && ann.getIdentifier().equals( identifier )){
					return ann;
				}
			}
		}
		return null;
	}
	
	/**
	* Returns the {@link SequenceAnnotation} no. <code>idx</code> of this {@link de.jstacs.data.sequences.Sequence} that has type <code>type</code>
	* 
	* @param type the chosen type of a subset of {@link SequenceAnnotation}s
	* @param idx the index of the returned {@link SequenceAnnotation} within this subset.
	* @return the {@link SequenceAnnotation} no. <code>idx</code> with type <code>type</code>
	*/
	public SequenceAnnotation getSequenceAnnotationByType(String type, int idx){
		if( annotation != null ) {
			int i=0;
			for(SequenceAnnotation ann : annotation){
				if(ann.getType().equals( type )){
					if(i==idx){
						return ann;
					}
					i++;
				}
			}
		}
		return null;
	}
	
	/**
	 * Returns the number of {@link SequenceAnnotation}s of type <code>type</code> for this {@link de.jstacs.data.sequences.Sequence}.
	 * @param type the type
	 * @return the number of annotations
	 */
	public int getNumberOfSequenceAnnotationsByType(String type){
		int i=0;
		if( annotation != null ) {
			for(SequenceAnnotation ann : annotation){
				if(ann.getType().equals( type )){
					i++;
				}
			}
		}
		return i;
	}

	/**
	 * This method should be used if one wants to create a {@link de.jstacs.data.DataSet} of
	 * {@link Sequence.CompositeSequence}s. With this constructor you are enabled to
	 * create a {@link de.jstacs.data.DataSet} where every {@link de.jstacs.data.sequences.Sequence} has the same
	 * {@link AlphabetContainer} instance.
	 * 
	 * <br>
	 * <br>
	 * 
	 * Internally it is checked that the {@link AlphabetContainer} matches with
	 * the one of the {@link Sequence.CompositeSequence}.
	 * 
	 * @param abc
	 *            the new {@link AlphabetContainer}
	 * @param starts
	 *            the start positions of the junks
	 * @param lengths
	 *            the length of each junk
	 * 
	 * @return the {@link Sequence.CompositeSequence}
	 * 
	 * @see Sequence.CompositeSequence#Sequence.CompositeSequence(de.jstacs.data.AlphabetContainer,
	 * 			de.jstacs.data.sequences.Sequence, int[], int[])
	 */
	public Sequence<T> getCompositeSequence( AlphabetContainer abc, int[] starts, int[] lengths ) {
		return new CompositeSequence( abc, this, starts, lengths );
	}

	/**
	 * This is a very efficient way to create a {@link Sequence.CompositeSequence} for
	 * sequences with a simple {@link AlphabetContainer}.
	 * 
	 * @param starts
	 *            the start positions of the junks
	 * @param lengths
	 *            the length of each junk
	 * 
	 * @return the {@link Sequence.CompositeSequence}
	 * 
	 * @see Sequence.CompositeSequence#Sequence.CompositeSequence(de.jstacs.data.sequences.Sequence, int[], int[])
	 */
	public Sequence getCompositeSequence( int[] starts, int[] lengths ) {
		return new CompositeSequence( this, starts, lengths );
	}

	/**
	 * This method should be used if one wants to create a {@link de.jstacs.data.DataSet} of
	 * subsequences of defined length. With this constructor you are enabled to
	 * create a {@link de.jstacs.data.DataSet} where every {@link de.jstacs.data.sequences.Sequence} has the same
	 * {@link AlphabetContainer} instance.
	 * 
	 * <br>
	 * <br>
	 * 
	 * Internally it is checked that the {@link AlphabetContainer} matches with
	 * the one of the subsequence.
	 * 
	 * @param abc
	 *            the new {@link AlphabetContainer}
	 * @param start
	 *            the index of the start position
	 * 
	 * @return the subsequence
	 * 
	 * @see Sequence#getSubSequence(de.jstacs.data.AlphabetContainer, int, int)
	 */
	public final Sequence getSubSequence( AlphabetContainer abc, int start ) {
		return getSubSequence( abc, start, getLength() - start );
	}

	/**
	 * This method should be used if one wants to create a {@link de.jstacs.data.DataSet} of
	 * subsequences of defined length. With this constructor you are enabled to
	 * create a {@link de.jstacs.data.DataSet} where every {@link de.jstacs.data.sequences.Sequence} has the same
	 * {@link AlphabetContainer} instance.
	 * 
	 * <br>
	 * <br>
	 * 
	 * Internally it is checked that the {@link AlphabetContainer} matches with
	 * the one of the subsequence.
	 * 
	 * @param abc
	 *            the new {@link AlphabetContainer}
	 * @param start
	 *            the index of the start position
	 * @param length
	 *            the length of the new {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @return the subsequence
	 * 
	 * @see Sequence.SubSequence#Sequence.SubSequence(de.jstacs.data.AlphabetContainer, de.jstacs.data.sequences.Sequence, int, int) SubSequence#SubSequence(de.jstacs.data.AlphabetContainer, de.jstacs.data.Sequence, int, int)
	 */
	public Sequence getSubSequence( AlphabetContainer abc, int start, int length ) {
		if( start == 0 && length == getLength() ) {
			return this;
		} else {
			return new SubSequence( abc, this, start, length );
		}
	}

	/**
	 * This is a very efficient way to create a subsequence/suffix for
	 * {@link de.jstacs.data.sequences.Sequence}s with a simple {@link AlphabetContainer}.
	 * 
	 * @param start
	 *            the index of the start position
	 * 
	 * @return the subsequence
	 * 
	 * @see Sequence#getSubSequence(int, int)
	 */
	public final Sequence getSubSequence( int start ) {
		return getSubSequence( start, getLength() - start );
	}

	/**
	 * This is a very efficient way to create a subsequence of defined length
	 * for {@link de.jstacs.data.sequences.Sequence}s with a simple {@link AlphabetContainer}.
	 * 
	 * @param start
	 *            the index of the start position
	 * @param length
	 *            the length of the new {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @return the subsequence
	 * 
	 * @see Sequence.SubSequence#Sequence.SubSequence(Sequence, int, int) SubSequence#SubSequence(Sequence, int, int)
	 */
	public Sequence getSubSequence( int start, int length ) {
		if( start == 0 && length == getLength() ) {
			return this;
		} else {
			return new SubSequence( this, start, length );
		}
	}

	/**
	 * This method allows to append annotation to a {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * @param add
	 *            indicates whether to add the new annotation to the existing or
	 *            not
	 * @param annotation
	 *            the new annotation
	 * 
	 * @return the new annotated {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @see Sequence#flatCloneWithoutAnnotation()
	 */
	public Sequence annotate( boolean add, SequenceAnnotation... annotation ) {
		Sequence seq = this.flatCloneWithoutAnnotation();
		if( add && annotation != null ) {
			int num = annotation.length;
			if( this.annotation != null ) {
				num += this.annotation.length;
			}
			SequenceAnnotation[] temp = new SequenceAnnotation[num];
			if( this.annotation != null ) {
				num = this.annotation.length;
				System.arraycopy( this.annotation, 0, temp, 0, num );
			} else {
				num = 0;
			}
			System.arraycopy( annotation, 0, temp, num, annotation.length );
			seq.annotation = temp;
		} else {
			if( annotation != null ) {
				seq.annotation = annotation.clone();
			} else {
				seq.annotation = null;
			}
		}
		return seq;
	}

	/**
	 * Works in analogy to {@link Sequence#clone()}, but does not clone the
	 * annotation. This method is used in
	 * {@link #annotate(boolean, SequenceAnnotation...)}.
	 * 
	 * @return the cloned {@link de.jstacs.data.sequences.Sequence} without annotation
	 */
	protected abstract Sequence flatCloneWithoutAnnotation();

	/**
	 * Returns the length of the {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * @return the length of the {@link de.jstacs.data.sequences.Sequence}
	 */
	public abstract int getLength();

	/**
	 * Returns a {@link String} representation of the {@link de.jstacs.data.sequences.Sequence} (normally
	 * the {@link de.jstacs.data.sequences.Sequence} in its original {@link de.jstacs.data.alphabets.Alphabet}).
	 * 
	 * @return the {@link de.jstacs.data.sequences.Sequence} as {@link String}
	 * 
	 * @see Sequence#toString(String, int, int)
	 */
	@Override
	public String toString() {
		return toString( alphabetCon.getDelim(), 0, getLength() );
	}

	/**
	 * Returns a {@link String} representation of the {@link de.jstacs.data.sequences.Sequence} (normally
	 * the {@link de.jstacs.data.sequences.Sequence} in its original {@link de.jstacs.data.alphabets.Alphabet}) beginning at
	 * position <code>start</code> with a default delimiter as separator.
	 * 
	 * @param start
	 *            the start index (inclusive)
	 * 
	 * @return the {@link de.jstacs.data.sequences.Sequence} as {@link String}
	 * 
	 * @see Sequence#toString(String, int, int)
	 */
	public String toString( int start ) {
		return toString( alphabetCon.getDelim(), start, getLength() );
	}

	/**
	 * Returns a {@link String} representation of the {@link de.jstacs.data.sequences.Sequence} (normally
	 * the {@link de.jstacs.data.sequences.Sequence} in its original {@link de.jstacs.data.alphabets.Alphabet}) between
	 * <code>start</code> and <code>end</code> with a default delimiter as
	 * separator.
	 * 
	 * @param start
	 *            the start index (inclusive)
	 * @param end
	 *            the end index (exclusive)
	 * 
	 * @return the {@link de.jstacs.data.sequences.Sequence} as {@link String}
	 * 
	 * @see Sequence#toString(String, int, int)
	 */
	public String toString( int start, int end ) {
		return toString( alphabetCon.getDelim(), start, end );
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo( Sequence<T> s ) {
		//old return this.toString().compareTo( s.toString() );
		int c = alphabetCon.compareTo( s.alphabetCon );
		if( c == 0 ) {
			int l = getLength(), seqL = s.getLength();
			if( l == seqL ) {
				int i = 0;
				T t1 = getEmptyContainer();
				T t2 = s.getEmptyContainer();
				
				while( i < l ) {
					fillContainer( t1, i );
					s.fillContainer( t2, i );
					c = compareTo( t1, t2 );
					if( c != 0 ) {
						break;
					}
					i++;
				}
				return i < l ? c : 0;
			} else {
				return l - seqL;
			}
		} else {
			return c;
		}
	}
	
	/**
	 * This method compares to container and is used in {@link #compareTo(Sequence)}.
	 * 
	 * @param t1 the first container
	 * @param t2 the second container
	 * 
	 * @return zero if arguments are equal
	 * 
	 * @see #getEmptyContainer()
	 * @see #fillContainer(Object, int)
	 * @see Comparable#compareTo(java.lang.Object)
	 */
	protected abstract int compareTo( T t1, T t2 );

	/**
	 * This method converts a continuous value at position <code>pos</code> of
	 * the {@link de.jstacs.data.sequences.Sequence} into a discrete one.
	 * 
	 * @param pos
	 *            the position of the {@link de.jstacs.data.sequences.Sequence}
	 * @param content
	 *            the value at this position
	 * 
	 * @return the discrete value for this position
	 * 
	 * @see AlphabetContainer#toDiscrete(int, double)
	 */
	protected int toDiscrete( int pos, double content ) {
		return alphabetCon.toDiscrete( pos, content );
	}

	/**
	 * Returns a {@link String} representation of the {@link de.jstacs.data.sequences.Sequence} (normally
	 * the {@link de.jstacs.data.sequences.Sequence} in its original alphabet) between <code>start</code>
	 * and <code>end</code> with <code>delim</code> as separator.
	 * 
	 * @param delim
	 *            the delimiter/separator
	 * @param start
	 *            the start index (inclusive)
	 * @param end
	 *            the end index (exclusive)
	 * 
	 * @return the {@link de.jstacs.data.sequences.Sequence} as {@link String}
	 * 
	 * @see #getEmptyRepresentation()
	 * @see #addToRepresentation(Object, int, String)
	 * @see #getStringRepresentation(Object)
	 */
	public String toString( String delim, int start, int end ) {
		Object representation = getEmptyRepresentation();
		for( int i = start; i < end; i++ ) {
			addToRepresentation( representation, i, i<end-1 ? delim : "" );
		}
		return getStringRepresentation( representation );
	}
	
	/**
	 * Returns an empty representation which is used to create the {@link String} representation of this instance in the method {@link #toString(String, int, int)}.
	 * 
	 * @return an empty representation which is used to create the {@link String} representation
	 * 
	 * @see #toString(String, int, int)
	 */
	protected abstract Object getEmptyRepresentation();
	
	/**
	 * This method adds the information of one position to the representation using the specified delimiter  
	 * 
	 * @param representation the representation
	 * @param pos the position
	 * @param delim the delimiter separating the information for different positions
	 * 
	 * @see #getEmptyRepresentation()
	 * @see #toString(String, int, int)
	 */
	protected abstract void addToRepresentation( Object representation, int pos, String delim );
	
	/**
	 * This method creates a String representation from the given representation.
	 * 
	 * @param representation the representation instance (which should be created by {@link #getEmptyContainer()} and filled by {@link #addToRepresentation(Object, int, String)})
	 * 
	 * @return a String representation
	 * 
	 * @see #getEmptyRepresentation()
	 * @see #addToRepresentation(Object, int, String)
	 * @see #toString(String, int, int)
	 */
	protected abstract String getStringRepresentation( Object representation );

	/**
	 * Creates a {@link de.jstacs.data.sequences.Sequence} from a {@link String} based on the given
	 * {@link AlphabetContainer} using the standard delimiter for this
	 * {@link AlphabetContainer}.
	 * 
	 * @param con
	 *            the {@link AlphabetContainer}
	 * @param sequence
	 *            the {@link String} containing the {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @return a new {@link de.jstacs.data.sequences.Sequence} instance
	 * 
	 * @throws WrongAlphabetException
	 *             if <code>sequence</code> is not defined over <code>con</code>
	 * @throws IllegalArgumentException
	 *             if the delimiter is empty and the {@link AlphabetContainer}
	 *             is not discrete
	 * 
	 * @see Sequence#create(AlphabetContainer, String, String)
	 */
	public static Sequence create( AlphabetContainer con, String sequence ) throws WrongAlphabetException, IllegalArgumentException {
		return create( con, sequence, con.getDelim() );
	}

	/**
	 * Creates a {@link de.jstacs.data.sequences.Sequence} from a {@link String} based on the given
	 * {@link AlphabetContainer} using the given delimiter <code>delim</code>.
	 * 
	 * @param con
	 *            the {@link AlphabetContainer}
	 * @param sequence
	 *            the {@link String} containing the {@link de.jstacs.data.sequences.Sequence}
	 * @param delim
	 *            the given delimiter
	 * 
	 * @return a new {@link de.jstacs.data.sequences.Sequence} instance
	 * 
	 * @throws WrongAlphabetException
	 *             if <code>sequence</code> is not defined over <code>con</code>
	 * @throws IllegalArgumentException
	 *             if the delimiter is empty and the {@link AlphabetContainer}
	 *             is not discrete
	 * 
	 * @see Sequence#create(AlphabetContainer, SequenceAnnotation[], String,
	 *      String)
	 */
	public static Sequence create( AlphabetContainer con, String sequence, String delim ) throws WrongAlphabetException,
			IllegalArgumentException {
		return create( con, null, sequence, delim );
	}

	/**
	 * Creates a {@link de.jstacs.data.sequences.Sequence} from a {@link String} based on the given
	 * {@link AlphabetContainer} using the given delimiter <code>delim</code>
	 * and some <code>annotation</code> for the {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * @param con
	 *            the {@link AlphabetContainer}
	 * @param annotation
	 *            the annotation for the {@link de.jstacs.data.sequences.Sequence}
	 * @param sequence
	 *            the {@link String} containing the {@link de.jstacs.data.sequences.Sequence}
	 * @param delim
	 *            the given delimiter
	 * 
	 * @return a new {@link de.jstacs.data.sequences.Sequence} instance
	 * 
	 * @throws WrongAlphabetException
	 *             if <code>sequence</code> is not defined over <code>con</code>
	 * @throws IllegalArgumentException
	 *             if the delimiter is empty and the {@link AlphabetContainer}
	 *             is not discrete
	 */
	public static Sequence create( AlphabetContainer con, SequenceAnnotation[] annotation, String sequence, String delim ) throws WrongAlphabetException,
			IllegalArgumentException {
		try {
			if( con.isDiscrete() ) {
				// create pure discrete sequence
				int l = (int)con.getMaximalAlphabetLength();
				if( l <= Byte.MAX_VALUE ) {
					return new ByteSequence( con, annotation, sequence, delim );
				} else if( l <= Short.MAX_VALUE ) {
					return new ShortSequence( con, annotation, sequence, delim );
				} else if( l <= Integer.MAX_VALUE ) {
					return new IntSequence( con, annotation, sequence, delim );
				} else {
					throw new WrongAlphabetException( "Could not encode. Too many symbols." );
				}
			} else {
				// create hybrid or pure continuous sequence
				return new ArbitrarySequence( con, annotation, sequence, delim );
			}
		} catch ( WrongSequenceTypeException e ) {
			RuntimeException doesNotHappen = new RuntimeException( e.getMessage() );
			doesNotHappen.setStackTrace( e.getStackTrace() );
			throw doesNotHappen;
		}
	}

	/**
	 * This method returns a new instance of {@link de.jstacs.data.sequences.Sequence} containing the
	 * reverse current {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * <br>
	 * 
	 * So invoking this method, for instance, on the sequence &quot;TAATA&quot;
	 * returns &quot;ATAAT&quot;.
	 * 
	 * @return the reverse {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @throws OperationNotSupportedException
	 *             if the current {@link de.jstacs.data.sequences.Sequence} is based on an
	 *             {@link AlphabetContainer} that is not simple
	 * 
	 * @see Sequence#reverse(int, int)
	 */
	public final Sequence reverse() throws OperationNotSupportedException {
		return reverse( 0, getLength() );
	}

	/**
	 * This method returns a new instance of {@link de.jstacs.data.sequences.Sequence} containing a part
	 * of the reverse current {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * @param start
	 *            the start position (inclusive) in the original
	 *            {@link de.jstacs.data.sequences.Sequence}
	 * @param end
	 *            the end position (exclusive) in the original {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @return the reverse {@link de.jstacs.data.sequences.Sequence} of the part
	 * 
	 * @throws OperationNotSupportedException
	 *             if the current {@link de.jstacs.data.sequences.Sequence} is based on an
	 *             {@link AlphabetContainer} that is not simple
	 */
	public Sequence reverse( int start, int end ) throws OperationNotSupportedException {
		if( alphabetCon.isSimple() ) {
			int i = 0, j = end;
			try {
				if( this instanceof SimpleDiscreteSequence || (this instanceof SubSequence && ((SubSequence)this).content instanceof SimpleDiscreteSequence ) ) {
					int[] erg = new int[j - start];
					for( j--; j >= start; j--, i++ ) {
						erg[i] = discreteVal( j );
					}
					return new IntSequence( alphabetCon, erg );
				} else {
					double[] erg = new double[j - start];
					for( j--; j >= start; j--, i++ ) {
						erg[i] = continuousVal( j );
					}
					return new ArbitrarySequence( alphabetCon, erg );
				}
			} catch ( Exception e ) {
				// the current sequence is defined on the correct alphabet => so no WrongAlphabetException can occur
				RuntimeException doesNotHappen = new RuntimeException( e.getMessage() );
				doesNotHappen.setStackTrace( e.getStackTrace() );
				throw doesNotHappen;
			}
		} else {
			throw new OperationNotSupportedException( "The sequence has to be simple." );
		}
	}

	/**
	 * This method returns a new instance of {@link de.jstacs.data.sequences.Sequence} containing the
	 * complementary current {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * <br>
	 * 
	 * So invoking this method, for instance, on the sequence &quot;TAATA&quot;
	 * with an {@link AlphabetContainer} on
	 * {@link de.jstacs.data.alphabets.DNAAlphabet} returns &quot;ATTAT&quot;.
	 * 
	 * @return the complementary {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @throws OperationNotSupportedException
	 *             if the current {@link de.jstacs.data.sequences.Sequence} is not based on a
	 *             {@link ComplementableDiscreteAlphabet}
	 * 
	 * @see ComplementableDiscreteAlphabet
	 * @see Sequence#complement(int, int)
	 */
	public Sequence complement() throws OperationNotSupportedException {
		return complement( 0, getLength() );
	}

	/**
	 * This method returns a new instance of {@link de.jstacs.data.sequences.Sequence} containing the
	 * reverse complementary current {@link de.jstacs.data.sequences.Sequence}. For more details see the
	 * methods {@link #reverse()} and {@link #complement()}.
	 * 
	 * @return the reverse complementary {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @throws OperationNotSupportedException
	 *             if the current {@link de.jstacs.data.sequences.Sequence} is not discrete and simple
	 *             (not based on a {@link ComplementableDiscreteAlphabet})
	 * 
	 * @see Sequence#reverse()
	 * @see Sequence#complement()
	 * @see Sequence#reverseComplement(int, int)
	 * @see ComplementableDiscreteAlphabet
	 */
	public Sequence reverseComplement() throws OperationNotSupportedException {
		return reverseComplement( 0, getLength() );
	}
	
	/**
	 * This method returns a new instance of {@link de.jstacs.data.sequences.Sequence} containing a part
	 * of the complementary current {@link de.jstacs.data.sequences.Sequence}.
	 * 
	 * <br>
	 * 
	 * So invoking this method, for instance, on the sequence &quot;TAATA&quot;
	 * with an {@link AlphabetContainer} on
	 * {@link de.jstacs.data.alphabets.DNAAlphabet} returns &quot;ATTAT&quot;.
	 * 
	 * @param start
	 *            the start position (inclusive) in the original
	 *            {@link de.jstacs.data.sequences.Sequence}
	 * @param end
	 *            the end position (exclusive) in the original {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @return the complementary {@link de.jstacs.data.sequences.Sequence} of the part
	 * 
	 * @throws OperationNotSupportedException
	 *             if the current {@link de.jstacs.data.sequences.Sequence} is not based on a
	 *             {@link ComplementableDiscreteAlphabet}
	 * 
	 * @see ComplementableDiscreteAlphabet
	 */
	public Sequence complement( int start, int end ) throws OperationNotSupportedException {
		if( alphabetCon.isReverseComplementable() ) {
			ComplementableDiscreteAlphabet cda = (ComplementableDiscreteAlphabet)alphabetCon.getAlphabetAt( 0 );
			try {
				if( cda.length() > 127 ) {
					int[] erg = new int[end - start];
					for( int i = 0, j = start; j < end; i++, j++ ) {
						erg[i] = cda.getComplementaryCode( discreteVal( j ) );
					}

					return new IntSequence( alphabetCon, erg );
				} else {
					byte[] erg = new byte[end - start];
					for( int i = 0, j = start; j < end; i++, j++ ) {
						erg[i] = (byte)cda.getComplementaryCode( discreteVal( j ) );
					}

					return new ByteSequence( alphabetCon, erg );
				}
			} catch ( Exception e ) {
				// the current sequence is defined on the correct alphabet => so no WrongAlphabetException can occur
				RuntimeException doesNotHappen = new RuntimeException( e.getMessage() );
				doesNotHappen.setStackTrace( e.getStackTrace() );
				throw doesNotHappen;
			}
		} else {
			throw new OperationNotSupportedException( "The alphabet of sequence has to be complementable." );
		}
	}

	/**
	 * This method returns a new instance of {@link de.jstacs.data.sequences.Sequence} containing a
	 * reverse part of the complementary current {@link de.jstacs.data.sequences.Sequence}. For more
	 * details see the methods {@link #reverse()} and {@link #complement()}.
	 * 
	 * @param start
	 *            the start position (inclusive) in the original
	 *            {@link de.jstacs.data.sequences.Sequence}
	 * @param end
	 *            the end position (exclusive) in the original {@link de.jstacs.data.sequences.Sequence}
	 * 
	 * @return the reverse complementary {@link de.jstacs.data.sequences.Sequence} of the part
	 * 
	 * @throws OperationNotSupportedException
	 *             if the current {@link de.jstacs.data.sequences.Sequence} is not discrete and simple
	 *             ((not based on a {@link ComplementableDiscreteAlphabet})
	 * 
	 * @see Sequence#reverse()
	 * @see Sequence#complement()
	 * @see ComplementableDiscreteAlphabet
	 */
	public Sequence reverseComplement( int start, int end ) throws OperationNotSupportedException {
		if( rc != null && start == 0 && end == getLength() ) {
			return rc;
		} else if( alphabetCon.isReverseComplementable() ) {
			SimpleDiscreteSequence revComp;
			try {
				ComplementableDiscreteAlphabet cda = (ComplementableDiscreteAlphabet)alphabetCon.getAlphabetAt( 0 );
				int i = 0, j = end;
				if( cda.length() > 127 ) {
					int[] erg = new int[end - start];
					for( j--; j >= start; j--, i++ ) {
						erg[i] = cda.getComplementaryCode( discreteVal( j ) );
					}
			
					revComp = new IntSequence( alphabetCon, erg );
			
				} else {
					byte[] erg = new byte[end - start];
					for( j--; j >= start; j--, i++ ) {
						erg[i] = (byte)cda.getComplementaryCode( discreteVal( j ) );
					}
			
					revComp = new ByteSequence( alphabetCon, erg );
				}
				if( start == 0 && end == getLength() ) {
					rc = (Sequence<T>) revComp;
					rc.rc = this;
				}
				return revComp;
			} catch ( Exception e ) {
				// the current sequence is defined on the correct alphabet => so no WrongAlphabetException can occur
				RuntimeException doesNotHappen = new RuntimeException( e.getMessage() );
				doesNotHappen.setStackTrace( e.getStackTrace() );
				throw doesNotHappen;
			}
		} else {
			throw new OperationNotSupportedException( "The alphabet of sequence has to be reverse-complementable." );
		}
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return hashCode==null ? computeHashCode() : hashCode;
	}
	
	private int computeHashCode() {
		int len = getLength(), hashCode = 0;
		for( int i = 0; i < len; i++ ) {
			hashCode = 31 * hashCode + hashCodeForPos( i );
		}
		if( annotation != null ) {
			hashCode = 31*hashCode + Arrays.hashCode( annotation );
		}
		this.hashCode = hashCode;
		return hashCode;
	}
	
	/**
	 * This method is used in {@link #hashCode()} and the hash code for one specific position.
	 *   
	 * @param pos the position
	 * 
	 * @return the hash code for the position
	 */
	protected abstract int hashCodeForPos( int pos );
	
	/**
	 * This method returns the Hamming distance between the current {@link de.jstacs.data.sequences.Sequence} and <code>seq</code>.
	 * If the sequence have different length -1 is returned.
	 * 
	 * @param seq the sequence to be compared
	 * 
	 * @return the Hamming distance
	 * 
	 * @throws WrongAlphabetException it the sequences have different {@link AlphabetContainer}
	 */
	public int getHammingDistance( Sequence seq ) throws WrongAlphabetException {
		if( !alphabetCon.checkConsistency( seq.alphabetCon ) ) {
			throw new WrongAlphabetException();
		}
		if( getLength() != seq.getLength() ) {
			return -1;
		} else {
			return getHammingDistance( 0, seq, Integer.MAX_VALUE );
		}
	}
	
	private int getHammingDistance( int startPos, Sequence shortSequence, int hMax ) {
		int l = shortSequence.getLength();
		int h = 0;
		for( int p = 0; p < l && h < hMax; p++ ) {
			h += (continuousVal( startPos+p ) == shortSequence.continuousVal( p ) ? 0 : 1);
		}
		return h;		
	}
	
	//TODO more efficient;
	/**
	 * This method allows to answer the question whether there is a
	 * similar pattern find a match with a given maximal number of mismatches. 
	 * 
	 * @param maxHammingDistance the maximal Hamming distance
	 * @param shortSequence the short sequence
	 * 
	 * @return <code>true</code> if a match with maximal Hamming distance smaller than <code>maxHammingDistance</code> exists, otherwise <code>false</code>
	 * 
	 * @throws WrongAlphabetException if the sequence have different {@link AlphabetContainer}
	 */
	public boolean matches( int maxHammingDistance, Sequence shortSequence ) throws WrongAlphabetException {
		if( !alphabetCon.checkConsistency( shortSequence.alphabetCon ) ) {
			throw new WrongAlphabetException();
		}
		int l = getLength() - shortSequence.getLength();
		if( l < 0 ) {
			return false;
		} else {
			for( int h, p = 0; p <= l; p++ ) {
				 h = getHammingDistance(p, shortSequence,maxHammingDistance+1);
				 if( h <= maxHammingDistance ) {
					 return true;
				 }
			}
			return false;
		}
	}
	
	/**
	 * The method returns <code>true</code> if the sequence is multidimensional, otherwise <code<false</code>.
	 * 
	 * @return <code>true</code> if the sequence is multidimensional, otherwise <code<false</code>
	 */
	public abstract boolean isMultiDimensional();
	
	/**
	 * The method returns a container that can be used for accessing the symbols for each position.
	 * This is especially of interest for multidimensional sequences.
	 * 
	 * @return a container that can be used for accessing the symbols for each position
	 * 
	 * @see #fillContainer(Object, int)
	 * @see #isMultiDimensional()
	 */
	public abstract T getEmptyContainer();
	
	/**
	 * The method fills the content of a specific position in to the container.
	 * This is especially of interest for multidimensional sequences.
	 * 
	 * @param container the container which is used for filling the content.
	 * @param pos the position
	 * 
	 * @see #getEmptyContainer()
	 * @see #isMultiDimensional()
	 */
	public abstract void fillContainer( T container, int pos );

	/**
	 * This is the main class for subsequences, composite sequences, ... . All these
	 * sequences are defined on an existing {@link de.jstacs.data.sequences.Sequence}. After creating an
	 * {@link Sequence.RecursiveSequence}, no {@link SequenceAnnotation} of the internally {@link de.jstacs.data.sequences.Sequence}
	 * is returned by {@link Sequence#getAnnotation()}, ...
	 * 
	 * @param <T> the type of each position
	 * 
	 * @author Jens Keilwagen
	 */
	public static abstract class RecursiveSequence<T> extends Sequence<T> {

		/**
		 * The internal sequence.
		 */
		protected Sequence<T> content;

		/**
		 * Creates a new {@link Sequence.RecursiveSequence} on the {@link de.jstacs.data.sequences.Sequence}
		 * <code>seq</code> with the {@link AlphabetContainer} <code>alphabet</code>
		 * and the annotation <code>annotation</code>.
		 * 
		 * @param alphabet
		 *            the {@link AlphabetContainer}
		 * @param annotation
		 *            the annotation of the {@link Sequence.RecursiveSequence}
		 * @param seq
		 *            the sequence
		 * 
		 * @see Sequence#Sequence(AlphabetContainer, SequenceAnnotation[])
		 */
		public RecursiveSequence( AlphabetContainer alphabet, SequenceAnnotation[] annotation, Sequence<T> seq ) {
			super( alphabet, annotation );
			content = seq;
		}

		/**
		 * Creates a new {@link Sequence.RecursiveSequence} on the {@link de.jstacs.data.sequences.Sequence}
		 * <code>seq</code> with the {@link AlphabetContainer} <code>alphabet</code>
		 * using the annotation of the given {@link de.jstacs.data.sequences.Sequence}.
		 * 
		 * @param alphabet
		 *            the {@link AlphabetContainer}
		 * @param seq
		 *            the sequence
		 */
		public RecursiveSequence( AlphabetContainer alphabet, Sequence<T> seq ) {
			this( alphabet, null, seq );
		}
		
		/* (non-Javadoc)
		 * @see de.jstacs.data.Sequence#continuousVal(int)
		 */
		@Override
		public double continuousVal( int pos ) {
			return content.continuousVal( getIndex( pos ) );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.data.Sequence#discreteVal(int)
		 */
		@Override
		public int discreteVal( int pos ) {
			return content.discreteVal( getIndex( pos ) );
		}

		/**
		 * Returns the index in the internal sequence.
		 * 
		 * @param pos
		 *            the index in the external sequence
		 * 
		 * @return the index in the internal sequence
		 */
		protected abstract int getIndex( int pos );
		
		/*
		 * (non-Javadoc)
		 * @see de.jstacs.data.Sequence#isMultiDimensional()
		 */
		public boolean isMultiDimensional()  {
			return content.isMultiDimensional();
		}

		/*
		 * (non-Javadoc)
		 * @see de.jstacs.data.Sequence#getEmptyContainer()
		 */
		public T getEmptyContainer() {
			return content.getEmptyContainer();
		}
		
		/*
		 * (non-Javadoc)
		 * @see de.jstacs.data.Sequence#fillContainer(java.lang.Object, int)
		 */
		public void fillContainer( T container, int pos ) {
			content.fillContainer( container, getIndex( pos ) );
		}
		
		public int compareTo( T t1, T t2 ) {
			return content.compareTo( t1, t2 );
		}
		
		protected Object getEmptyRepresentation() {
			return content.getEmptyRepresentation();
		}
		
		protected void addToRepresentation( Object representation, int pos, String delim ) {
			content.addToRepresentation( representation, getIndex( pos ), delim );
		}
		protected String getStringRepresentation( Object representation ) {
			return content.getStringRepresentation( representation );
		}
		
		protected int hashCodeForPos( int pos ) {
			return content.hashCodeForPos( getIndex( pos ) );
		}
	}
	
	/**
	 * The class handles composite {@link de.jstacs.data.sequences.Sequence}s. A
	 * {@link Sequence.CompositeSequence} consists of several (partial) {@link de.jstacs.data.sequences.Sequence}
	 * s. A biological example are promoters like in eukaryotes (-10 and -35
	 * box).
	 * 
	 * @param <T> the type of each position
	 * 
	 * @author Jens Keilwagen
	 */
	protected static class CompositeSequence<T> extends RecursiveSequence<T> {

		private static final long serialVersionUID = 1L;

		private int[] starts, lengths;

		private CompositeSequence( AlphabetContainer abc, SequenceAnnotation[] annotation, Sequence<T> seq, int[] starts, int[] lengths,
									boolean check ) {
			super( abc, annotation, seq );
			// TODO make this faster, no new object
			if( check && !abc.checkConsistency( seq.getAlphabetContainer().getCompositeContainer( starts, lengths ) ) ) {
				throw new IllegalArgumentException( "Wrong AlphabetContainer." );
			}
			if( starts.length != lengths.length ) {
				throw new IllegalArgumentException( "starts and lengths have to be from the same dimension" );
			}
			this.starts = starts.clone();
			this.lengths = lengths.clone();
			for( int i = 0; i < lengths.length; i++ ) {
				if( starts[i] < 0 || seq.getLength() <= starts[i] ) {
					throw new IllegalArgumentException( i + "-th start has to be from range [0," + ( seq.getLength() - 1 ) + "]" );
				}
				if( lengths[i] < 0 || seq.getLength() < starts[i] + lengths[i] ) {
					throw new IllegalArgumentException( i + "-th length has to be from range [0," + ( seq.getLength() - starts[i] ) + "]" );
				}
				this.starts[i] = starts[i];
				this.lengths[i] = lengths[i];
			}
		}

		/**
		 * This is a very efficient way to create a {@link Sequence.CompositeSequence}
		 * for {@link de.jstacs.data.sequences.Sequence}s with a simple {@link AlphabetContainer}.
		 * 
		 * @param seq
		 *            the original {@link de.jstacs.data.sequences.Sequence}
		 * @param starts
		 *            the start positions of the junks
		 * @param lengths
		 *            the length of each junk
		 */
		public CompositeSequence( Sequence seq, int[] starts, int[] lengths ) {
			this( seq.getAlphabetContainer().getCompositeContainer( starts, lengths ), null, seq, starts, lengths, false );
		}

		/**
		 * This constructor should be used if one wants to create a
		 * {@link de.jstacs.data.DataSet} of {@link Sequence.CompositeSequence}s. With this constructor
		 * you are enabled to create a {@link de.jstacs.data.DataSet} where every
		 * {@link de.jstacs.data.sequences.Sequence} has the same {@link AlphabetContainer} instance.
		 * 
		 * <br>
		 * <br>
		 * 
		 * Internally it is checked that the {@link AlphabetContainer} matches
		 * with the one of the subsequence.
		 * 
		 * @param abc
		 *            the new {@link AlphabetContainer}
		 * @param seq
		 *            the original {@link de.jstacs.data.sequences.Sequence}
		 * @param starts
		 *            the start positions of the junks
		 * @param lengths
		 *            the length of each junk
		 */
		public CompositeSequence( AlphabetContainer abc, Sequence<T> seq, int[] starts, int[] lengths ) {
			this( abc, null, seq, starts, lengths, true );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.data.sequences.RecursiveSequence#getIndex(int)
		 */
		@Override
		protected int getIndex( int pos ) {
			int sum = 0, i = 0;
			while( i < lengths.length && sum + lengths[i] <= pos ) {
				sum += lengths[i++];
			}
			return starts[i] + ( pos - sum );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.data.Sequence#getLength()
		 */
		@Override
		public int getLength() {
			int sum = 0;
			for( int i = 0; i < lengths.length; i++ ) {
				sum += lengths[i];
			}
			return sum;
		}

		/* (non-Javadoc)
		 * @see de.jstacs.data.Sequence#flatCloneWithoutAnnotation()
		 */
		@Override
		protected Sequence flatCloneWithoutAnnotation() {
			return new CompositeSequence( alphabetCon, null, content, starts, lengths, false );
		}
	}

	/**
	 * This class handles subsequences. {@link Sequence.SubSequence}s are often used to
	 * extract the {@link de.jstacs.data.sequences.Sequence} of a sliding window on a long
	 * {@link de.jstacs.data.sequences.Sequence}. The class is implemented in such a way that it avoids
	 * chains of {@link Sequence.SubSequence}s.
	 * 
	 * @param <T> the type of each position
	 * 
	 * @author Jens Keilwagen
	 */
	public static class SubSequence<T> extends RecursiveSequence<T> {

		private static final long serialVersionUID = 1L;

		private int start, length;

		private SubSequence( AlphabetContainer abc, SequenceAnnotation[] annotation, Sequence<T> seq, int start, int length, boolean check ) {
			super( abc, annotation, ( seq instanceof SubSequence ) ? ( (SubSequence)seq ).content : seq );
			// TODO make this faster, no new object
			if( check && !abc.checkConsistency( seq.getAlphabetContainer().getSubContainer( start, length ) ) ) {
				throw new IllegalArgumentException( "Wrong AlphabetContainer." );
			}
			int sl = seq.getLength();
			if( start < 0 || start > sl ) {
				throw new IllegalArgumentException( "Illegal start position: start=" + start + " not in [0," + sl + "]" );
			}
			if( length < 0 || start + length > sl ) {
				throw new IllegalArgumentException( "Illegal length: length=" + length + " not in [0," + ( sl - start ) + "]" );
			}
			if( seq instanceof SubSequence ) {
				this.start = ((SubSequence)seq).start + start;
			} else {
				this.start = start;
			}
			this.length = length;
		}

		/**
		 * This constructor should be used if one wants to create a
		 * {@link de.jstacs.data.DataSet} of {@link Sequence.SubSequence}s of defined length. With this
		 * constructor you are enabled to create a {@link de.jstacs.data.DataSet} where every
		 * {@link de.jstacs.data.sequences.Sequence} has the same {@link AlphabetContainer} instance.
		 * 
		 * <br>
		 * <br>
		 * 
		 * Internally it is checked that the {@link AlphabetContainer} matches
		 * with the one of the {@link Sequence.SubSequence}.
		 * 
		 * @param abc
		 *            the new {@link AlphabetContainer}
		 * @param seq
		 *            the original {@link de.jstacs.data.sequences.Sequence}
		 * @param start
		 *            the index of the start position
		 * @param length
		 *            the length of the new sequence
		 */
		public SubSequence( AlphabetContainer abc, Sequence seq, int start, int length ) {
			this( abc, null, seq, start, length, true );
		}

		/**
		 * This is a very efficient way to create a {@link Sequence.SubSequence} of
		 * defined length for {@link de.jstacs.data.sequences.Sequence}s with a simple
		 * {@link AlphabetContainer}.
		 * 
		 * @param seq
		 *            the original {@link de.jstacs.data.sequences.Sequence}
		 * @param start
		 *            the index of the start position
		 * @param length
		 *            the length of the new {@link de.jstacs.data.sequences.Sequence}
		 */
		public SubSequence( Sequence seq, int start, int length ) {
			this( seq.getAlphabetContainer().getSubContainer( start, length ), null, seq, start, length, false );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.data.Sequence#reverseComplement(int, int)
		 */
		@Override
		public Sequence reverseComplement( int start, int end ) throws OperationNotSupportedException {
			if( rc != null && start == 0 && end == getLength() ) {
				return rc;
			} else if( content.rc != null ) {
				Sequence revComp = new SubSequence( content.rc, content.rc.getLength() - this.start - end, end - start );
				if( start == 0 && end == getLength() ) {
					rc = revComp;
					rc.rc = this;
				}
				return revComp;
			}
			return super.reverseComplement( start, end );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.data.sequences.RecursiveSequence#getIndex(int)
		 */
		@Override
		protected final int getIndex( int pos ) {
			if( pos < 0 || pos >= length ) {
				throw new ArrayIndexOutOfBoundsException( pos );
			}
			return start + pos;
		}

		public Sequence<T> getOriginal(){
			return content;
		}
		
		public int getStart(){
			return start;
		}
		
		/* (non-Javadoc)
		 * @see de.jstacs.data.sequences.RecursiveSequence#discreteVal(int)
		 */
		@Override
		public final int discreteVal( int pos ) {
			if( pos < 0 || pos >= length ) {
				throw new ArrayIndexOutOfBoundsException( pos );
			}
			return content.discreteVal( start + pos );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.data.Sequence#getLength()
		 */
		@Override
		public int getLength() {
			return length;
		}

		/* (non-Javadoc)
		 * @see de.jstacs.data.Sequence#flatCloneWithoutAnnotation()
		 */
		@Override
		protected Sequence flatCloneWithoutAnnotation() {
			return new SubSequence( alphabetCon, null, content, start, length, false );
		}
	}
}