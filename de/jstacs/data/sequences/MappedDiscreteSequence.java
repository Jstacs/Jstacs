/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.data.sequences;

import de.jstacs.WrongAlphabetException;
import de.jstacs.data.Alphabet;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.DiscreteAlphabetMapping;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;


/**
 * This class allows to map a discrete {@link Sequence} to an new {@link Sequence} using some {@link DiscreteAlphabetMapping}s.
 * 
 * @author Jens Keilwagen
 */
public class MappedDiscreteSequence extends SimpleDiscreteSequence {

	/**
	 * The original {@link Sequence}.
	 */
	protected SimpleDiscreteSequence original;
	
	/**
	 * The original {@link AlphabetContainer}.
	 */
	protected AlphabetContainer originalAlphabetContainer;
	
	/**
	 * The {@link DiscreteAlphabetMapping}s used for mapping the original {@link Sequence}.
	 */
	protected DiscreteAlphabetMapping[] transformation;
	
	/**
	 * This method allows to create a new {@link AlphabetContainer} given an old {@link AlphabetContainer} and some {@link DiscreteAlphabetMapping}s.
	 * 
	 * @param original the original {@link AlphabetContainer}
	 * @param transformation the {@link DiscreteAlphabetMapping}s defining the mapping
	 * 
	 * @return a new {@link AlphabetContainer}
	 */
	public static final AlphabetContainer getNewAlphabetContainer( AlphabetContainer original, DiscreteAlphabetMapping... transformation ) {
		int n = original.getNumberOfAlphabets();
		if( n != 1 && n != transformation.length ) {
			throw new IllegalArgumentException( "Check number of transformations." );
		}
		int[] index = new int[transformation.length];
		for( int i = 0; i < index.length; i++ ) {
			index[i] = i;
		}
		Alphabet[] abc = new Alphabet[transformation.length];
		for( int i = 0; i < abc.length; i++ ) {
			abc[i] = transformation[i].getNewAlphabet();
		}
		return new AlphabetContainer(abc, index);
	}
	
	/**
	 * This method allows to create an empty {@link MappedDiscreteSequence}. That is a {@link MappedDiscreteSequence} that contains no original {@link Sequence}.
	 * 
	 * @param originalAlphabetContainer the original {@link AlphabetContainer}
	 * @param seqAn the {@link SequenceAnnotation}s for the {@link MappedDiscreteSequence}
	 * @param transformation the {@link DiscreteAlphabetMapping} defining the mapping
	 * 
	 * @throws WrongAlphabetException if the mapped {@link AlphabetContainer} is not discrete (should never happen)
	 */
	protected MappedDiscreteSequence( AlphabetContainer originalAlphabetContainer, SequenceAnnotation[] seqAn, DiscreteAlphabetMapping... transformation ) throws WrongAlphabetException {
		super( getNewAlphabetContainer( originalAlphabetContainer, transformation ), seqAn );
		this.originalAlphabetContainer = null;
		if( originalAlphabetContainer.getNumberOfAlphabets() > 1 ) {
			this.originalAlphabetContainer = originalAlphabetContainer;
		}
		this.transformation = transformation.clone();
	}
	
	/**
	 * This method allows to create a {@link MappedDiscreteSequence} from a given {@link Sequence} and some given {@link DiscreteAlphabetMapping}s.
	 * 
	 * @param original the original {@link Sequence}
	 * @param transformation the {@link DiscreteAlphabetMapping} defining the mapping
	 * 
	 * @throws WrongAlphabetException if the mapped {@link AlphabetContainer} is not discrete (should never happen)
	 */
	public MappedDiscreteSequence( SimpleDiscreteSequence original, DiscreteAlphabetMapping... transformation  ) throws WrongAlphabetException {
		this( original.getAlphabetContainer(), original.getAnnotation(), transformation );
		this.original = original;
	}
	
	@Override
	public int discreteVal( int pos ) {
		return transformation[getIndex( pos )].getNewDiscreteValue( original.discreteVal( pos ) );
	}

	@Override
	protected Sequence flatCloneWithoutAnnotation() {
		// TODO Auto-generated method stub: XXX JAN?
		return null;
	}

	@Override
	public int getLength() {
		return original.getLength();
	}
	
	/**
	 * This method returns the logarithm of the number of original {@link Sequence}s that yield the same mapped {@link Sequence}.
	 * 
	 * @return the logarithm of the number of original {@link Sequence}s that yield the same mapped {@link Sequence}
	 */
	public double getLogNumberOfPossibleOriginalSequences() {
		return getLogNumberOfPossibleOriginalSequences( 0, getLength() );
	}

	/**
	 * This method returns the logarithm of the number of original {@link Sequence}s that yield the same mapped {@link Sequence}.
	 *
	 * @param start the start position (inclusive)
	 * @param end the end position (exclusive)
	 * 
	 * @return the logarithm of the number of original {@link Sequence}s that yield the same mapped {@link Sequence}
	 */
	public double getLogNumberOfPossibleOriginalSequences( int start, int end ) {
		double res = 0;
		while( start < end ) {
			res += transformation[getIndex( start )].getLogNumberOfSimilarSymbols(original.discreteVal(start));
			start++;
		}
		return res;
	}
	
	private int getIndex( int pos ) {
		return originalAlphabetContainer == null ? pos : originalAlphabetContainer.getAlphabetIndexForPosition( pos );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getEmptyContainer()
	 */
	public final int[] getEmptyContainer() {
		return new int[1];
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#fillContainer(java.lang.Object, int)
	 */
	public final void fillContainer( int[] container, int pos ) {
		container[0] = discreteVal( pos );
	}
}
