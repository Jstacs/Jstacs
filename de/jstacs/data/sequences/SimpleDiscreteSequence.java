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

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;

/**
 * This is the main class for any discrete sequence.
 * 
 * @author Jens Keilwagen
 */
public abstract class SimpleDiscreteSequence extends Sequence<int[]> {

	/**
	 * This constructor creates a new {@link SimpleDiscreteSequence} with the
	 * {@link AlphabetContainer} <code>container</code> and the annotation
	 * <code>annotation</code> but without the content. The content has to be
	 * set by the constructor of the subclass.
	 * 
	 * @param container
	 *            the {@link AlphabetContainer} of the sequence
	 * @param annotation
	 *            the annotation of the sequence
	 * 
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} is not discrete
	 * 
	 * @see de.jstacs.data.sequences.Sequence#Sequence(AlphabetContainer, SequenceAnnotation[])
	 */
	public SimpleDiscreteSequence( AlphabetContainer container, SequenceAnnotation[] annotation ) throws WrongAlphabetException {
		super( container, annotation );
		if( !container.isDiscrete() ) {
			throw new WrongAlphabetException( "The alphabet is not discrete." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#continuousVal(int)
	 */
	@Override
	public final double continuousVal( int pos ) {
		return discreteVal( pos );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#isMultiDimensional()
	 */
	public boolean isMultiDimensional()  {
		return false;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getEmptyContainer()
	 */
	public int[] getEmptyContainer() {
		return new int[1];
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#fillContainer(java.lang.Object, int)
	 */
	public void fillContainer( int[] container, int pos ) {
		container[0] = discreteVal( pos );
	}
	
	protected Object getEmptyRepresentation() {
		return new StringBuffer();
	}
	
	protected void addToRepresentation( Object representation, int pos, String delim ) {
		((StringBuffer)representation).append( alphabetCon.getSymbol( pos, discreteVal( pos ) ) + delim );
	}
	
	protected String getStringRepresentation( Object representation ) {
		return representation.toString();
	}
	
	protected int hashCodeForPos( int pos ) {
		return discreteVal( pos );
	}
	
	public int compareTo( int[] t1, int[] t2 ) {
		if( t1.length == t2.length ) {
			for( int i = 0; i < t1.length; i++ ) {
				if( t1[i] != t2[i] ) {
					return t1[i]-t2[i];
				}
			}
			return 0;
		} else {
			return t1.length - t2.length;
		}
	}
}
