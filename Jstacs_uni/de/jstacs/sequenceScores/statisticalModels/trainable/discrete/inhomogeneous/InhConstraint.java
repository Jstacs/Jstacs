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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint;

/**
 * This class is the superclass for all inhomogeneous constraints.
 * 
 * @author Jens Keilwagen
 */
public abstract class InhConstraint extends Constraint implements Cloneable {

	/**
	 * This array is used to find the start indices of the conditional
	 * distributions.
	 */
	protected int[] offset;

	/**
	 * Computes the product of some selected values.
	 * 
	 * @param selected
	 *            the indices of the selected values
	 * @param values
	 *            the values to be computed
	 * 
	 * @return the product of the selected values
	 */
	private static int product( int[] selected, int[] values ) {
		int erg = 1, i = 0;
		while( i < selected.length ) {
			erg *= values[selected[i++]];
		}
		return erg;
	}

	/**
	 * Creates a new {@link InhConstraint} instance.
	 * 
	 * @param pos
	 *            the positions
	 * @param alphabetLength
	 *            the length of each alphabet (not only the used positions)
	 * 
	 * @see Constraint#Constraint(int[], int)
	 */
	protected InhConstraint( int[] pos, int[] alphabetLength ) {
		super( pos, product( pos, alphabetLength ) );
		int i = pos.length;
		offset = new int[i--];
		offset[i--] = 1;
		for( ; i >= 0; i-- ) {
			offset[i] = offset[i + 1] * alphabetLength[pos[i + 1]];
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link InhConstraint} instance out of its XML
	 * representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link InhConstraint} instance could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see Constraint#Constraint(StringBuffer)
	 */
	protected InhConstraint( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint#extractAdditionalInfo(java.lang.StringBuffer)
	 */
	@Override
	protected void extractAdditionalInfo( StringBuffer xml ) throws NonParsableException {
		offset = XMLParser.extractObjectForTags( xml, "offset", int[].class );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint#clone()
	 */
	@Override
	public InhConstraint clone() throws CloneNotSupportedException {
		InhConstraint clone = (InhConstraint)super.clone();
		clone.offset = (int[])offset.clone();
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint#satisfiesSpecificConstraint(de.jstacs.data.Sequence, int)
	 */
	@Override
	public int satisfiesSpecificConstraint( Sequence s, int start ) {
		int erg = 0, counter = 0;
		for( ; counter < usedPositions.length; counter++ ) {
			erg += offset[counter] * s.discreteVal( start + usedPositions[counter] );
		}
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint#appendAdditionalInfo(java.lang.StringBuffer)
	 */
	@Override
	protected void appendAdditionalInfo( StringBuffer xml ) {
		XMLParser.appendObjectWithTags( xml, offset, "offset" );
	}
	
	public String getDescription( AlphabetContainer con, int i ) {
		String res = null, s;
		DiscreteAlphabet d;
		for( int j = 0; j < offset.length; j++ ) {
			d = (DiscreteAlphabet) con.getAlphabetAt( usedPositions[j] );
			s = "X_" + usedPositions[j] + "=" + d.getSymbolAt( i / offset[j] );
			if( res == null ) {
				res = s;
			} else {
				res = s + ", " + res;
			}
			i = i % offset[j];
		}
		return "P(" + res + ")"; 
	}
}