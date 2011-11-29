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

package de.jstacs.models.discrete.inhomogeneous;

import java.io.OutputStream;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.data.Sequence;
import de.jstacs.models.discrete.DGMParameterSet;
import de.jstacs.models.discrete.DiscreteGraphicalModel;
import de.jstacs.models.discrete.inhomogeneous.parameters.IDGMParameterSet;
import de.jstacs.utils.SafeOutputStream;

/**
 * This class is the main class for all inhomogeneous <b>d</b>iscrete
 * <b>g</b>raphical <b>m</b>odels ({@link InhomogeneousDGM}).
 * 
 * @author Jens Keilwagen
 * 
 * @see IDGMParameterSet
 */
public abstract class InhomogeneousDGM extends DiscreteGraphicalModel {

	/**
	 * The lengths of the alphabets.
	 */
	int[] alphabetLength;

	/**
	 * This stream is used for comments, computation steps/results or any other
	 * kind of output during the training, ... etc.
	 */
	protected SafeOutputStream sostream;

	/**
	 * The default {@link OutputStream}.
	 */
	protected static final OutputStream DEFAULT_STREAM = System.out;

	/**
	 * Creates a new {@link InhomogeneousDGM} from a given
	 * {@link IDGMParameterSet}.
	 * 
	 * @param params
	 *            the given parameter set
	 * 
	 * @throws CloneNotSupportedException
	 *             if the parameter set could not be cloned
	 * @throws IllegalArgumentException
	 *             if the parameter set is not instantiated
	 * @throws NonParsableException
	 *             if the parameter set is not parsable
	 * 
	 * @see DiscreteGraphicalModel#DiscreteGraphicalModel(DGMParameterSet)
	 */
	public InhomogeneousDGM( IDGMParameterSet params ) throws CloneNotSupportedException, IllegalArgumentException, NonParsableException {
		super( params );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link InhomogeneousDGM} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link InhomogeneousDGM} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see DiscreteGraphicalModel#DiscreteGraphicalModel(StringBuffer)
	 */
	public InhomogeneousDGM( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#clone()
	 */
	@Override
	public InhomogeneousDGM clone() throws CloneNotSupportedException {
		InhomogeneousDGM i = (InhomogeneousDGM)super.clone();
		i.alphabetLength = alphabetLength.clone();
		i.setOutputStream( sostream.doesNothing() ? null : SafeOutputStream.DEFAULT_STREAM );
		return i;
	}

	/**
	 * Returns a {@link String} representation of the underlying graph.
	 * 
	 * @return a {@link String} representation of the underlying graph
	 * 
	 * @throws NotTrainedException
	 *             if the structure is not set, this can only be the case if the
	 *             model is not trained
	 */
	public abstract String getStructure() throws NotTrainedException;

	/**
	 * Sets the {@link OutputStream} for the model. This stream will sometimes
	 * be used to write some information about the
	 * training/progress/computation... etc. to the screen, a file ... etc.
	 * 
	 * @param stream
	 *            the {@link OutputStream}
	 * 
	 * @see SafeOutputStream
	 */
	public void setOutputStream( OutputStream stream ) {
		sostream = SafeOutputStream.getSafeOutputStream( stream );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#check(de.jstacs.data.Sequence, int, int)
	 */
	@Override
	protected void check( Sequence sequence, int startpos, int endpos ) throws NotTrainedException, IllegalArgumentException {
		super.check( sequence, startpos, endpos );
		if( startpos > endpos || endpos >= sequence.getLength() ) {
			throw new IllegalArgumentException( "This endposition is impossible. Try: startposition <= endposition < sequence.length" );
		} else if( endpos - startpos + 1 != length ) {
			throw new IllegalArgumentException( "This sequence has not length " + length + "." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#set(de.jstacs.models.discrete.DGMParameterSet, boolean)
	 */
	@Override
	protected void set( DGMParameterSet parameter, boolean trained ) throws CloneNotSupportedException, NonParsableException {
		super.set( parameter, trained );
		alphabetLength = new int[length];
		for( int i = 0; i < length; i++ ) {
			alphabetLength[i] = (int)alphabets.getAlphabetLengthAt( i );
		}
		setOutputStream( SafeOutputStream.DEFAULT_STREAM );
	}
}
