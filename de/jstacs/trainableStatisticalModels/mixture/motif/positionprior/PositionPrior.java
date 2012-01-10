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

package de.jstacs.trainableStatisticalModels.mixture.motif.positionprior;

import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.io.XMLParser;

/**
 * This is the main class for any position prior that can be used in a motif
 * discovery.
 * 
 * @author Jens Keilwagen
 */
public abstract class PositionPrior implements Storable, Cloneable {

	/**
	 * The length of the motif.
	 */
	protected int motifLength;

	/**
	 * This empty constructor creates an instance with motif length -1.
	 */
	protected PositionPrior() {
		motifLength = -1;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link PositionPrior} out of a {@link StringBuffer}.
	 * 
	 * @param rep
	 *            the {@link StringBuffer} containing the model
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} can not be parsed
	 */
	protected PositionPrior( StringBuffer rep ) throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag( rep, getInstanceName() );
		extractAdditionalInformation( xml );
		setMotifLength( XMLParser.extractObjectForTags( xml, "motifLength", int.class ) );
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	public PositionPrior clone() throws CloneNotSupportedException {
		return (PositionPrior)super.clone();
	}

	/**
	 * The logarithmic value of the prior for specified start positions of the
	 * part motifs.
	 * 
	 * @param seqLength
	 *            the length of the sequence
	 * @param starts
	 *            the start positions of the part motifs
	 * 
	 * @return the value of the prior for given start positions
	 * 
	 * @throws IllegalArgumentException
	 *             if something went wrong, e.g. the motif length is not set.
	 */
	public abstract double getLogPriorForPositions( int seqLength, int... starts ) throws IllegalArgumentException;

	/**
	 * Returns the length that is supported by this prior. If all lengths
	 * greater than the minimal length are supported, this method returns 0.
	 * 
	 * @return the length that is supported or 0 for all lengths
	 */
	public abstract int getLength();

	/**
	 * Returns the instance name.
	 * 
	 * @return the instance name
	 */
	public abstract String getInstanceName();

	/**
	 * Sets the length of the current motif.
	 * 
	 * @param motifLength
	 *            the motif length
	 * 
	 * @throws IllegalArgumentException
	 *             if the length is not positive or to big
	 */
	public void setMotifLength( int motifLength ) throws IllegalArgumentException {
		int l = getLength();
		if( motifLength > 0 && ( l == 0 || motifLength <= l ) ) {
			this.motifLength = motifLength;
		} else {
			throw new IllegalArgumentException( "The motif length has to be positive and at most as big as the sequence length." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 1000 );
		XMLParser.appendObjectWithTags( xml, motifLength, "motifLength" );
		xml.append( getAdditionalInformation() );
		XMLParser.addTags( xml, getInstanceName() );
		return xml;
	}

	/**
	 * This method returns a {@link StringBuffer} containing additional
	 * information for the XML representation.
	 * 
	 * @return a {@link StringBuffer} containing additional information for the
	 *         XML representation
	 */
	protected abstract StringBuffer getAdditionalInformation();

	/**
	 * This method extracts additional information from a {@link StringBuffer}.
	 * 
	 * @param xml
	 *            the {@link StringBuffer} containing the additional information
	 *            of the XML representation
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} is not parsable
	 */
	protected abstract void extractAdditionalInformation( StringBuffer xml ) throws NonParsableException;
}
