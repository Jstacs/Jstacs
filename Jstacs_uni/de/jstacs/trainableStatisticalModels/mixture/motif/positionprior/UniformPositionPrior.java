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

/**
 * This prior implements a uniform distribution for the start position.
 * 
 * @author Jens Keilwagen
 */
public class UniformPositionPrior extends PositionPrior {

	/**
	 * This empty constructor creates an instance with motif length -1.
	 */
	public UniformPositionPrior() {
		super();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link UniformPositionPrior} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} can not be parsed
	 */
	public UniformPositionPrior( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.mixture.motif.positionprior.PositionPrior#getLogPriorForPositions(int, int[])
	 */
	@Override
	public double getLogPriorForPositions( int seqLength, int... starts ) throws IllegalArgumentException {
		if( motifLength > 0 ) {
			return -Math.log( seqLength - motifLength + 1 );
		} else {
			throw new IllegalArgumentException( "The motif length has to be set." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.mixture.motif.positionprior.PositionPrior#getLength()
	 */
	@Override
	public final int getLength() {
		return 0;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.mixture.motif.positionprior.PositionPrior#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return getClass().getSimpleName();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.mixture.motif.positionprior.PositionPrior#getAdditionalInformation()
	 */
	@Override
	protected StringBuffer getAdditionalInformation() {
		return new StringBuffer();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.mixture.motif.positionprior.PositionPrior#extractAdditionalInformation(java.lang.StringBuffer)
	 */
	@Override
	protected void extractAdditionalInformation( StringBuffer xml ) {}
}
