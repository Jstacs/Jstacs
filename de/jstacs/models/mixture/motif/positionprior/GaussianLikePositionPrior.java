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

package de.jstacs.models.mixture.motif.positionprior;

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * This class implements a gaussian like discrete truncated prior.
 * 
 * @author Jens Keilwagen
 */
public class GaussianLikePositionPrior extends PositionPrior {

	private int length;

	private double max, prec, z;

	/**
	 * This constructor creates an instance with given sequence length, maximal
	 * value and sigma for the Gaussian.
	 * 
	 * @param length
	 *            the sequence length
	 * @param max
	 *            the maximal position
	 * @param sigma
	 *            the standard deviation (sigma) of the Gaussian
	 */
	public GaussianLikePositionPrior( int length, double max, double sigma ) {
		super();
		if( length < 1 ) {
			throw new IllegalArgumentException( "The length has to be positive." );
		}
		this.length = length;
		if( max < 0 ) {
			throw new IllegalArgumentException( "The maximum has to be non-negative." );
		}
		this.max = max;
		if( sigma <= 0 ) {
			throw new IllegalArgumentException( "The sigma has to be non-negative." );
		}
		this.prec = 1d / ( 2d * sigma * sigma );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link GaussianLikePositionPrior} out of a
	 * {@link StringBuffer}.
	 * 
	 * @param rep
	 *            the {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} can not be parsed
	 */
	public GaussianLikePositionPrior( StringBuffer rep ) throws NonParsableException {
		super( rep );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.motif.positionprior.PositionPrior#setMotifLength(int)
	 */
	@Override
	public void setMotifLength( int motifLength ) {
		super.setMotifLength( motifLength );
		z = 0;
		for( int m = length - motifLength, i = 0; i <= m; i++ ) {
			z += Math.exp( -prec * ( max - i ) * ( max - i ) );
		}
		z = Math.log( z );
	}

	/**
	 * Returns only the important part and leaving the logarithm of the
	 * normalization constant out.
	 */
	@Override
	public double getLogPriorForPositions( int seqLength, int... starts ) throws IllegalArgumentException {
		if( seqLength == length && starts[0] <= length - motifLength ) {
			return -prec * ( max - starts[0] ) * ( max - starts[0] ) - z;
		} else {
			throw new IllegalArgumentException( "This sequence length could not be modeled." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.motif.positionprior.PositionPrior#getLength()
	 */
	@Override
	public int getLength() {
		return length;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.motif.positionprior.PositionPrior#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return getClass().getSimpleName();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.motif.positionprior.PositionPrior#getAdditionalInformation()
	 */
	@Override
	protected StringBuffer getAdditionalInformation() {
		StringBuffer xml = new StringBuffer( 1000 );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		XMLParser.appendObjectWithTags( xml, max, "max" );
		XMLParser.appendObjectWithTags( xml, prec, "prec" );
		return xml;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.motif.positionprior.PositionPrior#extractAdditionalInformation(java.lang.StringBuffer)
	 */
	@Override
	protected void extractAdditionalInformation( StringBuffer xml ) throws NonParsableException {
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		max = XMLParser.extractObjectForTags( xml, "max", double.class );
		prec = XMLParser.extractObjectForTags( xml, "prec", double.class );
	}

}
