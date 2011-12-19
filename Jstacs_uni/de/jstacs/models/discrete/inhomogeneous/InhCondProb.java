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

import javax.naming.OperationNotSupportedException;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;

/**
 * This class handles (conditional) probabilities of sequences for
 * inhomogeneous models.
 * 
 * @author Jens Keilwagen
 */
public class InhCondProb extends InhConstraint {

	private boolean conditional;

	private double[] lnFreq;

	/**
	 * Creates a new {@link InhCondProb} instance.
	 * 
	 * @param pos
	 *            the position
	 * @param alphabetLength
	 *            the length of each alphabet (not only the used position)
	 * 
	 * @see InhCondProb#InhCondProb(int[], int[], boolean)
	 */
	public InhCondProb( int pos, int... alphabetLength ) {
		this( new int[]{ pos }, alphabetLength, false );
	}

	/**
	 * Creates a new {@link InhCondProb} instance.
	 * 
	 * @param pos
	 *            the positions
	 * @param alphabetLength
	 *            the length of each alphabet (not only the used positions)
	 * @param cond
	 *            indicates if the instance has to use conditional probabilities
	 * 
	 * @see InhConstraint#InhConstraint(int[], int[])
	 */
	public InhCondProb( int[] pos, int[] alphabetLength, boolean cond ) {
		super( pos, alphabetLength );
		conditional = pos.length > 1 && cond;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link InhCondProb} instance out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link InhCondProb} instance could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see InhConstraint#InhConstraint(StringBuffer)
	 */
	public InhCondProb( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.inhomogeneous.InhConstraint#clone()
	 */
	@Override
	public InhCondProb clone() throws CloneNotSupportedException {
		InhCondProb clone = (InhCondProb)super.clone();
		if( lnFreq != null ) {
			clone.lnFreq = null;
			clone.lnFreq( 0, freq.length );
		}
		return clone;
	}

	/**
	 * Draws the parameters from a Dirichlet distribution using the counts and
	 * the given <code>ess</code> (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * as hyperparameters.
	 * 
	 * @param ess
	 *            the given ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * 
	 * @see InhCondProb#drawUnConditional(int, int, double)
	 */
	public void drawParameters( double ess ) {
		double pc = ess / (double)getNumberOfSpecificConstraints();
		if( !conditional || usedPositions.length == 1 ) {
			drawUnConditional( 0, freq.length, pc );
		} else {
			// conditional
			int l = offset[offset.length - 2];
			for( int counter1 = 0; counter1 < freq.length; counter1 += l ) {
				drawUnConditional( counter1, counter1 + l, pc );
			}
		}
		lnFreq( 0, freq.length );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.Constraint#estimate(double)
	 */
	@Override
	public void estimate( double ess ) {
		double pc = ess / (double)getNumberOfSpecificConstraints();
		if( !conditional || usedPositions.length == 1 ) {
			estimateUnConditional( 0, freq.length, pc, false );
		} else {
			// conditional
			int l = offset[offset.length - 2];
			for( int counter1 = 0; counter1 < freq.length; counter1 += l ) {
				estimateUnConditional( counter1, counter1 + l, pc, false );
			}
		}
	}

	/**
	 * Estimates the unconditional frequencies using the ess (<b>e</b>quivalent
	 * <b>s</b>ample <b>s</b>ize).
	 * 
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * @param all
	 *            the sum of all weights used to fill the counts
	 */
	public void estimateUnConditional( double ess, double all ) {
		conditional = false;
		all = ess + all;
		double pc = ess / (double)getNumberOfSpecificConstraints();
		for( int i = 0; i < freq.length; i++ ) {
			freq[i] = ( counts[i] + pc ) / all;
		}
		lnFreq( 0, freq.length );
	}

	/**
	 * Returns the logarithm of the relative frequency (=probability) at
	 * position <code>index</code> in the distribution.
	 * 
	 * @param index
	 *            the index of the entry in the distribution
	 * 
	 * @return the logarithm of the relative frequency (=probability)
	 */
	public double getLnFreq( int index ) {
		return lnFreq[index];
	}

	/**
	 * Returns the logarithm of the relative frequency (=probability) with the
	 * position in the distribution given by the index of the specific
	 * constraint that is fulfilled by the {@link Sequence} <code>s</code>
	 * beginning at <code>start</code>.
	 * 
	 * @param s
	 *            the sequence
	 * @param start
	 *            the index of the start position
	 * 
	 * @return the logarithm of the relative frequency (=probability)
	 * 
	 * @see InhCondProb#satisfiesSpecificConstraint(Sequence, int)
	 */
	public double getLnFreq( Sequence s, int start ) {
		return lnFreq[satisfiesSpecificConstraint( s, start )];
	}

	/**
	 * This method is used to create random sequences.
	 * 
	 * @param content
	 *            the content of the random sequence as far as it is known
	 * @param p
	 *            a random number in (0,1)
	 * 
	 * @throws OperationNotSupportedException
	 *             if this instance models a joint probability for more than one
	 *             position (shall be implemented in the future)
	 * 
	 * @see de.jstacs.StatisticalModel#emitDataSet(int, int[])
	 */
	public void getOutput( byte[] content, double p ) throws OperationNotSupportedException {
		int off = 0;
		byte i = 0;
		if( conditional ) {
			for( ; i < usedPositions.length - 1; i++ ) {
				off += content[usedPositions[i]] * offset[i];
			}
		} else if( usedPositions.length > 1 ) {
			throw new OperationNotSupportedException();
		}

		i = 0;
		while( p > freq[off] ) {
			p -= freq[off++];
			i++;
		}
		content[usedPositions[usedPositions.length - 1]] = i;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.Constraint#toString()
	 */
	@Override
	public String toString() {
		String erg = "";
		int i = 1, l = usedPositions.length - 1;
		if( l > 0 ) {
			erg += usedPositions[0];
			while( i < l ) {
				erg += ", " + usedPositions[i++];
			}
			if( conditional ) {
				erg += " -> ";
			}
		}
		return erg + usedPositions[l];
	}

	/**
	 * This method draws the parameters for a part of this constraint. It is
	 * used to draw from a distribution with fixed context.
	 * 
	 * @param start
	 *            the start index
	 * @param end
	 *            the end index
	 * @param pc
	 *            the pseudocount/hyperparameter
	 */
	protected void drawUnConditional( int start, int end, double pc ) {
		double[] alpha = new double[end - start];
		for( int index = 0; index < alpha.length; index++ ) {
			alpha[index] = counts[start + index] + pc;
		}
		DirichletMRG.DEFAULT_INSTANCE.generate( freq, start, alpha.length, new DirichletMRGParams( alpha ) );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.Constraint#estimateUnConditional(int, int, double, boolean)
	 */
	@Override
	protected void estimateUnConditional( int start, int end, double pc, boolean exceptionWhenNoData ) {
		super.estimateUnConditional( start, end, pc, exceptionWhenNoData );
		lnFreq( start, end );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.inhomogeneous.InhConstraint#appendAdditionalInfo(java.lang.StringBuffer)
	 */
	@Override
	protected void appendAdditionalInfo( StringBuffer xml ) {
		super.appendAdditionalInfo( xml );
		XMLParser.appendObjectWithTags( xml, conditional, "conditional" );
	}

	private static final String XML_TAG = "InhCondProb";

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.Constraint#getXMLTag()
	 */
	@Override
	protected String getXMLTag() {
		return XML_TAG;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.inhomogeneous.InhConstraint#extractAdditionalInfo(java.lang.StringBuffer)
	 */
	@Override
	protected void extractAdditionalInfo( StringBuffer xml ) throws NonParsableException {
		super.extractAdditionalInfo( xml );
		conditional = XMLParser.extractObjectForTags( xml, "conditional", boolean.class );
		lnFreq( 0, freq.length );
	}

	/**
	 * Computes the logarithm for all frequencies.
	 */
	private void lnFreq( int start, int end ) {
		if( lnFreq == null ) {
			lnFreq = new double[freq.length];
		}
		for( int i = start; i < end; i++ ) {
			lnFreq[i] = Math.log( freq[i] );
		}
	}

	/**
	 * This method is used to restore the values of a Gibbs Sampling run. The
	 * parameters are encoded in {@link String}s. Each {@link String} is one
	 * frequency. The index <code>start</code> is used to begin at a specific
	 * position in the array.
	 * 
	 * @param array
	 *            the array of {@link String} chunks to be parsed
	 * @param start
	 *            the start index
	 * 
	 * @throws IllegalArgumentException
	 *             if something is wrong with the frequencies
	 */
	public void setFreqs( String[] array, int start ) throws IllegalArgumentException {
		int l;
		if( !conditional || usedPositions.length == 1 ) {
			l = freq.length;
		} else {
			l = offset[offset.length - 2];
		}

		double p = 0;
		for( int k = 0, i = 0; i < freq.length; i++ ) {
			freq[i] = Double.parseDouble( array[start++] );
			if( freq[i] < 0 ) {
				throw new IllegalArgumentException( "this is no proper frequency" );
			}
			p += freq[i];
			k++;
			if( k == l ) {
				if( Math.abs( 1 - p ) > 1E-7 ) {
					throw new IllegalArgumentException( "The frequencies does not sum to 1." );
				}
				k = 0;
				p = 0;
			}
		}
		lnFreq( 0, freq.length );
	}
	
	@Override
	public String getDescription( AlphabetContainer con, int i ) {
		String res = super.getDescription( con, i );
		if( conditional ) {
			res = res.replaceFirst( ", ", " | " );
		}
		return res;
	}
}