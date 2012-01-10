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

package de.jstacs.trainableStatisticalModels.discrete;

import java.util.Arrays;

import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;
import de.jstacs.io.XMLParser;

/**
 * The main class for all constraints on models. This class centralizes all
 * specific constraints (e.g. P(X<sub>1</sub>=A,X<sub>2</sub>=A)) from a
 * partitioning constraint (e.g. all specific constraints for some given positions).
 * 
 * @author Jens Keilwagen
 */
public abstract class Constraint implements Storable {

	/**
	 * The counts for each specific constraint.
	 */
	protected double[] counts;

	/**
	 * The frequencies estimated from the counts.
	 */
	protected double[] freq;

	/**
	 * The used positions.
	 */
	protected int[] usedPositions;

	/**
	 * The main constructor. Creates a new {@link Constraint} and checks that
	 * each position is used maximally once.
	 * 
	 * @param usedPositions
	 *            the used positions (will be cloned), have to be non-negative
	 * @param n
	 *            the number of specific constraints
	 * @throws IllegalArgumentException if there are some fundamental errors in the arguments for the positions
	 */
	protected Constraint( int[] usedPositions, int n ) throws IllegalArgumentException {
		// is unique??
		int max = -1, i = 0;
		for( ; i < usedPositions.length; i++ ) {
			if( usedPositions[i] < 0 ) {
				throw new IllegalArgumentException( "The positions have to be non-negative." );
			}
			if( max < usedPositions[i] ) {
				max = usedPositions[i];
			}
		}
		boolean[] used = new boolean[max + 1];
		for( i = 0; i < usedPositions.length; i++ ) {
			if( used[usedPositions[i]] ) {
				throw new IllegalArgumentException( "Each position can be used only once, corrupted at " + usedPositions[i] + "." );
			} else {
				used[usedPositions[i]] = true;
			}
		}
		this.usedPositions = usedPositions.clone();
		counts = new double[n];
		freq = new double[n];
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link Constraint} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link Constraint} could not be reconstructed out of
	 *             the XML representation (the {@link StringBuffer} could not be
	 *             parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see Constraint#extractAdditionalInfo(StringBuffer)
	 */
	protected Constraint( StringBuffer xml ) throws NonParsableException {
		StringBuffer erg = XMLParser.extractForTag( xml, getXMLTag() );
		usedPositions = XMLParser.extractObjectForTags( erg, "usedPositions", int[].class );
		counts = XMLParser.extractObjectForTags( erg, "counts", double[].class );
		freq = XMLParser.extractObjectForTags( erg, "freq", double[].class );
		extractAdditionalInfo( erg );
	}

	/**
	 * Adds the given <code>weight</code> to the count with index
	 * <code>index</code>.
	 * 
	 * @param index
	 *            the index
	 * @param weight
	 *            the weight
	 */
	public final void add( int index, double weight ) {
		counts[index] += weight;
	}

	/**
	 * This method determines the specific constraint that is fulfilled by the
	 * {@link Sequence} <code>seq</code> and adds the <code>weight</code> to the
	 * specific counter.
	 * 
	 * @param seq
	 *            the sequence
	 * @param start
	 *            the start position
	 * @param weight
	 *            the weight for the sequence
	 * 
	 * @see Constraint#satisfiesSpecificConstraint(Sequence, int)
	 * @see Constraint#add(int, double)
	 */
	public void add( Sequence seq, int start, double weight ) {
		add( satisfiesSpecificConstraint( seq, start ), weight );
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected Constraint clone() throws CloneNotSupportedException {
		Constraint p = (Constraint)super.clone();
		p.counts = counts.clone();
		p.freq = freq.clone();
		p.usedPositions = usedPositions.clone();
		return p;
	}

	/**
	 * Estimates the (smoothed) relative frequencies using the ess
	 * (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize).
	 * 
	 * @param ess
	 *            the ess
	 */
	public abstract void estimate( double ess );

	/**
	 * Estimates unconditionally.
	 * 
	 * @param start
	 *            the start index
	 * @param end
	 *            the end index
	 * @param pc
	 *            the pseudocount for each parameter
	 * @param exceptionWhenNoData
	 *            indicates if an (runtime) exception is thrown if no data was
	 *            available to estimate the parameters
	 */
	protected void estimateUnConditional( int start, int end, double pc, boolean exceptionWhenNoData ) {
		double sum = 0;
		for( int index = start; index < end; index++ ) {
			sum += counts[index];
		}
		sum += ( end - start ) * pc;
		if( sum <= 0 ) {
			if( exceptionWhenNoData ) {
				throw new IllegalArgumentException( "A marginal distribution returned an illegal value (" + sum + ")." );
			} else {
				Arrays.fill( freq, start, end, 1d / (double)( end - start ) );
			}
		} else {
			while( start < end ) {
				freq[start] = ( counts[start++] + pc ) / sum;
			}
		}
	}

	/**
	 * This method appends additional information that is not stored in the base
	 * class to the {@link StringBuffer}.
	 * 
	 * @param xml
	 *            the {@link StringBuffer} that is used for appending additional
	 *            information
	 */
	protected abstract void appendAdditionalInfo( StringBuffer xml );

	/**
	 * Returns the current count with index <code>index</code>.
	 * 
	 * @param index
	 *            the index
	 * 
	 * @return the current count
	 */
	public double getCount( int index ) {
		return counts[index];
	}

	/**
	 * Returns the current frequency with index <code>index</code>.
	 * 
	 * @param index
	 *            the index
	 * 
	 * @return the current frequency
	 */
	public double getFreq( int index ) {
		return freq[index];
	}

	/**
	 * This method determines the specific constraint that is fulfilled by the
	 * {@link Sequence} <code>seq</code> beginning at position
	 * <code>start</code>.
	 * 
	 * @param seq
	 *            the {@link Sequence}
	 * @param start
	 *            the start position
	 * 
	 * @return the according frequency
	 * 
	 * @see Constraint#satisfiesSpecificConstraint(Sequence, int)
	 * @see Constraint#getFreq(int)
	 */
	public double getFreq( Sequence seq, int start ) {
		return getFreq( satisfiesSpecificConstraint( seq, start ) );
	}

	/**
	 * Returns the marginal order, i.e. the number of used random variables.
	 * 
	 * @return the marginal order
	 */
	public int getMarginalOrder() {
		return usedPositions.length;
	}

	/**
	 * Returns the number of specific constraints.
	 * 
	 * @return the number of specific constraints
	 */
	public int getNumberOfSpecificConstraints() {
		return counts.length;
	}

	/**
	 * Returns the position with index <code>index</code>.
	 * 
	 * @param index
	 *            the index
	 * 
	 * @return the position with index <code>index</code>
	 */
	public int getPosition( int index ) {
		return usedPositions[index];
	}

	/**
	 * Returns a clone of the array of used positions.
	 * 
	 * @return a clone of the array of used positions
	 */
	public int[] getPositions() {
		return usedPositions.clone();
	}

	/**
	 * Returns the XML tag that is used for the class to en- or decode.
	 * 
	 * @return the XML tag that is used for the class to en- or decode
	 */
	protected abstract String getXMLTag();

	/**
	 * This method resets all member variables that are in some way counters,
	 * frequencies, ...
	 */
	public void reset() {
		resetCounts();
		Arrays.fill( freq, 0 );
	}

	/**
	 * This method returns the index of the specific constraint that is
	 * fulfilled by the {@link Sequence} <code>seq</code> beginning at position
	 * <code>start</code>.
	 * 
	 * @param seq
	 *            the sequence
	 * @param start
	 *            the start position
	 * 
	 * @return the index of the fulfilled, specific constraint
	 */
	public abstract int satisfiesSpecificConstraint( Sequence seq, int start );

	/**
	 * This method parses additional information from the {@link StringBuffer}
	 * that is not parsed in the base class.
	 * 
	 * @param xml
	 *            the {@link StringBuffer} to be parsed
	 * 
	 * @throws NonParsableException
	 *             if something with the parsing went wrong
	 */
	protected abstract void extractAdditionalInfo( StringBuffer xml ) throws NonParsableException;

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public abstract String toString();

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer erg = new StringBuffer( 10000 );
		XMLParser.appendObjectWithTags( erg, usedPositions, "usedPositions" );
		XMLParser.appendObjectWithTags( erg, counts, "counts" );
		XMLParser.appendObjectWithTags( erg, freq, "freq" );
		appendAdditionalInfo( erg );
		XMLParser.addTags( erg, getXMLTag() );
		return erg;
	}
	
	/**
	 * Returns the decoded symbol for the encoded symbol <code>i</code>.
	 * 
	 * @param con the {@link AlphabetContainer}
	 * @param i the encoded symbol
	 * 
	 * @return the decoded symbol for the encoded symbol <code>i</code>
	 */
	public abstract String getDescription( AlphabetContainer con, int i );
	
	/**
	 * Returns an information about the stored frequencies.
	 * 
	 * @param con the {@link AlphabetContainer}
	 * 
	 * @return an information about the stored frequencies
	 */
	public String getFreqInfo( AlphabetContainer con ) {
		StringBuffer sb = new StringBuffer();
		int l = (int) con.getAlphabetLengthAt( usedPositions[usedPositions.length-1] );
		for( int i = 0; i < freq.length; ) {
			for( int j = 0; j < l; j++, i++ ) {
				sb.append( getDescription( con, i ) + " = " + freq[i] + "\t" );	
			}
			sb.append( "\n" );
		}
		
		return sb.toString();
	}
	
	/**
	 * This method allows to reset the internal field {@link Constraint#counts} that is used for estimating the probabilities.
	 * 
	 * @see #add(int, double)
	 * @see #add(Sequence, int, double)
	 * @see #estimate(double)
	 */
	public void resetCounts() {
		Arrays.fill( counts, 0 );
	}
}
