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

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.FastDirichletMRGParams;

/**
 * This constraint can be used for any <b>m</b>aximum <b>e</b>ntropy
 * <b>m</b>odel (MEM) application.
 * 
 * @author Jens Keilwagen
 */
public class MEMConstraint extends InhConstraint {

	private double[] expLambda;

	private double[] lambda;

	private int[] corrected_positions;

	private static int[] isSorted( int[] pos ) throws IllegalArgumentException {
		int i = 1;
		while( i < pos.length && pos[i - 1] < pos[i] ) {
			i++;
		}
		if( i < pos.length ) {
			throw new IllegalArgumentException( "The position array is not unique." );
		}
		return pos;
	}

	/**
	 * Creates a {@link MEMConstraint} as part of a (whole) model.
	 * 
	 * @param pos
	 *            the used positions (have to be sorted)
	 * @param alphabetLength
	 *            an array containing the length of the the alphabet for each
	 *            position (of the whole model)
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>pos</code> is not sorted
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint#Constraint(int[], int)
	 */
	public MEMConstraint( int[] pos, int[] alphabetLength ) throws IllegalArgumentException {
		super( isSorted( pos ), alphabetLength );
		expLambda = new double[counts.length];
		Arrays.fill( expLambda, 1 );
		lambda = new double[counts.length];
		corrected_positions = usedPositions;
	}

	/**
	 * Creates a {@link MEMConstraint} as part of a model.
	 * 
	 * @param pos
	 *            the used positions (have to be sorted)
	 * @param alphabetLength
	 *            an array containing the length of the the alphabet for each
	 *            position (of the whole model)
	 * @param corrected_positions
	 *            an array containing the corrected positions if the maximum
	 *            entropy model is decomposed into parts
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>pos</code> is not sorted or
	 *             <code>pos.length != corrected_positions.length</code>
	 * 
	 * @see InhConstraint#InhConstraint(int[], int[])
	 */
	public MEMConstraint( int[] pos, int[] alphabetLength, int[] corrected_positions ) throws IllegalArgumentException {
		super( isSorted( pos ), alphabetLength );
		if( pos.length != corrected_positions.length ) {
			throw new IllegalArgumentException( "The length of pos and corrected_positions is not equal." );
		}
		expLambda = new double[counts.length];
		lambda = new double[counts.length];
		this.corrected_positions = corrected_positions.clone();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link MEMConstraint} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MEMConstraint} could not be reconstructed out
	 *             of the XML representation (the {@link StringBuffer} could not
	 *             be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see InhConstraint#InhConstraint(StringBuffer)
	 */
	public MEMConstraint( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.InhConstraint#clone()
	 */
	@Override
	public MEMConstraint clone() throws CloneNotSupportedException {
		MEMConstraint clone = (MEMConstraint)super.clone();
		clone.expLambda = expLambda.clone();
		clone.lambda = lambda.clone();
		if( corrected_positions == usedPositions ) {
			clone.corrected_positions = clone.usedPositions;
		} else {
			clone.corrected_positions = corrected_positions.clone();
		}
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint#estimate(double)
	 */
	@Override
	public void estimate( double ess ) {
		estimateUnConditional( 0, freq.length, ess / (double)freq.length, true );
	}
	
	/**
	 * Draws the parameters from a Dirichlet.
	 * 
	 * @param ess the equivalent sample size
	 */
	public void draw( double ess ) {
		FastDirichletMRGParams alpha = new FastDirichletMRGParams( ess / (double) freq.length );
		DirichletMRG.DEFAULT_INSTANCE.generate( expLambda, 0, expLambda.length, alpha );
		for( int i = 0; i < expLambda.length; i++ ) {
			lambda[i] = Math.log( expLambda[i] );
		}
	}

	/**
	 * Returns the value of the corrected position <code>index</code>.
	 * 
	 * @param index
	 *            the index of the position
	 * 
	 * @return the value of the corrected position
	 */
	public int getCorrectedPosition( int index ) {
		return corrected_positions[index];
	}

	/**
	 * Returns the exponential value of {@latex.inline $\\lambda_{index}$} at position
	 * <code>index</code>: {@latex.inline $\\exp(\\lambda_{index})$}.
	 * 
	 * @param index
	 *            the index
	 * 
	 * @return {@latex.inline $\\exp(\\lambda_{index})$}
	 */
	public double getExpLambda( int index ) {
		return expLambda[index];
	}

	/**
	 * Returns the value of {@latex.inline $\\lambda_{index}$}.
	 * 
	 * @param index
	 *            the index
	 * 
	 * @return {@latex.inline $\\lambda_{index}$}
	 */
	public double getLambda( int index ) {
		return lambda[index];
	}

	/**
	 * Multiplies the exponential value of {@latex.inline $\\lambda_{index}$} with the factor <code>val</code>:
	 * {@latex.inline $\\exp(\\lambda_{index}) \\cdot val$}. <br>
	 * (Additionally it adds the logarithmic value of <code>val</code> to the
	 * value of {@latex.inline $\\lambda_{index}$} at position <code>index</code>:
	 * {@latex.inline $\\lambda_{index} + \\log(val)$}.)
	 * 
	 * @param index
	 *            the index
	 * @param val
	 *            the factor/value
	 */
	public void multiplyExpLambdaWith( int index, double val ) {
		expLambda[index] *= val;
		lambda[index] += Math.log( val );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint#reset()
	 */
	@Override
	public void reset() {
		super.reset();
		Arrays.fill( expLambda, 1d );
		Arrays.fill( lambda, 0d );
	}

	/**
	 * Returns the index of the constraint that is satisfied by
	 * <code>sequence</code>.
	 * 
	 * @param sequence
	 *            the {@link SequenceIterator}
	 * 
	 * @return the index of the fulfilled constraint
	 */
	public int satisfiesSpecificConstraint( SequenceIterator sequence ) {
		int erg = 0, i = 0;
		for( ; i < corrected_positions.length; i++ ) {
			//erg += offset[i] * sequence.getValueAt( corrected_positions[i] );
			erg += offset[i] * sequence.seq[corrected_positions[i]];
		}
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint#getFreq(int)
	 */
	@Override
	public double getFreq( int index ) {
		return freq[index];
	}

	/**
	 * Sets the exponential value of {@latex.inline $\\exp(\\lambda_{index}) = val$}.<br>
	 * (Additionally it sets the value of {@latex.inline $\\lambda_{index}$} to the logarithmic value of <code>val</code>:
	 * {@latex.inline $\\lambda_{index} = \\log(val)$}.)
	 * 
	 * @param index
	 *            the index
	 * @param val
	 *            the value to be set
	 */
	public void setExpLambda( int index, double val ) {
		expLambda[index] = val;
		lambda[index] = Math.log( val );
	}

	/**
	 * Sets the value of {@latex.inline $\\lambda_{index} = val$}.<br>
	 * (Additionally it sets the exponential value of {@latex.inline $\\lambda_{index} = val$}
	 * to the exponential value of <code>val</code>:
	 * {@latex.inline $\\exp(\\lambda_{index}) = \\exp(val)$}.
	 * 
	 * @param index
	 *            the index
	 * @param val
	 *            the value
	 */
	public void setLambda( int index, double val ) {
		expLambda[index] = Math.exp( val );
		lambda[index] = val;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint#toString()
	 */
	@Override
	public String toString() {
		String erg = "" + usedPositions[0];
		for( int i = 1; i < usedPositions.length; i++ ) {
			erg += ", " + usedPositions[i];
		}
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.InhConstraint#appendAdditionalInfo(java.lang.StringBuffer)
	 */
	@Override
	protected void appendAdditionalInfo( StringBuffer xml ) {
		super.appendAdditionalInfo( xml );
		XMLParser.appendObjectWithTags( xml, lambda, "lambda" );
		if( corrected_positions != usedPositions ) {
			StringBuffer b = new StringBuffer( 500 );
			XMLParser.appendObjectWithTags( b, corrected_positions, "corrected_positions" );
			XMLParser.addTags( b, "corrected" );
			xml.append( b );
		}
	}

	private static String XML_TAG = "MEMConstraint";

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.Constraint#getXMLTag()
	 */
	@Override
	protected String getXMLTag() {
		return XML_TAG;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.InhConstraint#extractAdditionalInfo(java.lang.StringBuffer)
	 */
	@Override
	protected void extractAdditionalInfo( StringBuffer xml ) throws NonParsableException {
		super.extractAdditionalInfo( xml );
		lambda = XMLParser.extractObjectForTags( xml, "lambda", double[].class );
		expLambda = new double[lambda.length];
		for( int i = 0; i < lambda.length; i++ ) {
			expLambda[i] = Math.exp( lambda[i] );
		}
		StringBuffer corrected = XMLParser.extractForTag( xml, "corrected" );
		if( corrected == null ) {
			corrected_positions = usedPositions;
		} else {
			corrected_positions = XMLParser.extractObjectForTags( corrected, "corrected_positions", int[].class );
		}
	}
	
	public int comparePosition( int offset, MEMConstraint constr ) {
		int matching = 0;
		for( int u1 = 0, u2 = 0; u1 < usedPositions.length; u1++ ) {
			while( u2 < constr.usedPositions.length && (constr.usedPositions[u2]-offset) < usedPositions[u1] ) {
				u2++;
			}
			if( u2 < constr.usedPositions.length && (constr.usedPositions[u2]-offset) == usedPositions[u1] ) {
				matching++;
			}
		}
		return matching;
	}

	public void addParameters( int offset, IntList list, MEMConstraint[] constraint, double[] params, int[] start ) {
		int u, n, idx;
		HashSet<Integer> hash = new HashSet<Integer>();
		for( u = 0; u < usedPositions.length; u++ ) {
			hash.add( usedPositions[u] );
		}
		for( n = 0; n < list.length(); n++ ) {
			idx = list.get(n);
			for( u = 0; u < constraint[idx].usedPositions.length; u++ ) {
				hash.add( constraint[idx].usedPositions[u]-offset );
			}
		}
		int[] pos = new int[hash.size()];
		Iterator<Integer> it = hash.iterator();
		u = 0;
		while( it.hasNext() ) {
			pos[u] = it.next();
			u++;
		}
		Arrays.sort( pos );
		Arrays.fill( expLambda, Double.NEGATIVE_INFINITY );
		
		//assume all alphabets have same length
		int alphLength = (lambda.length / this.offset[0])-1;
		int[] assignment = new int[pos.length+1];
		int index, i;
		double p;
		while( assignment[pos.length]==0 ) {
			//compute the "marginalized" parameter
			p = 0;
			for( n = 0; n < list.length(); n++ ) {
				idx = list.get(n);
				index = 0;
				for( i = 0, u = 0; u < pos.length && i < constraint[idx].usedPositions.length; u++ ) {
					if( pos[u] == constraint[idx].usedPositions[i]-offset ) {
						index += constraint[idx].offset[i] * assignment[u];
						i++;
					}
				}
				p += params[start[idx]+index];
			}
			
			//compute the index for this constraint
			index = 0;
			for( i = 0, u = 0; u < pos.length && i < usedPositions.length; u++ ) {
				if( pos[u] == usedPositions[i] ) {
					index += this.offset[i] * assignment[u];
					i++;
				}
			}
			
			//System.out.println( Arrays.toString(assignment) + "\t" + index + "\t" + expLambda[index] + "\t" + p + "\t" + Normalisation.getLogSum( expLambda[index], p ) );			
			
			//add
			expLambda[index] = Normalisation.getLogSum( expLambda[index], p );
			
			//next
			u = 0;
			while( u < pos.length && assignment[u] == alphLength ) {
				assignment[u]=0;
				u++;
			}
			assignment[u]++;
		}
		
		for( u = 0; u < lambda.length; u++ ) {
			lambda[u] += expLambda[u];
			expLambda[u] = Math.exp( lambda[u] );
		}
		//System.out.println(Arrays.toString(lambda));
	}
}