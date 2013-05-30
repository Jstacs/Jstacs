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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous;

import java.text.NumberFormat;
import java.util.AbstractList;
import java.util.Arrays;

import de.jstacs.Storable;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.StatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.ConstraintManager;
import de.jstacs.utils.SafeOutputStream;

/**
 * This class represents a maximum entropy model. It should not be used without the MEManager.
 * 
 * @author Jens Keilwagen
 */
public final class MEM implements Storable, Cloneable
{
	private static final String XML_TAG = "MEM";

	/**
	 * The specific constraints of this MEM.
	 */
	protected MEMConstraint[] constraints;

	// normalizing constant
	private double lnZ, Z;

	// the positions in the sequence that are used in this MEM
	private int[] indices;

	// the positions in the sequence that are used in this MEM
	private int[] alphLen;

	// the conditional fixed positions (or null)
	private int[][] cond;

	private MEM( int[] indices, int[] allAlphLen )
	{
		this.indices = indices.clone();
		this.alphLen = new int[indices.length];
		for( int i = 0; i < indices.length; i++ )
		{
			this.alphLen[i] = allAlphLen[indices[i]];
		}
	}

	/**
	 * The main constructor of a MEM.
	 * 
	 * @param constr
	 *            the constraints
	 * @param allAlphLen
	 *            the alphabet length of all alphabets
	 * @param indices
	 *            the indices of the positions to be used
	 */
	public MEM( AbstractList<int[]> constr, int[] allAlphLen, int[] indices )
	{
		this( indices, allAlphLen );
		this.constraints = ConstraintManager.createConstraints( constr, allAlphLen, indices );
	}

	/**
	 * The main constructor of a MEM.
	 * 
	 * @param constraint
	 *            the constraint
	 * @param allAlphLen
	 *            the alphabet length of all alphabets
	 * @param cond
	 *            the conditional fixed positions
	 */
	public MEM( int[] constraint, int[] allAlphLen, int[][] cond )
	{
		this( constraint, allAlphLen );
		int i = 0, j;
		int[] corrected = new int[constraint.length];
		while( i < constraint.length )
		{
			corrected[i] = i++;
		}
		this.constraints = new MEMConstraint[]{ new MEMConstraint( constraint, allAlphLen, corrected ) };
		if( cond != null )
		{
			i = 0;
			corrected = new int[allAlphLen.length];
			while( i < indices.length )
			{
				corrected[indices[i]] = i++;
			}
			this.cond = new int[cond.length][];
			for( i = 0; i < cond.length; i++ )
			{
				this.cond[i] = new int[cond[i].length];
				for( j = 0; j < cond[i].length; j++ )
				{
					this.cond[i][j] = corrected[cond[i][j]];
				}
			}
		}
	}

	/**
	 * The constructor for the <code>Storable</code> interface.
	 * 
	 * @param representation
	 *            the representation
	 * 
	 * @throws NonParsableException
	 *             if the the StringBuffer could not be parsed.
	 */
	public MEM( StringBuffer representation ) throws NonParsableException
	{
		StringBuffer xml = XMLParser.extractForTag( representation, XML_TAG );
		indices = XMLParser.extractObjectForTags( xml, "used position", int[].class );
		alphLen = XMLParser.extractObjectForTags( xml, "alphLen", int[].class );
		Z = XMLParser.extractObjectForTags( xml, "normalization constant", double.class );
		lnZ = Math.log( Z ); 
		constraints = XMLParser.extractObjectForTags( xml, "constraints", MEMConstraint[].class );
		xml = XMLParser.extractForTag( xml, "conditions" );
		if( xml == null )
		{
			cond = null;
		}
		else
		{
			cond = XMLParser.extractObjectForTags( xml, "cond", int[][].class );
		}
	}

	public MEM clone() throws CloneNotSupportedException
	{
		MEM clone = (MEM) super.clone();
		int i = 0;
		clone.constraints = ArrayHandler.clone( constraints );
		clone.indices = indices.clone();
		clone.alphLen = alphLen.clone();
		if( cond != null )
		{
			clone.cond = new int[cond.length][];
			for( i = 0; i < constraints.length; i++ )
			{
				clone.cond[i] = cond[i].clone();
			}
		}
		return clone;
	}

	/**
	 * Returns the score for the <code>sequence<code> beginning at <code>start</code>.
	 * 
	 * @param seq the sequence
	 * @param start the start position
	 * 
	 * @return the score for the <code>sequence<code> beginning at <code>start</code>
	 */
	public double getScoreFor( Sequence seq, int start )
	{
		double erg = constraints[0].getExpLambda( constraints[0].satisfiesSpecificConstraint( seq, start ) );
		for( int i = 1; i < constraints.length; i++ )
		{
			erg *= constraints[i].getExpLambda( constraints[i].satisfiesSpecificConstraint( seq, start ) );
		}
		return erg / Z;
	}
	
	/**
	 * Returns the logarithmic score for the <code>sequence<code> beginning at <code>start</code>.
	 * 
	 * @param seq the sequence
	 * @param start the start position
	 * 
	 * @return the logarithmic score for the <code>sequence<code> beginning at <code>start</code>
	 */
	public double getLogScoreFor( Sequence seq, int start )
	{
		double erg = constraints[0].getLambda( constraints[0].satisfiesSpecificConstraint( seq, start ) );
		for( int i = 1; i < constraints.length; i++ )
		{
			erg += constraints[i].getLambda( constraints[i].satisfiesSpecificConstraint( seq, start ) );
		}
		return erg - lnZ;
	}
	
	/**
	 * This method compute the prior for the current parameter ignoring some constants.
	 * 
	 * @param ess the ESS to be used
	 * @return the prior for the current parameter ignoring some constants
	 * 
	 *@see StatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorPart( double ess )
	{
		double res = 0, hyper;
		for( int k, j, i = 0; i < constraints.length; i++ )
		{
			k = constraints[i].getNumberOfSpecificConstraints();
			hyper = ess / (double) k;
			for( j = 0; j < k; j++ )
			{
				res += hyper*constraints[i].getLambda(j);
			}
		}
		return res - ess*Math.log(Z);
	}

	/**
	 * This method returns some {@link String} using a given {@link NumberFormat}.
	 * @param nf the {@link NumberFormat}.
	 * @return some String representation
	 */
	public String toString( NumberFormat nf )
	{
		return Arrays.toString( indices );
	}

	public StringBuffer toXML()
	{
		StringBuffer erg = new StringBuffer( 500 * constraints.length );
		XMLParser.appendObjectWithTags( erg, indices, "used position" );
		XMLParser.appendObjectWithTags( erg, alphLen, "alphLen" );
		XMLParser.appendObjectWithTags( erg, Z, "normalization constant" );
		XMLParser.appendObjectWithTags( erg, constraints, "constraints" );
		if( cond != null )
		{
			StringBuffer b = new StringBuffer( 500 );
			XMLParser.appendObjectWithTags( b, cond, "cond" );
			XMLParser.addTags( b, "conditions" );
			erg.append( b );
		}
		XMLParser.addTags( erg, XML_TAG );
		return erg;
	}

	/**
	 * This method approximates the distribution either analytically or numerically.
	 * 
	 * @param s
	 *            the SequenceIterator used in normalization and numerical approximation
	 * @param algo
	 *            the choice of numerical approximation
	 * @param condition
	 *            the {@link TerminationCondition} for stopping the iterative algorithm
	 * @param sostream
	 *            a possibility for writing some information
	 * 
	 * @throws Exception
	 *             if something went wrong inside the algorithms
	 *
	 * @see MEMTools
	 */
	public void train(SequenceIterator s, byte algo, TerminationCondition condition, SafeOutputStream sostream ) throws Exception {
		Z = MEMTools.train(constraints, cond, s, algo, condition, sostream, alphLen);
		lnZ = Math.log(Z);
	}
}
