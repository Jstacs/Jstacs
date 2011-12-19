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

package de.jstacs.scoringFunctions.mix.motifSearch;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.io.XMLParser;
import de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class implements a position scoring function that enables the user to get a score without using a Sequence
 * object. It furthermore implements a kind of iterator that allows the user to iterate over all reasonable possibilities. 
 * 
 * @see Sequence
 * @see PositionScoringFunction#reset()
 * @see PositionScoringFunction#next()
 * @see PositionScoringFunction#getLogScoreForInternal()
 * @see PositionScoringFunction#getLogScoreAndPartialDerivationForInternal(IntList, DoubleList)
 * 
 * @author Jens Keilwagen
 */
public abstract class PositionScoringFunction extends AbstractNormalizableScoringFunction
{
	/**
	 * This array is used for some method of {@link DurationScoringFunction} that use an internal memory
	 * 
	 * @see DurationScoringFunction#reset()
	 * @see DurationScoringFunction#next()
	 * @see DurationScoringFunction#getInternalPosition(int[])
	 */
	protected int[] internal;	
	
	/**
	 * This constructor allows create instance with more than one dimension.
	 * This constructor can be used if the {@link AlphabetContainer} already exists.
	 * 
	 * @param con the {@link AlphabetContainer}
	 * @param length the number of dimensions, e.g. the length of the modeled sequences
	 */
	protected PositionScoringFunction( AlphabetContainer con, int length )
	{
		super( con, length );
		internal = new int[length];
		reset();
	}
	
	/**
	 * This is the main constructor that creates the {@link AlphabetContainer} internally.
	 * 
	 * @param min the minimal value
	 * @param max the maximal value
	 */
	protected PositionScoringFunction( int min, int max )
	{
		this( new AlphabetContainer( new DiscreteAlphabet( min, max ) ), 1 );
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link DurationScoringFunction} out of a {@link StringBuffer}.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	protected PositionScoringFunction( StringBuffer source ) throws NonParsableException
	{
		super( source );
	}
	
	public PositionScoringFunction clone() throws CloneNotSupportedException
	{
		PositionScoringFunction clone = (PositionScoringFunction) super.clone();
		clone.internal = internal.clone();
		return clone;
	}

	/**
	 * This method resets the iterator to the initial state (first reasonable output) so that it can be used again.
	 */
	public abstract void reset();
	
	/**
	 * This method steps to the next reasonable outcome if possible.
	 * 
	 * @return <code>true</code> if a next reasonable outcome could be set, otherwise <code>false</code>
	 */
	public abstract boolean next();
	
	/**
	 * Copies the current value of the internal iterator in the given array.
	 * 
	 * @param positions
	 *            the array containing the internal positions after calling this method
	 */
	public void getInternalPosition( int[] positions )
	{
		for( int i = 0; i < length; i++ )
		{
			positions[i] = internal[i]; 
		}
	}
	
	/**
	 * This method enables the user to get the log-score without using a sequence object by using the internal iterator.
	 * 
	 * @return the score
	 */
	public double getLogScoreForInternal()
	{
		return getLogScore( internal );
	}

	/**
	 * This method enables the user to get the log-score and the partial derivations without using a sequence object by using the internal iterator.
	 * 
	 * @param indices
	 *            a list for the indices of the parameters
	 * @param partialDer
	 *            a list of the partial derivations
	 * 
	 * @return the score
	 */
	public double getLogScoreAndPartialDerivationForInternal( IntList indices, DoubleList partialDer )
	{
		return getLogScoreAndPartialDerivation( indices, partialDer, internal );
	}
	
	/**
	 * This method enables the user to get the log-score without using a sequence object.
	 * 
	 * @param values
	 *            the values
	 * 
	 * @return the score
	 */
	public abstract double getLogScore( int... values );

	/**
	 * This method enables the user to get the log-score and the partial derivations without using a sequence object.
	 * 
	 * @param indices
	 *            a list for the indices of the parameters
	 * @param partialDer
	 *            a list of the partial derivations
	 * @param values
	 *            the values
	 * 
	 * @return the score
	 */
	public abstract double getLogScoreAndPartialDerivation( IntList indices, DoubleList partialDer, int... values );

	public double getLogScoreFor( Sequence seq, int start )
	{
		return getLogScore( getValuesFromSequence( seq, start ) );
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer )
	{
		return getLogScoreAndPartialDerivation( indices, partialDer, getValuesFromSequence( seq, start ) );
	}

	/**
	 * This method extracts the values form a sequence.
	 * 
	 * @param seq
	 *            the sequence
	 * @param start
	 *            the index of the startposition
	 * 
	 * @return the values of the sequence
	 */
	protected int[] getValuesFromSequence( Sequence seq, int start )
	{
		int[] res = new int[length];
		for( int i = 0; i < length; i++ )
		{
			res[i] = seq.discreteVal( start + i );
		}
		return res;
	}

	/**
	 * Generates a new set of positions in <code>positions</code> from the underlying density.
	 * 
	 * @param positions
	 *            the array containing the drawn positions after calling this method
	 */
	//public abstract void drawPosition( int[] positions ) throws Exception;
	
	
	/**
	 * This method returns <code>true</code> if the given <code>positions</code> are in the domain of the
	 * PositionScoringFunction.
	 * 
	 * @param positions
	 *            the positions to be tested
	 * 
	 * @return <code>true</code> if the given <code>positions</code> are in the domain of the
	 *         PositionScoringFunction
	 */
	public abstract boolean isPossible( int... positions );
	
	
	public StringBuffer toXML()
	{
		StringBuffer xml = new StringBuffer( 10000 );
		XMLParser.appendObjectWithTags( xml, alphabets, "AlphabetContainer" );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		XMLParser.appendObjectWithTags( xml, internal, "internal" );
		return xml;
	}
	
	protected void fromXML( StringBuffer xml ) throws NonParsableException
	{
		alphabets = XMLParser.extractObjectForTags( xml, "AlphabetContainer", AlphabetContainer.class );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		internal = XMLParser.extractObjectForTags( xml, "internal", int[].class );
	}
}
