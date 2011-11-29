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
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.io.XMLParser;

/**
 * This class is the super class for all one dimensional position scoring functions that can be used as durations for semi Markov models.
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.scoringFunctions.mix.motifSearch.HiddenMotifsMixture
 */
public abstract class DurationScoringFunction extends PositionScoringFunction
{
	/**
	 * The equivalent sample size.
	 */
	protected double ess;
	/**
	 * The minimal value.
	 */
	protected int min;
	/**
	 * The maximal value.
	 */
	protected int max;
	/**
	 * The difference of maximal and minimal value.
	 */
	protected int delta;
	
	/**
	 * The default constructor.
	 * 
	 * @param min the minimal value
	 * @param max the maximal value
	 * @param ess the equivalent sample size
	 */
	protected DurationScoringFunction(int min, int max, double ess )
	{
		super(min, max);
		setMinMax( min, max );
		if( ess < 0 )
		{
			throw new IllegalArgumentException( "The given ess has to be non-negative." );
		}
		this.ess = ess;
		reset();
	}
	
	private void setMinMax( int min, int max ) {
		this.min = min;
		if( min < 0 )
		{
			throw new IllegalArgumentException( "The given minimum is below 0." );
		}	
		this.max = max;
		delta = max - min;
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
	protected DurationScoringFunction(StringBuffer source) throws NonParsableException
	{
		super(source);
		reset();
	}
	
	private static final String XML_TAG = "DurationScoringFunction";
	
	public StringBuffer toXML()
	{
		StringBuffer b = super.toXML();
		XMLParser.appendObjectWithTags( b, max, "max" );
		XMLParser.appendObjectWithTags( b, min, "min" );
		XMLParser.appendObjectWithTags( b, ess, "ess" );
		XMLParser.addTags( b, XML_TAG );
		return b;		
	}
	
	protected void fromXML( StringBuffer xml ) throws NonParsableException
	{
		StringBuffer b = XMLParser.extractForTag( xml, XML_TAG );
		super.fromXML(b);
		setMinMax( XMLParser.extractObjectForTags( b, "min", int.class ), XMLParser.extractObjectForTags( b, "max", int.class ) );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		ess = XMLParser.extractObjectForTags( b, "ess", double.class );
	}

	public void reset()
	{
		internal[0] = min;		
	}
	
	public boolean next()
	{
		internal[0]++;
		return internal[0] <= max;
	}

	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index)
	{
		// TODO Auto-generated method stub
		// return 0;
		throw new RuntimeException( "???" );
	}

	public final double getEss()
	{
		return ess;
	}

	public boolean isPossible( int... positions )
	{
		return positions != null && positions.length == 1 && positions[0] >= min && positions[0] <= max;
	}
	
	/**
	 * Returns the minimal value that can be scored.
	 * 
	 * @return the minimal value that can be scored.
	 */
	public final int getMin()
	{
		return min;
	}
	
	/**
	 * Returns the maximal value that can be scored.
	 * 
	 * @return the maximal value that can be scored.
	 */
	public final int getMax()
	{
		return max;
	}
	
	/**
	 * Returns the number of different possibilities that can be scored.
	 * 
	 * @return the number of different possibilities that can be scored.
	 */
	public int getNumberOfPossibilities()
	{
		return max-min+1;
	}
	
	/**
	 * This method set special parameters that lead to an uniform distribution.
	 */
	public abstract void initializeUniformly();

	/**
	 * This method adjust the parameter based on the given statistic.
	 * 
	 * @param length an array containing length values
	 * @param weight an array containing corresponding weight values
	 */
	public abstract void adjust( int[] length, double[] weight );
	
	/**
	 * This method modifies the underlying {@link AlphabetContainer}. This might be necessary if the motif length changed. 
	 * 
	 * @param delta the change
	 * 
	 * @see de.jstacs.motifDiscovery.Mutable#modify(int, int)
	 * @see de.jstacs.motifDiscovery.MutableMotifDiscoverer#modifyMotif(int, int, int)
	 */
	public void modify( int delta ) {
		if( delta != 0 ) {
			setMinMax( min, max + delta );
			alphabets = new AlphabetContainer( new DiscreteAlphabet( min, max ) );
		}
	}
	
	public final double getLogNormalizationConstant()
	{
		return 0;
	}

	public final double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception
	{
		return Double.NEGATIVE_INFINITY;
	}
	
	/**
	 * This method returns the distribution in <a href="http://www.r-project.org/">R</a> notation.
	 * 
	 * @param distributionName the name of the distribution, e.g., &quot;p&quot;
	 * 
	 * @return the distribution in R notation
	 * 
	 * @see de.jstacs.utils.REnvironment
	 */
	protected abstract String getRNotation( String distributionName );

	public String toString() {
		return getRNotation( "p" );
	}
	
}
