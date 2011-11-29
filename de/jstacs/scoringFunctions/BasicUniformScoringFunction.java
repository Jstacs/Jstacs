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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.scoringFunctions;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This {@link ScoringFunction} does nothing. So it is possible to save
 * parameters in an optimization.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class BasicUniformScoringFunction extends AbstractScoringFunction {
	
	/**
	 * This is the main constructor that creates an instance of a
	 * {@link BasicUniformScoringFunction} that models each sequence uniformly.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param length
	 *            the length of the modeled sequences
	 */
	public BasicUniformScoringFunction(AlphabetContainer alphabets, int length) {
		super( alphabets, length );
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link BasicUniformScoringFunction} out of its XML
	 * representation as returned by {@link #fromXML(StringBuffer)}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public BasicUniformScoringFunction(StringBuffer xml) throws NonParsableException {
		super( xml );
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getInstanceName()
	 */
	public String getInstanceName() {
		return "uniform";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#getLogScore(de.jstacs.data
	 * .Sequence, int)
	 */
	public double getLogScoreFor(Sequence seq, int start) {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#getLogScoreAndPartialDerivation
	 * (de.jstacs.data.Sequence, int, de.jstacs.utils.IntList,
	 * de.jstacs.utils.DoubleList)
	 */
	public double getLogScoreAndPartialDerivation(Sequence seq, int start,
			IntList indices, DoubleList dList) {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getNumberOfParameters()
	 */
	public int getNumberOfParameters() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#setParameters(double[],
	 * int)
	 */
	public void setParameters( double[] params, int start ) {
	}

	/**
	 * This method is used to append further information of the instance to the
	 * XML representation. This method is designed to allow subclasses to add
	 * information to the XML representation.
	 * 
	 * @return the further information as XML code in a {@link StringBuffer}
	 * 
	 * @see BasicUniformScoringFunction#extractFurtherInformation(StringBuffer)
	 */
	protected StringBuffer getFurtherInformation() {
		return new StringBuffer( 1 );
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer b = new StringBuffer( 1000 );
		XMLParser.appendObjectWithTags( b, length, "length" );
		XMLParser.appendObjectWithTags( b, alphabets, "alphabets" );
		XMLParser.addTags( b, getClass().getSimpleName() );
		return b;
	}
	
	/**
	 * This method is the opposite of {@link #getFurtherInformation()}. It
	 * extracts further information of the instance from a XML representation.
	 * 
	 * @param xml
	 *            the {@link StringBuffer} containing the information to be
	 *            extracted as XML code
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} could not be parsed
	 * 
	 * @see BasicUniformScoringFunction#getFurtherInformation()
	 */
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {}

	protected void fromXML(StringBuffer xml) throws NonParsableException {
		StringBuffer b = XMLParser.extractForTag(xml, getClass().getSimpleName());
		length = XMLParser.extractObjectForTags(b, "length", int.class );
		alphabets = XMLParser.extractObjectForTags(b, "alphabets", AlphabetContainer.class );
		extractFurtherInformation( b );
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "value 0 (zero) for each sequence ";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#getCurrentParameterValues()
	 */
	public double[] getCurrentParameterValues() throws Exception {
		return new double[0];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#isInitialized()
	 */
	public boolean isInitialized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#initializeFunctionRandomly
	 * (boolean)
	 */
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunction(int,
	 * boolean, de.jstacs.data.Sample[], double[][])
	 */
	public void initializeFunction( int index, boolean meila, Sample[] data, double[][] weights ) {
		// does nothing
	}
}
