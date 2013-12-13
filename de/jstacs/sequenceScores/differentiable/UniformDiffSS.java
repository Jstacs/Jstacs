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

package de.jstacs.sequenceScores.differentiable;

import java.text.NumberFormat;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This {@link DifferentiableSequenceScore} does nothing. So it is possible to save
 * parameters in an optimization.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class UniformDiffSS extends AbstractDifferentiableSequenceScore {
	
	/**
	 * This is the main constructor that creates an instance of a
	 * {@link UniformDiffSS} that models each sequence uniformly.
	 * 
	 * @param alphabets
	 *            the {@link AlphabetContainer}
	 * @param length
	 *            the length of the modeled sequences
	 */
	public UniformDiffSS(AlphabetContainer alphabets, int length) {
		super( alphabets, length );
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link UniformDiffSS} out of its XML
	 * representation as returned by {@link #fromXML(StringBuffer)}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public UniformDiffSS(StringBuffer xml) throws NonParsableException {
		super( xml );
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getInstanceName()
	 */
	public String getInstanceName() {
		return "uniform";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScore(de.jstacs.data
	 * .Sequence, int)
	 */
	public double getLogScoreFor(Sequence seq, int start) {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScoreAndPartialDerivation
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
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getNumberOfParameters()
	 */
	public int getNumberOfParameters() {
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#setParameters(double[],
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
	 * @see UniformDiffSS#extractFurtherInformation(StringBuffer)
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
		b.append( getFurtherInformation() );
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
	 * @see UniformDiffSS#getFurtherInformation()
	 */
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {}

	protected void fromXML(StringBuffer xml) throws NonParsableException {
		StringBuffer b = XMLParser.extractForTag(xml, getClass().getSimpleName());
		length = XMLParser.extractObjectForTags(b, "length", int.class );
		alphabets = (AlphabetContainer) XMLParser.extractObjectForTags(b, "alphabets");
		extractFurtherInformation( b );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.SequenceScore#toString(java.text.NumberFormat)
	 */
	@Override
	public String toString( NumberFormat nf ) {
		return "value 0 (zero) for each sequence ";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getCurrentParameterValues()
	 */
	public double[] getCurrentParameterValues() throws Exception {
		return new double[0];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#isInitialized()
	 */
	public boolean isInitialized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#initializeFunctionRandomly
	 * (boolean)
	 */
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#initializeFunction(int,
	 * boolean, de.jstacs.data.DataSet[], double[][])
	 */
	public void initializeFunction( int index, boolean meila, DataSet[] data, double[][] weights ) {
		// does nothing
	}
}
