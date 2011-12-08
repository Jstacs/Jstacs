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
package de.jstacs.scoringFunctions;

import java.util.Arrays;

import de.jstacs.NonParsableException;
import de.jstacs.data.MultiDimensionalSequence;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class implements a simple wrapper for multidimensional sequences.
 * 
 * The logarithm of the score given by {@link MultiDimensionalSequenceWrapperScoringFunction#getLogScoreFor(Sequence)} is defined as
 * {@latex.ilb \\[\\frac{1}{N}\\sum_{n=0}^N \\log p(\\underline{x_n} | \\underline{\\lambda}).\\]}
 * 
 * @author Jens Keilwagen
 * 
 * @see MultiDimensionalSequence
 */
public class MultiDimensionalSequenceWrapperScoringFunction extends AbstractScoringFunction {

	private ScoringFunction function;
	private IntList iList;
	private DoubleList dList;
	private double[] gradient;
	
	/**
	 * The main constructor.
	 * 
	 * @param function the internally used function
	 * 
	 * @throws IllegalArgumentException            
	 *            if the length is negative or does not match with {@link de.jstacs.data.AlphabetContainer#getPossibleLength()}
	 * @throws CloneNotSupportedException
	 *            if the function can not be cloned properly
	 */
	public MultiDimensionalSequenceWrapperScoringFunction( ScoringFunction function ) throws IllegalArgumentException, CloneNotSupportedException {
		super( function.getAlphabetContainer(), function.getLength() );
		this.function = function.clone();
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link MultiDimensionalSequenceWrapperScoringFunction} out of a {@link StringBuffer}
	 * .
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public MultiDimensionalSequenceWrapperScoringFunction( StringBuffer xml ) throws NonParsableException {
		super(xml);
	}

	private static final String XML_TAG = MultiDimensionalSequenceWrapperScoringFunction.class.getSimpleName();
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.AbstractScoringFunction#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, XML_TAG );
		function = XMLParser.extractObjectForTags( xml, "function", ScoringFunction.class );
		alphabets = function.getAlphabetContainer();
		length = function.getLength();
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, function, "function" );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}

	public MultiDimensionalSequenceWrapperScoringFunction clone() throws CloneNotSupportedException {
		MultiDimensionalSequenceWrapperScoringFunction clone = (MultiDimensionalSequenceWrapperScoringFunction) super.clone();
		clone.function = function.clone();
		if( gradient != null ) {
			clone.gradient = gradient.clone();
			clone.iList = iList.clone();
			clone.dList = dList.clone();
		}
		return clone;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getCurrentParameterValues()
	 */
	@Override
	public double[] getCurrentParameterValues() throws Exception {
		return function.getCurrentParameterValues();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "multidimensional wrapper of " + function.getInstanceName();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getLogScore(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor( Sequence seq, int start ) {
		double res;
		if( seq instanceof MultiDimensionalSequence ) {
			MultiDimensionalSequence mdSeq = (MultiDimensionalSequence) seq;
			int n = mdSeq.getNumberOfSequences();
			res = 0;
			for( int i = 0; i < n; i++ ) {
				res += function.getLogScoreFor( mdSeq.getSequence( i ), start );
			}
			res /= n;
		} else {
			res = function.getLogScoreFor( seq, start );
		}
		return res;
	}
	
	private void init() {
		gradient = new double[function.getNumberOfParameters()];
		iList = new IntList();
		dList = new DoubleList();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getLogScoreAndPartialDerivation(de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer ) {
		double res;
		if( seq instanceof MultiDimensionalSequence ) {
			MultiDimensionalSequence mdSeq = (MultiDimensionalSequence) seq;
			int n = mdSeq.getNumberOfSequences();
			res = 0;
			if( gradient == null ) {
				init();
			}
			Arrays.fill( gradient, 0 );
			for( int i = 0; i < n; i++ ) {
				iList.clear();
				dList.clear();
				res += function.getLogScoreAndPartialDerivation( mdSeq.getSequence( i ), start, iList, dList );
				
				for( int j = 0; j < iList.length(); j++ ) {
					gradient[ iList.get( j ) ] += dList.get( j );
				}
			}
			res /= n;
			for( int i = 0; i < gradient.length; i++ ) {
				if( gradient[i] != 0 ) {
					indices.add( i );
					partialDer.add( gradient[i]/n );
				}
			}
		} else {
			res = function.getLogScoreAndPartialDerivation( seq, start, indices, partialDer );
		}
		return res;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters() {
		return function.getNumberOfParameters();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunction(int, boolean, de.jstacs.data.Sample[], double[][])
	 */
	@Override
	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		function.initializeFunction( index, freeParams, data, weights );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunctionRandomly(boolean)
	 */
	@Override
	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {
		function.initializeFunctionRandomly( freeParams );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#isInitialized()
	 */
	@Override
	public boolean isInitialized() {
		return function.isInitialized();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#setParameters(double[], int)
	 */
	@Override
	public void setParameters(double[] params, int start) {
		function.setParameters( params, start );

	}
	
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "wrapper of " + function.getInstanceName() + ":\n" + function;
	}
}
