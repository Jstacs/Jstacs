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
package de.jstacs.sequenceScores.differentiable.logistic;

import java.util.Arrays;
import java.util.Random;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.differentiable.AbstractDifferentiableSequenceScore;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class implements a logistic function. The score is computed by the following formula
 * {@latex.ilb \\[ s(\\underline{x}) = \\exp\\left(\\sum_i \\lambda_i \\cdot f_i(\\underline{x})\\right).\\]}
 * The constraints {@latex.inline $f_i(\\underline{x})$} are defined by {@link LogisticConstraint}s.  
 * 
 * @author Jens Keilwagen
 */
public class LogisticDiffSS extends AbstractDifferentiableSequenceScore {

	private LogisticConstraint[] constraint;
	private double[] parameter;
	
	/**
	 * This is the main constructor to create {@link LogisticDiffSS} instance.
	 * 
	 * @param con
	 *            the {@link AlphabetContainer} of this instance
	 * @param length
	 *            the length of this instance, i.e. the length of
	 *            the modeled sequences
	 * @param constraint the constraints used in this instance
	 * 
	 * @throws CloneNotSupportedException if the constraints could not be cloned
	 */
	public LogisticDiffSS( AlphabetContainer con, int length, LogisticConstraint... constraint ) throws CloneNotSupportedException {
		super( con, length );
		
		this.constraint = ArrayHandler.clone( constraint );
		this.parameter = new double[constraint.length];
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link LogisticDiffSS} out of a {@link StringBuffer}
	 * .
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public LogisticDiffSS( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	public LogisticDiffSS clone() throws CloneNotSupportedException {
		LogisticDiffSS clone = (LogisticDiffSS) super.clone();
		clone.parameter = parameter.clone();
		clone.constraint = ArrayHandler.clone( constraint );
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getCurrentParameterValues()
	 */
	@Override
	public double[] getCurrentParameterValues() throws Exception {
		return parameter.clone();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "logistic function";
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScore(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor( Sequence seq, int start ) {
		double res = 0;
		for( int i = 0; i < constraint.length; i++ ) {
			res += parameter[i]*constraint[i].getValue( seq, start );
		}
		return res;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getLogScoreAndPartialDerivation(de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
	 */
	@Override
	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer ) {
		double res = 0, f;
		for( int i = 0; i < constraint.length; i++ ) {
			f = constraint[i].getValue( seq, start );
			res += parameter[i]*f;
			indices.add( i );
			partialDer.add( f );
		}
		return res;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters() {
		return parameter.length;
	}

	@Override
	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		initializeFunctionRandomly( freeParams );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#initializeFunctionRandomly(boolean)
	 */
	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		Random r = new Random();
		for( int i = 0; i < parameter.length; i++ ) {
			parameter[i] = r.nextGaussian();
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#isInitialized()
	 */
	@Override
	public boolean isInitialized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore#setParameters(double[], int)
	 */
	@Override
	public void setParameters( double[] params, int start ) {
		System.arraycopy( params, start, parameter, 0, parameter.length );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, alphabets, "alphabetContainer" );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		XMLParser.appendObjectWithTags( xml, constraint, "constraint" );
		XMLParser.appendObjectWithTags( xml, parameter, "parameter" );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}
	
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, XML_TAG );
		alphabets = (AlphabetContainer) XMLParser.extractObjectForTags( xml, "alphabetContainer" );
		length = XMLParser.extractObjectForTags( xml, "length", Integer.class );
		constraint = XMLParser.extractObjectForTags( xml, "constraint", LogisticConstraint[].class );
		parameter = XMLParser.extractObjectForTags( xml, "parameter", double[].class );
	}
	
	private static final String XML_TAG = LogisticDiffSS.class.getSimpleName(); 
	
	public String toString() {
		return Arrays.toString( parameter );
	}
}
