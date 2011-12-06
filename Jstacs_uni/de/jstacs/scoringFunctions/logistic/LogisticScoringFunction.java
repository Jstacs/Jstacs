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
package de.jstacs.scoringFunctions.logistic;

import java.util.Arrays;
import java.util.Random;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;
import de.jstacs.scoringFunctions.AbstractScoringFunction;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class implements a logistic function. The score is computed by the following formula
 * {@latex.ilb \\[ s(\\underline{x}) = \\exp\\left(\\sum_i \\lambda_i \\cdot f_i(\\underline{x})\\right).\\]}
 * The constraints {@latex.inline $f_i(\\underline{x})$} are defined by {@link LogisticConstraint}s.  
 * 
 * @author Jens Keilwagen
 */
public class LogisticScoringFunction extends AbstractScoringFunction {

	protected LogisticConstraint[] constraint;
	protected double[] parameter;
	
	public LogisticScoringFunction( AlphabetContainer con, int length, LogisticConstraint... constraint ) throws CloneNotSupportedException {
		super( con, length );
		
		this.constraint = ArrayHandler.clone( constraint );
		this.parameter = new double[constraint.length];
	}

	public LogisticScoringFunction( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	public LogisticScoringFunction clone() throws CloneNotSupportedException {
		LogisticScoringFunction clone = (LogisticScoringFunction) super.clone();
		clone.parameter = parameter.clone();
		clone.constraint = ArrayHandler.clone( constraint );
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getCurrentParameterValues()
	 */
	@Override
	public double[] getCurrentParameterValues() throws Exception {
		return parameter.clone();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "logistic function";
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getLogScore(de.jstacs.data.Sequence, int)
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
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getLogScoreAndPartialDerivation(de.jstacs.data.Sequence, int, de.jstacs.utils.IntList, de.jstacs.utils.DoubleList)
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
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters() {
		return parameter.length;
	}

	@Override
	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		// TODO Auto-generated method stub
		initializeFunctionRandomly( freeParams );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunctionRandomly(boolean)
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
	 * @see de.jstacs.scoringFunctions.ScoringFunction#isInitialized()
	 */
	@Override
	public boolean isInitialized() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.ScoringFunction#setParameters(double[], int)
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
		alphabets = XMLParser.extractObjectForTags( xml, "alphabetContainer", AlphabetContainer.class );
		length = XMLParser.extractObjectForTags( xml, "length", Integer.class );
		constraint = XMLParser.extractObjectForTags( xml, "constraint", LogisticConstraint[].class );
		parameter = XMLParser.extractObjectForTags( xml, "parameter", double[].class );
	}
	
	private static final String XML_TAG = LogisticScoringFunction.class.getSimpleName(); 
	
	public String toString() {
		//TODO
		return Arrays.toString( parameter );
	}
}
