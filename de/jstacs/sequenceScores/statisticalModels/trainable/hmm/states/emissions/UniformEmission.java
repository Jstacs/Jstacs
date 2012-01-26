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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions;

import java.text.NumberFormat;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class implements a simple uniform emission.
 * 
 * @author Jens Keilwagen
 */
public class UniformEmission implements DifferentiableEmission {

	private AlphabetContainer con;
	private double logP;
	
	/**
	 * The main constructor for a {@link UniformEmission}.
	 * 
	 * @param con the {@link AlphabetContainer} for this emission
	 */
	public UniformEmission( AlphabetContainer con ) {
		this.con = con;
		logP = -Math.log( con.getAlphabetLengthAt(0) );
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link UniformEmission} out of an XML representation.
	 * 
	 * @param xml the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link UniformEmission} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public UniformEmission( StringBuffer xml ) throws NonParsableException {
		this( new AlphabetContainer( XMLParser.extractForTag( xml, XML_TAG ) ) );
	}
	
	public UniformEmission clone() throws CloneNotSupportedException {
		return (UniformEmission)super.clone();
	}

	private static final String XML_TAG = "UniformEmission";
	
	@Override
	public void addToStatistic(boolean forward, int startPos, int endPos, double weight, Sequence seq) throws OperationNotSupportedException {}

	@Override
	public void estimateFromStatistic() {}

	@Override
	public AlphabetContainer getAlphabetContainer() {
		return con;
	}

	@Override
	public double getLogPriorTerm() {
		return 0;
	}

	@Override
	public double getLogProbFor(boolean forward, int startPos, int endPos, Sequence seq) throws OperationNotSupportedException {
		return (endPos-startPos+1)*logP;
	}

	@Override
	public String getNodeLabel(double weight, String name, NumberFormat nf) {
		if(weight < 0){
			return "\""+name+"\"";
		}else{
			String res = "<<table border=\"0\" cellspacing=\"0\"><tr><td\">";
			if(weight < 0.5){
				res += "<font color=\"white\">"+name+"</font>";
			} else {
				res += name;
			}
			res += "</td></tr></table>>";
			return res;
		}		
	}

	@Override
	public String getNodeShape(boolean forward) {
		String res = "";
		if( getAlphabetContainer().isReverseComplementable() ) {
			res += "\"house\", orientation=";
			if( forward ) {
				res += "-";
			}
			res += "90";
		} else {
			res += "\"diamond\"";
		}
		return res;
	}

	@Override
	public void initializeFunctionRandomly() {}

	@Override
	public void joinStatistics(Emission... emissions) {}

	@Override
	public void resetStatistic() {}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = con.toXML();
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}

	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int offset) {}

	@Override
	public void fillCurrentParameter(double[] params) {}

	@Override
	public void fillSamplingGroups(int parameterOffset, LinkedList<int[]> list) {}

	@Override
	public double getLogProbAndPartialDerivationFor(boolean forward,
			int startPos, int endPos, IntList indices, DoubleList partDer,
			Sequence seq) throws OperationNotSupportedException {
		return getLogProbFor(forward, startPos, endPos, seq);
	}

	@Override
	public int getNumberOfParameters() {
		return 0;
	}

	@Override
	public int getSizeOfEventSpace() {
		return 0;
	}

	@Override
	public void setParameter(double[] params, int offset) {
	}

	@Override
	public int setParameterOffset(int offset) {
		return offset;
	}
	
	@Override
	public void setParameters(Emission t) throws IllegalArgumentException {
		if( !t.getClass().equals( getClass() ) || ((UniformEmission)t).logP != this.logP ) {
			throw new IllegalArgumentException( "The transitions are not comparable." );
		}		
	}
}