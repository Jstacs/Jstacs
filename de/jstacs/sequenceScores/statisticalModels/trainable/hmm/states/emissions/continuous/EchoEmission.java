package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.continuous;

import java.text.NumberFormat;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * A wrapper class for {@link Emission}s echoing the values of the sequence.
 * 
 * @author Jens Keilwagen
 * 
 * @see Sequence#continuousVal(int)
 */
public class EchoEmission implements DifferentiableEmission {

	private static String XML_TAG = "ECHO";
	
	private AlphabetContainer con;
	
	/**
	 * The main constructor.
	 * 
	 * @param con the alphabet container
	 * 
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} does not allow continuous {@link Sequence}s of arbitrary length 
	 */
	public EchoEmission( AlphabetContainer con ) throws WrongAlphabetException {
		if( con.getPossibleLength()!= 0 || con.isDiscreteAt(0) ) {
			throw new WrongAlphabetException("The Alphabet needs to be continuous.");
		}
		this.con = con;
	}
	
	public EchoEmission( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, XML_TAG);
		con = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "con");
	}

	public EchoEmission clone() throws CloneNotSupportedException {
		return (EchoEmission) super.clone();
	}
	
	@Override
	public AlphabetContainer getAlphabetContainer() {
		return con;
	}

	@Override
	public void initializeFunctionRandomly() {
	}

	@Override
	public double getLogProbFor(boolean forward, int startPos, int endPos, Sequence seq)
			throws OperationNotSupportedException {
		double res = 0;
		while( startPos <= endPos ) {
				res += seq.continuousVal(startPos++);
		}
		return res;
	}

	@Override
	public double getLogPriorTerm() {
		return 0;
	}

	@Override
	public void resetStatistic() {
	}

	@Override
	public void addToStatistic(boolean forward, int startPos, int endPos, double weight, Sequence seq)
			throws OperationNotSupportedException {
	}

	@Override
	public void joinStatistics(Emission... emissions) {
	}

	@Override
	public void estimateFromStatistic() {
	}

	@Override
	public String getNodeShape(boolean forward) {
		return "\"box\"";
	}

	@Override
	public String getNodeLabel( double weight, String name, NumberFormat nf ) {
		return "\""+name+"\"";
	}

	@Override
	public void setParameters(Emission t) throws IllegalArgumentException {
	}

	@Override
	public String toString(NumberFormat nf) {
		return "EchoEmission";
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, con, "con");
		XMLParser.addTags(xml, XML_TAG);
		return xml;
	}

	@Override
	public void fillCurrentParameter(double[] params) {
	}

	@Override
	public void setParameter(double[] params, int offset) {
	}

	@Override
	public int setParameterOffset(int offset) {
		return offset;
	}

	@Override
	public double getLogProbAndPartialDerivationFor(boolean forward, int startPos, int endPos, IntList indices,
			DoubleList partDer, Sequence seq) throws OperationNotSupportedException {
		double res = 0;
		while( startPos < endPos ) {
				res += seq.continuousVal(startPos++);
		}
		return res;
	}

	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int offset) {
	}

	@Override
	public void fillSamplingGroups(int parameterOffset, LinkedList<int[]> list) {
	}

	@Override
	public int getNumberOfParameters() {
		return 0;
	}

	@Override
	public int getSizeOfEventSpace() {
		return 0;
	}
}