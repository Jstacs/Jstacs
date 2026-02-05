package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions;

import java.text.NumberFormat;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.MultiDimensionalSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * A wrapper class for {@link Emission}s using only one dimension of a {@link MultiDimensionalSequence}.
 * 
 * @author Jens Keilwagen
 * 
 * @see MultiDimensionalSequence
 */
public class MultiDimensionalSequenceWrapperEmission implements DifferentiableEmission {

	private static String XML_TAG = "MDSWE";
	
	private DifferentiableEmission em;
	private int index;
	
	/**
	 * The main constructor.
	 * 
	 * @param em the emission
	 * @param dimension the dimension
	 * 
	 * @throws CloneNotSupportedException if the emission cannot be cloned
	 */
	public MultiDimensionalSequenceWrapperEmission( DifferentiableEmission em, int dimension) throws CloneNotSupportedException {
		this.em = em.clone();		
		if( dimension<0 ) {
			throw new IllegalArgumentException("The index has to be positive.");
		}
		index=dimension;
	}
	
	public MultiDimensionalSequenceWrapperEmission( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, XML_TAG);
		em = (DifferentiableEmission) XMLParser.extractObjectForTags(xml, "emission");
		index = (Integer) XMLParser.extractObjectForTags(xml, "dimension");
	}

	public MultiDimensionalSequenceWrapperEmission clone() throws CloneNotSupportedException {
		MultiDimensionalSequenceWrapperEmission clone = (MultiDimensionalSequenceWrapperEmission) super.clone();
		clone.em = em.clone();
		return clone;
	}
	
	@Override
	public AlphabetContainer getAlphabetContainer() {
		return em.getAlphabetContainer();
	}

	@Override
	public void initializeFunctionRandomly() {
		em.initializeFunctionRandomly();
	}

	@Override
	public double getLogProbFor(boolean forward, int startPos, int endPos, Sequence seq)
			throws OperationNotSupportedException {
		return em.getLogProbFor(forward, startPos, endPos, ((MultiDimensionalSequence) seq).getSequence(index));
	}

	@Override
	public double getLogPriorTerm() {
		return em.getLogPriorTerm();
	}

	@Override
	public void resetStatistic() {
		em.resetStatistic();
	}

	@Override
	public void addToStatistic(boolean forward, int startPos, int endPos, double weight, Sequence seq)
			throws OperationNotSupportedException {
		em.addToStatistic(forward, startPos, endPos, weight, ((MultiDimensionalSequence)seq).getSequence(index));
	}

	@Override
	public void joinStatistics(Emission... emissions) {
		Emission[] e = new Emission[emissions.length];
		for( int i = 0; i < emissions.length; i++ ) {
			e[i] = ((MultiDimensionalSequenceWrapperEmission) emissions[i]).em;
		}
		em.joinStatistics(e);
	}

	@Override
	public void estimateFromStatistic() {
		em.estimateFromStatistic();
	}

	@Override
	public String getNodeShape(boolean forward) {
		return em.getNodeShape(forward);
	}

	@Override
	public String getNodeLabel(double weight, String name, NumberFormat nf) {
		return em.getNodeLabel(weight, name, nf);
	}

	@Override
	public void setParameters(Emission t) throws IllegalArgumentException {
		em.setParameters(((MultiDimensionalSequenceWrapperEmission)t).em);
	}

	@Override
	public String toString(NumberFormat nf) {
		return "MultiDimensionalSequenceWrapper dim="+ index+"\n"+em.toString();
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, em, "emission");
		XMLParser.appendObjectWithTags(xml, index, "dimension");
		XMLParser.addTags(xml, XML_TAG);
		return xml;
	}

	@Override
	public void fillCurrentParameter(double[] params) {
		em.fillCurrentParameter(params);
	}

	@Override
	public void setParameter(double[] params, int offset) {
		em.setParameter(params, offset);
	}

	@Override
	public int setParameterOffset(int offset) {
		return em.setParameterOffset(offset);
	}

	@Override
	public double getLogProbAndPartialDerivationFor(boolean forward, int startPos, int endPos, IntList indices,
			DoubleList partDer, Sequence seq) throws OperationNotSupportedException {
		return em.getLogProbAndPartialDerivationFor(forward, startPos, endPos, indices, partDer, ((MultiDimensionalSequence)seq).getSequence(index));
	}

	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int offset) {
		em.addGradientOfLogPriorTerm(grad, offset);
	}

	@Override
	public void fillSamplingGroups(int parameterOffset, LinkedList<int[]> list) {
		em.fillSamplingGroups(parameterOffset, list);
	}

	@Override
	public int getNumberOfParameters() {
		return em.getNumberOfParameters();
	}

	@Override
	public int getSizeOfEventSpace() {
		return em.getSizeOfEventSpace();
	}

	@Override
	public boolean isNormalized() {
		return em.isNormalized();
	}
}