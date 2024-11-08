package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * A wrapper class for emissions based on {@link DifferentiableStatisticalModel} with a fixed length. 
 * This emission can be used in numerical optimization as it used the methods from {@link DifferentiableStatisticalModel}.
 * However, some other methods that are related to sufficient statistics are currently not implemented and will throw an Exception.
 * Hence, Baum-Welch- as well as Viterbi-training (without numerical optimization) are currently not supported. 
 * 
 * @author Jens Keilwagen
 */
public class DifferentiableSMWrapperEmission implements DifferentiableEmission {

	private DifferentiableStatisticalModel model;
	private int pOffset, offset, l;
	private double logUniform;
	
	/**
	 * Create a {@link DifferentiableSMWrapperEmission}.
	 * 
	 * @param offset the offset for the start
	 * @param model a {@link DifferentiableStatisticalModel} with a fixed length
	 * 
	 * @throws CloneNotSupportedException if the model can not be cloned
	 */
	public DifferentiableSMWrapperEmission( int offset, DifferentiableStatisticalModel model ) throws CloneNotSupportedException {
		if( !model.getAlphabetContainer().isSimple() ) {
			throw new IllegalArgumentException("Possible only for simple AlphabetContainer");
		} 
		if( model.getLength()==0 ) {
			throw new IllegalArgumentException("Possible only for models with fixed length");
		}
		this.model=(DifferentiableStatisticalModel) model.clone();
		this.offset=offset;
		computeUniform();
	}

	private void computeUniform() {
		l = model.getLength();
		logUniform = - l * Math.log( model.getAlphabetContainer().getAlphabetLengthAt(0));
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link DifferentiableSMWrapperEmission} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DifferentiableSMWrapperEmission} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public DifferentiableSMWrapperEmission( StringBuffer xml) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, "DifferentiableSMWrapperEmission");
		model = (DifferentiableStatisticalModel) XMLParser.extractObjectForTags(xml, "model");
		offset = (Integer) XMLParser.extractObjectForTags(xml, "offset");
		pOffset = (Integer) XMLParser.extractObjectForTags(xml, "pOffset");
		computeUniform();
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, model, "model");
		XMLParser.appendObjectWithTags(xml, offset, "offset");
		XMLParser.appendObjectWithTags(xml, pOffset, "pOffset");
		XMLParser.addTags(xml, "DifferentiableSMWrapperEmission");
		return xml;
	}

	@Override
	public DifferentiableSMWrapperEmission clone() throws CloneNotSupportedException {
		DifferentiableSMWrapperEmission clone = (DifferentiableSMWrapperEmission) super.clone();
		clone.model = (DifferentiableStatisticalModel) model.clone();
		return clone;
	}
	
	@Override
	public AlphabetContainer getAlphabetContainer() {
		return model.getAlphabetContainer();
	}

	@Override
	public void initializeFunctionRandomly() {
		try {
			model.initializeFunctionRandomly(false);
		} catch (Exception e) {
			throw new RuntimeException( e );
		}
	}

	@Override
	public double getLogProbFor(boolean forward, int startPos, int endPos, Sequence seq)
			throws OperationNotSupportedException {
		if( !forward ) throw new OperationNotSupportedException("Only forward allowed");
		int s = startPos + offset;//, e=endPos+offset+l-1;
		//System.out.println(startPos+".."+endPos + "\t" + s + ".." + e); //System.exit(1);
		try {
			double logScore;
			if( s < 0 || s+l > seq.getLength() ) logScore = logUniform;
			else logScore = model.getLogScoreFor(seq, s);
			return logScore;
		} catch (Exception ex) {
			ex.printStackTrace();
			throw new RuntimeException(ex);
		}
	}
	
	@Override
	public double getLogProbAndPartialDerivationFor(boolean forward, int startPos, int endPos, IntList indices,
			DoubleList partDer, Sequence seq) throws OperationNotSupportedException {
		if( !forward ) throw new OperationNotSupportedException("Only forward allowed");
		int s = startPos + offset;//, e=endPos+offset+l-1;
		//System.out.println(startPos+".."+endPos + " -> " + s + ".." + e); //System.exit(1);
		try {
			int indStart = indices.length();
			double logScore;
			if( s < 0 || s+l > seq.getLength() ) logScore = logUniform;
			else logScore = model.getLogScoreAndPartialDerivation(seq, s, indices, partDer);
			indices.addToValues(indStart, indices.length(), pOffset);
			return logScore;
		} catch (Exception ex) {
			throw new RuntimeException(ex);
		}
	}

	@Override
	public double getLogPriorTerm() {
		return model.getLogPriorTerm();
	}

	@Override
	public String getNodeShape(boolean forward) {
		return "\"box\"";
	}

	@Override
	public String getNodeLabel(double weight, String name, NumberFormat nf) {
		// TODO Auto-generated method stub
		return "\""+name+"\"";
	}

	@Override
	public String toString(NumberFormat nf) {
		return "offset=" + offset + "\nmodel:\n" + model.toString(nf);
	}

	@Override
	public void fillCurrentParameter(double[] params) {
		try {
			double[] p = model.getCurrentParameterValues();
			System.arraycopy(p, 0, params, pOffset, p.length);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	public void setParameter(double[] params, int offset) {
		model.setParameters(params, offset+pOffset);
	}

	@Override
	public int setParameterOffset(int offset) {
		pOffset=offset;
		return pOffset+getNumberOfParameters();
	}

	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int offset) {
		try {
			model.addGradientOfLogPriorTerm(grad, offset);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	public int getNumberOfParameters() {
		return model.getNumberOfParameters();
	}

	@Override
	public int getSizeOfEventSpace() {
		return (int) Math.pow(getAlphabetContainer().getAlphabetLengthAt(0),l);
	}
	
	@Override
	public void setParameters(Emission t) throws IllegalArgumentException {
		DifferentiableSMWrapperEmission d = (DifferentiableSMWrapperEmission) t;
		try {
			model.setParameters(d.model.getCurrentParameterValues(), 0);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	@Override
	public void fillSamplingGroups(int parameterOffset, LinkedList<int[]> list) {
		throw new RuntimeException("DifferentiableSMWrapperEmission");
	}

	@Override
	public void resetStatistic() {
		throw new RuntimeException("DifferentiableSMWrapperEmission");
	}

	@Override
	public void addToStatistic(boolean forward, int startPos, int endPos, double weight, Sequence seq)
			throws OperationNotSupportedException {
		throw new RuntimeException("DifferentiableSMWrapperEmission");
	}

	@Override
	public void joinStatistics(Emission... emissions) {
		throw new RuntimeException("DifferentiableSMWrapperEmission");
	}

	@Override
	public void estimateFromStatistic() {
		throw new RuntimeException("DifferentiableSMWrapperEmission");
	}
	
	public void initializeUniformly() throws Exception {
		model.initializeFunctionRandomly(false);
		double[] params = model.getCurrentParameterValues();
		Arrays.fill(params, 0);
		model.setParameters(params, 0);
	}
	
	/**
	 * Allows to initialize the internal model using some data.
	 * 
	 * @param index
	 *            the index of the class the {@link DifferentiableSequenceScore} models
	 * @param freeParams
	 *            indicates whether the (reduced) parameterization is used
	 * @param data
	 *            the data sets
	 * @param weights
	 *            the weights of the sequences in the data sets
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see DifferentiableStatisticalModel#initializeFunction(int, boolean, DataSet[], double[][])
	 */
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		model.initializeFunction( index, freeParams, data, weights );
	}

	/**
	 * Returns the length of sequences this internal model can score.
	 * 
	 * @return the length of sequences the internal model can score
	 *  
	 * @see DifferentiableStatisticalModel#getLength()
	 */
	public int getLength() {
		return l;
	}

	/**
	 * Returns the offset that is used for the internal model to extract the sequence to be scored.
	 * 
	 * @return the offset that is used for the internal model
	 */
	public int getOffset() {
		return offset;
	}

	@Override
	public boolean isNormalized() {
		return model.isNormalized();
	}
}