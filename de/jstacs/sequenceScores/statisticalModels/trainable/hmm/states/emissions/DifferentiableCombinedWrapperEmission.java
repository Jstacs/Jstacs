package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.CyclicMarkovModelDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.PeriodicHomogeneousModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * A wrapper class for emissions that combines multiple features. 
 * This emission can be used in numerical optimization as it used the methods from {@link DifferentiableStatisticalModel}.
 * However, some other methods that are related to sufficient statistics are currently not implemented and will throw an Exception.
 * Hence, Baum-Welch- as well as Viterbi-training (without numerical optimization) are currently not supported. 
 * 
 * @author Jens Keilwagen
 */
public class DifferentiableCombinedWrapperEmission implements DifferentiableEmission {

	private DifferentiableStatisticalModel[] model;
	private int[] offset, pOffset, length;
	private int l, minOffset;
	private double importance;
	private double[] logUniform;
	
	/**
	 * Create a {@link DifferentiableCombinedWrapperEmission}.
	 * 
	 * @param offset the offset for the start
	 * @param model the {@link DifferentiableStatisticalModel}s with a fixed length
	 * 
	 * @throws Exception if something goes wrong while initialization (e.g. not cloneable)
	 */
	public DifferentiableCombinedWrapperEmission( int[] offset, int[] length, DifferentiableStatisticalModel[] model ) throws Exception {
		this.model = new DifferentiableStatisticalModel[model.length];
		this.offset = new int[model.length];
		this.length = new int[model.length];
		this.pOffset = new int[model.length+1];
		AlphabetContainer con = model[0].getAlphabetContainer();
		if( !con.isSimple() ) {
			throw new IllegalArgumentException( "Possible only for simple AlphabetContainer" );
		} 

		minOffset = Integer.MAX_VALUE;
		int end = Integer.MIN_VALUE;
		for( int i = 0; i < model.length; i++ ) {
			if( !con.checkConsistency(model[i].getAlphabetContainer()) ) {
				throw new WrongAlphabetException( "Model " + i );
			}
			int len = model[i].getLength();
			if( length[i]==0 ) {
				throw new WrongLengthException( "Not possible to score variable length" );
			}
			if( len!=0 && len != length[i] ) {
				throw new WrongLengthException( "Model " + i + " with fixed length " + len + " cannot be used to score sequences with length " + length[i] );
			}
			this.model[i]=(DifferentiableStatisticalModel) model[i].clone();
			this.length[i] = length[i];
			this.offset[i]=offset[i];
			if( minOffset > offset[i] ) {
				minOffset = offset[i];
			}
			int e = offset[i]+length[i];
			if( end < e ) {
				end=e;
			}
		}
		l = end - minOffset;
		computeUniform();
	}

	private void computeUniform() {
		double a = Math.log( model[0].getAlphabetContainer().getAlphabetLengthAt(0));
		logUniform = new double[model.length];
		for( int i = 0; i<model.length; i++ ) {
			logUniform[i] = - length[i] * a;
		}
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link DifferentiableCombinedWrapperEmission} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DifferentiableCombinedWrapperEmission} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public DifferentiableCombinedWrapperEmission( StringBuffer xml) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, "DifferentiableCombinedWrapperEmission");
		model = (DifferentiableStatisticalModel[]) XMLParser.extractObjectForTags(xml, "model");
		offset = (int[]) XMLParser.extractObjectForTags(xml, "offset");
		length = (int[]) XMLParser.extractObjectForTags(xml, "length");
		pOffset = (int[]) XMLParser.extractObjectForTags(xml, "pOffset");
		importance = (Double) XMLParser.extractObjectForTags(xml, "importance");
		computeUniform();
		minOffset = Integer.MAX_VALUE;
		int end = Integer.MIN_VALUE;
		for( int i = 0; i < model.length; i++ ) {
			if( minOffset > offset[i] ) {
				minOffset = offset[i];
			}
			int e = offset[i]+model[i].getLength();
			if( end < e ) {
				end=e;
			}
		}
		l = end - minOffset;
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, model, "model");
		XMLParser.appendObjectWithTags(xml, offset, "offset");
		XMLParser.appendObjectWithTags(xml, length, "length");
		XMLParser.appendObjectWithTags(xml, pOffset, "pOffset");
		XMLParser.appendObjectWithTags(xml, importance, "importance");
		XMLParser.addTags(xml, "DifferentiableCombinedWrapperEmission");
		return xml;
	}

	@Override
	public DifferentiableCombinedWrapperEmission clone() throws CloneNotSupportedException {
		DifferentiableCombinedWrapperEmission clone = (DifferentiableCombinedWrapperEmission) super.clone();
		clone.offset = offset.clone();
		clone.length = length.clone();
		clone.pOffset = pOffset.clone();
		clone.model = new DifferentiableStatisticalModel[model.length];
		for( int i = 0; i < model.length; i++ ) {
			clone.model[i] = (DifferentiableStatisticalModel) model[i].clone();
		}
		clone.logUniform = logUniform.clone();
		return clone;
	}
	
	@Override
	public AlphabetContainer getAlphabetContainer() {
		return model[0].getAlphabetContainer();
	}

	@Override
	public void initializeFunctionRandomly() {
		try {
			for( int i = 0; i < model.length; i++ ) {
				model[i].initializeFunctionRandomly(false);
			}
			importance=0;
		} catch (Exception e) {
			throw new RuntimeException( e );
		}
	}
	
	public boolean setStartPhase( int phase ) {
		boolean set = false;
		for( int i = 0; i < model.length; i++ ) {
			if( model[i] instanceof PeriodicHomogeneousModel ) {
				set |= ((PeriodicHomogeneousModel)model[i]).setStartPhase(phase);
			}
		}
		return set;
	}

	@Override
	public double getLogProbFor(boolean forward, int startPos, int endPos, Sequence seq)
			throws OperationNotSupportedException {
		if( !forward ) throw new OperationNotSupportedException("Only forward allowed");
		try {
			double logScore=0;
			for( int i = 0; i < model.length; i++ ) {
				int s = startPos + offset[i];
				if( s < 0 || s+length[i] > seq.getLength() ) {
					logScore += logUniform[i];
				} else {
					logScore += model[i].getLogScoreFor(seq, s, s+length[i]-1);
				}
			}
			return importance + logScore;
		} catch (Exception ex) {
			ex.printStackTrace();
			throw new RuntimeException(ex);
		}
	}
	
	@Override
	public double getLogProbAndPartialDerivationFor(boolean forward, int startPos, int endPos, IntList indices,
			DoubleList partDer, Sequence seq) throws OperationNotSupportedException {
		if( !forward ) throw new OperationNotSupportedException("Only forward allowed");
		try {
			double logScore=0;
			for( int i = 0; i < model.length; i++ ) {
				int s = startPos + offset[i];
				if( s < 0 || s+length[i] > seq.getLength() ) {
					logScore += logUniform[i];
				} else {
					int indStart = indices.length();
					logScore += model[i].getLogScoreAndPartialDerivation(seq, s, s+length[i]-1, indices, partDer);
					indices.addToValues(indStart, indices.length(), pOffset[i]);
				}
			}
			partDer.add(1);
			indices.add(pOffset[model.length]);
			return importance + logScore;
		} catch (Exception ex) {
			throw new RuntimeException(ex);
		}
	}

	@Override
	public double getLogPriorTerm() {
		throw new RuntimeException("DifferentiableCombinedWrapperEmission");
	}

	@Override
	public String getNodeShape(boolean forward) {
		return "\"box\"";
	}

	@Override
	public String getNodeLabel(double weight, String name, NumberFormat nf) {
		return "\""+name+"\"";
	}

	@Override
	public String toString(NumberFormat nf) {
		StringBuffer sb = new StringBuffer("importance=" + nf.format(importance) );
		for( int i = 0; i < model.length; i++ ) {
			sb.append( "\noffset=" + offset[i] + "\nmodel:\n" + model[i].toString(nf) );
		}
		return sb.toString();
	}

	@Override
	public void fillCurrentParameter(double[] params) {
		try {
			for( int i = 0; i < model.length; i++ ) {
				double[] p = model[i].getCurrentParameterValues();
				System.arraycopy(p, 0, params, pOffset[i], p.length);
			}
			params[pOffset[model.length]]=importance;
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	public void setParameter(double[] params, int offset) {
		for( int i = 0; i < model.length; i++ ) {
			model[i].setParameters(params, offset+pOffset[i]);
		}
		importance=params[offset+pOffset[model.length]];
	}

	@Override
	public int setParameterOffset(int offset) {
		for( int i = 0; i < model.length; i++ ) {
			pOffset[i]=offset;
			offset += model[i].getNumberOfParameters();
		}
		pOffset[model.length] = offset;
		return offset+1;
	}

	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int offset) {
		throw new RuntimeException("DifferentiableCombinedWrapperEmission");
	}

	@Override
	public int getNumberOfParameters() {
		int num = 0;
		for( int i = 0; i < model.length; i++ ) {
			num += model[i].getNumberOfParameters();
		}
		return 1+num;
	}

	@Override
	public int getSizeOfEventSpace() {
		return (int) Math.pow(getAlphabetContainer().getAlphabetLengthAt(0),l);
	}
	
	@Override
	public void setParameters(Emission t) throws IllegalArgumentException {
		DifferentiableCombinedWrapperEmission d = (DifferentiableCombinedWrapperEmission) t;
		try {
			importance = d.importance;
			for( int i = 0; i < model.length; i++ ) {
				model[i].setParameters(d.model[i].getCurrentParameterValues(), 0);
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	@Override
	public void fillSamplingGroups(int parameterOffset, LinkedList<int[]> list) {
		throw new RuntimeException("DifferentiableCombinedWrapperEmission");
	}

	@Override
	public void resetStatistic() {
		throw new RuntimeException("DifferentiableCombinedWrapperEmission");
	}

	@Override
	public void addToStatistic(boolean forward, int startPos, int endPos, double weight, Sequence seq)
			throws OperationNotSupportedException {
		throw new RuntimeException("DifferentiableCombinedWrapperEmission");
	}

	@Override
	public void joinStatistics(Emission... emissions) {
		throw new RuntimeException("DifferentiableCombinedWrapperEmission");
	}

	@Override
	public void estimateFromStatistic() {
		throw new RuntimeException("DifferentiableCombinedWrapperEmission");
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
		//initialize parameters
		DataSet[] part = new DataSet[data.length];
		for( int i = 0; i < model.length; i++ ) {
			int len = model[i].getLength();
			if( len!=0 ) {
				int start = offset[i]-minOffset, dataLen = length[i];
				boolean use = true;
				if( model[i] instanceof PeriodicHomogeneousModel ) {
					if( len == 1 ) {
						int o = ((PeriodicHomogeneousModel) model[i]).getMaximalMarkovOrder();; 
						start -= o;
						dataLen += o;
					} else {
						use = false;
					}
				}
				if( use ) {
					for( int j = 0; j < data.length; j++ ) {
						part[j] = data[j].getInfixDataSet( start, dataLen );
					}
					model[i].initializeFunction( index, freeParams, part, weights );
				}
			} else {
				if( model[i] instanceof CyclicMarkovModelDiffSM ) {
					model[i].initializeFunctionRandomly(false);
				} else {
					double[] zeros = new double[model[i].getNumberOfParameters()];
					model[i].setParameters(zeros, 0);
				}
			}
		}
		importance = 0;

		//idea is to have neutral score close to zero (0)
		double[]  logScore = new double[data[index].getNumberOfElements()];
		for( int j = 0; j < logScore.length; j++ ) {
			logScore[j] = getLogProbFor( true,-minOffset, 0, data[index].getElementAt(j) );
		}
		Arrays.sort(logScore);
		importance = -logScore[(int) (0.1*logScore.length)];
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
		return minOffset;
	}
	
	public DifferentiableStatisticalModel getModel( int i ) throws CloneNotSupportedException {
		return (DifferentiableStatisticalModel) model[i].clone();
	}

	@Override
	public boolean isNormalized() {
		return false;
	}
	
	public void normalize() {
		for( int i = 0; i < model.length; i++ ) {
			if( model[i] instanceof BayesianNetworkDiffSM ) {
				((BayesianNetworkDiffSM)model[i]).normalize();
			}
		}
	}
}