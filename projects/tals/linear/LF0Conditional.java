package projects.tals.linear;

import java.text.NumberFormat;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.differentiable.AbstractDifferentiableSequenceScore;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import projects.tals.RVDSequence;

public class LF0Conditional extends AbstractDifferentiableSequenceScore {

	private double[][] beta;
	private int[] indexMap;
	private AlphabetContainer rvds;
	private double[] prior;
	
	public LF0Conditional(AlphabetContainer rvds,String[] specific, double[] priorProbs) throws IllegalArgumentException {
		super(DNAAlphabetContainer.SINGLETON, 1);
		beta = new double[specific.length+1][(int) getAlphabetContainer().getAlphabetLengthAt(1)];
		int k=1;
		indexMap = new int[(int) rvds.getAlphabetLengthAt(0)];
		for(int i=0;i<rvds.getAlphabetLengthAt(0);i++){
			String curr = rvds.getSymbol(0, i);
			for(int j=0;j<specific.length;j++){
				if(curr.equals(specific[j])){
					indexMap[i] = k;
					k++;
					break;
				}
			}
		}
		this.rvds = rvds;
		this.prior = new double[priorProbs.length];
		for(int i=0;i<prior.length;i++){
			prior[i] = Math.log(priorProbs[i]);
		}
	}

	public LF0Conditional clone() throws CloneNotSupportedException{
		LF0Conditional clone = (LF0Conditional) super.clone();
		clone.beta = ArrayHandler.clone(beta);
		clone.indexMap = indexMap.clone();
		return clone;
	}
	
	public LF0Conditional(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		
		initializeFunctionRandomly(freeParams);
		
	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		for(int j=0;j<beta[0].length;j++){
			beta[0][j] = Math.random()*2.0 - 1.0;
		}
		for(int i=1;i<beta.length;i++){
			System.arraycopy(beta[0], 0, beta[i], 0, beta[0].length);
		}
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		int idx = seq.discreteVal(start);
		
		ReferenceSequenceAnnotation data_anno =(ReferenceSequenceAnnotation)seq.getSequenceAnnotationByType("reference", 0);
		RVDSequence rvd_seq=(RVDSequence)data_anno.getReferenceSequence();
		
		int index_tj=indexMap[rvd_seq.discreteVal(start)];
		
		
		double bet = beta[0][idx] + beta[index_tj][idx];
		
		indices.add(index_tj*beta[0].length + idx);
		partialDer.add(1);
		
		indices.add(idx);
		partialDer.add(1);
		
		
		return prior[seq.discreteVal(start)] + bet;
	}

	@Override
	public int getNumberOfParameters() {
		return beta.length*beta[0].length;
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		double[] params = new double[getNumberOfParameters()];
		for(int i=0,k=0;i<beta.length;i++){
			for(int j=0;j<beta[i].length;j++,k++){
				params[k] = beta[i][j];
			}
		}
		return params;
	}

	@Override
	public void setParameters(double[] params, int start) {
		for(int i=0,k=0;i<beta.length;i++){
			for(int j=0;j<beta[i].length;j++,k++){
				beta[i][j] = params[k+start];
			}
		}
	}

	@Override
	public String getInstanceName() {
		return "LF0Conditional";
	}

	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		
		ReferenceSequenceAnnotation data_anno =(ReferenceSequenceAnnotation)seq.getSequenceAnnotationByType("reference", 0);
		RVDSequence rvd_seq=(RVDSequence)data_anno.getReferenceSequence();
		
		int index_tj=indexMap[rvd_seq.discreteVal(start)];
		
		return prior[seq.discreteVal(start)] + beta[0][seq.discreteVal(start)] + beta[index_tj][seq.discreteVal(start)];
		
	}

	@Override
	public boolean isInitialized() {
		return true;
	}

	@Override
	public String toString(NumberFormat nf) {
		StringBuffer sb = new StringBuffer();
		sb.append("OT: "+Arrays.toString(beta[0])+"\n");
		for(int i=0;i<indexMap.length;i++){
			if(indexMap[i] > 0){
				sb.append(rvds.getSymbol(0, i)+": "+Arrays.toString(beta[indexMap[i]])+"\n");
			}
		}
		return sb.toString();
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer sb = new StringBuffer();
		XMLParser.appendObjectWithTags(sb, beta, "beta");
		XMLParser.appendObjectWithTags(sb, indexMap, "indexMap");
		XMLParser.appendObjectWithTags(sb, rvds, "rvds");
		XMLParser.appendObjectWithTags(sb, prior, "prior");
		XMLParser.addTags(sb, "LF0Conditional");
		return sb;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		this.alphabets = DNAAlphabetContainer.SINGLETON;
		this.length = 1;
		xml = XMLParser.extractForTag(xml, "LF0Conditional");
		this.beta = (double[][]) XMLParser.extractObjectForTags(xml, "beta");
		indexMap = (int[]) XMLParser.extractObjectForTags(xml, "indexMap");
		rvds = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "rvds");
		prior = (double[]) XMLParser.extractObjectForTags(xml, "prior");
	}

	public double[] getSpecs(Sequence rvds) {
		
		int index_tj=indexMap[rvds.discreteVal(0)];
		
		
		double[] specs = new double[(int) alphabets.getAlphabetLengthAt(0)];
		
		for(int i=0;i<specs.length;i++){
			specs[i] = prior[i] + beta[0][i] + beta[index_tj][i];
		}
		
		return specs;
	}

}
