package de.jstacs.sequenceScores.statisticalModels.differentiable.continuous;

import java.text.NumberFormat;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.sequences.ArbitrarySequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import umontreal.ssj.util.Num;

public class NegativeBinomial extends AbstractDifferentiableStatisticalModel {

	
	private double r,lambda;
	private double p,beta;
	private boolean fixR;
	private double offset;
	
	
	public NegativeBinomial(Double r, double offset) throws IllegalArgumentException {
		super(new AlphabetContainer(new ContinuousAlphabet()), 1);
		if(r != null) {
			this.r = r;
			this.lambda = Math.log(r);
			this.fixR = true;
		}
		this.offset = offset;
	}

	
	public NegativeBinomial clone() throws CloneNotSupportedException{
		NegativeBinomial clone = (NegativeBinomial) super.clone();
		return clone;
	}
	
	
	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		return 1;
	}

	@Override
	public double getLogNormalizationConstant() {
		return 0;
	}

	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex) throws Exception {
		return Double.NEGATIVE_INFINITY;
	}

	@Override
	public double getLogPriorTerm() {
		return 0;
	}

	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start) throws Exception {
		// TODO Auto-generated method stub

	}

	@Override
	public double getESS() {
		return 0;
	}

	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		initializeFunctionRandomly(freeParams);
	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		p = 0.9;
		beta = Math.log(p/(1.0-p));
		
		if(!fixR) {
			r=1;
		}
		lambda = Math.log(r);

	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		double k = seq.continuousVal(start) - offset;
		
		int i=0;
		
		if(!fixR) {
			indices.add(i);
			partialDer.add( r*( Num.digamma(k+r) - Num.digamma(r) + Math.log(1-p) ) );
			i++;
		}
		
		indices.add(i);
		partialDer.add( k*(1.0-p) - r*p );
		
		return getLogScoreFor(seq, start);
		
		
	}

	@Override
	public int getNumberOfParameters() {
		if(fixR) {
			return 1;
		}else {
			return 2;
		}
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		if(fixR) {
			return new double[] {beta};
		}else {
			return new double[]{lambda,beta};
		}
	}

	@Override
	public void setParameters(double[] params, int start) {
		if(!fixR) {
			lambda = params[start++];
			r = Math.exp(lambda);
		}
		beta = params[start];
		p = Math.exp(beta)/(1.0+Math.exp(beta));
	}

	@Override
	public String getInstanceName() {
		return "NB";
	}

	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		double k = seq.continuousVal(start) - offset;
		return Num.lnGamma(k+r) - Num.lnGamma(k+1) - Num.lnGamma(r) + k*Math.log(p) + r*Math.log1p(-p);
	}

	@Override
	public boolean isInitialized() {
		return true;
	}

	@Override
	public String toString(NumberFormat nf) {
		return "NB("+r+","+(p)+")";
	}

	@Override
	public StringBuffer toXML() {
		return null;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		// TODO Auto-generated method stub

	}

	
	public static void main(String[] args) throws Exception {
		DataSet ds = new DataSet(new AlphabetContainer(new ContinuousAlphabet()), new SparseStringExtractor("/Users/dev/Desktop/current_projects/Matthias/lens.txt", '#'));
		
		NegativeBinomial nb = new NegativeBinomial(3.0,0);
		
		nb.initializeFunctionRandomly(false);
		
		System.out.println(nb);
		System.out.println( nb.getLogScoreFor(new ArbitrarySequence(new AlphabetContainer(new ContinuousAlphabet()), 2.0)) );
		
		//System.exit(1);
		
		ConstantDiffSM bg = new ConstantDiffSM(1);
		
		GenDisMixClassifierParameterSet params = new GenDisMixClassifierParameterSet(new AlphabetContainer(new ContinuousAlphabet()), 1, Optimizer.CONJUGATE_GRADIENTS_PRP, 1E-10, 1E-10, 1E-4, false, KindOfParameter.LAST, true, 1);
		GenDisMixClassifier cl = new GenDisMixClassifier(params, DoesNothingLogPrior.defaultInstance, LearningPrinciple.ML, nb,bg);
		
		cl.train(ds,ds);
		
		System.out.println(cl);
		
		nb = (NegativeBinomial) cl.getDifferentiableSequenceScore(0);
		
		System.out.println(nb);
		
	}
	
	
}
