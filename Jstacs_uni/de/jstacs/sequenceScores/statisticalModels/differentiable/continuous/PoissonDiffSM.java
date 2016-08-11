package de.jstacs.sequenceScores.statisticalModels.differentiable.continuous;

import java.text.NumberFormat;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.random.RandomNumberGenerator;
import umontreal.iro.lecuyer.util.Num;

public class PoissonDiffSM extends AbstractDifferentiableStatisticalModel {

	private double llambda;
	private double lambda;
	
	private double ess;
	
	public PoissonDiffSM(double ess) throws IllegalArgumentException {
		super(new AlphabetContainer( new ContinuousAlphabet(0, Double.MAX_VALUE) ), 1);
		this.ess = ess;
		this.lambda = -1;
	}

	public PoissonDiffSM(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	
	public PoissonDiffSM clone() throws CloneNotSupportedException{
		PoissonDiffSM clone = (PoissonDiffSM) super.clone();
		return clone;
	}
	
	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		return 0;
	}

	@Override
	public double getLogNormalizationConstant() {
		return 0;
	}

	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex) throws Exception {
		return 0;
	}

	@Override
	public double getLogPriorTerm() {
		//TODO Gamma prior
		return 0;
	}

	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start) throws Exception {
		//TODO Gamma prior
	}

	@Override
	public double getESS() {
		return ess;
	}

	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		initializeFunctionRandomly(freeParams);
	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		RandomNumberGenerator rng = new RandomNumberGenerator();
		llambda = rng.nextGammaLog(1.0, 1.0);
		lambda = Math.exp(lambda);
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		double v = seq.continuousVal(start);//TODO discrete? ints?
		
		indices.add(0);
		partialDer.add(v - lambda);
		
		return v*llambda - Num.lnGamma(v+1) - lambda;
	}

	@Override
	public int getNumberOfParameters() {
		return 1;
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		return new double[]{llambda};
	}

	@Override
	public void setParameters(double[] params, int start) {
		llambda = params[start];
		lambda = Math.exp(llambda);
	}

	@Override
	public String getInstanceName() {
		return "Poisson";
	}

	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		double v = seq.continuousVal(start);//TODO discrete? ints?
		return v*llambda - Num.lnGamma(v+1) - lambda;
	}

	@Override
	public boolean isInitialized() {
		return lambda >= 0;
	}

	@Override
	public String toString(NumberFormat nf) {
		return "Poisson("+nf.format(lambda)+")";
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer sb = new StringBuffer();
		XMLParser.appendObjectWithTags(sb, llambda, "llambda");
		XMLParser.addTags(sb, "Poisson");
		return sb;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		this.alphabets = new AlphabetContainer( new ContinuousAlphabet(0, Double.MAX_VALUE) );
		this.length = 1;
		xml = XMLParser.extractForTag(xml, "Poisson");
		llambda = (Double) XMLParser.extractObjectForTags(xml, "llambda");
		lambda = Math.exp(llambda);
	}

}
