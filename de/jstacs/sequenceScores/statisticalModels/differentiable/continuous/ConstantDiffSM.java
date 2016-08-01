package de.jstacs.sequenceScores.statisticalModels.differentiable.continuous;

import java.text.NumberFormat;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

public class ConstantDiffSM extends AbstractDifferentiableStatisticalModel {

	public ConstantDiffSM(int length) throws IllegalArgumentException {
		super(new AlphabetContainer(new ContinuousAlphabet()), length);
	}

	public ConstantDiffSM(StringBuffer xml) throws NonParsableException {
		super(xml);
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
		return 0;
	}

	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start) throws Exception {

	}

	@Override
	public double getESS() {
		return 0;
	}

	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {

	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {

	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		return 0;
	}

	@Override
	public int getNumberOfParameters() {
		return 0;
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		return new double[0];
	}

	@Override
	public void setParameters(double[] params, int start) {

	}

	@Override
	public String getInstanceName() {
		return "const";
	}

	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		return 0;
	}

	@Override
	public boolean isInitialized() {
		return true;
	}

	@Override
	public String toString(NumberFormat nf) {
		return "const";
	}

	@Override
	public StringBuffer toXML() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		// TODO Auto-generated method stub

	}

}
