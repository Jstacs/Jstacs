package projects.dream2016.mix;

import java.text.NumberFormat;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

public class SubDiffSM extends AbstractDifferentiableStatisticalModel {

	DifferentiableStatisticalModel dsm;
	int start;
	
	public SubDiffSM(AlphabetContainer con, int length, DifferentiableStatisticalModel dsm, int start) throws IllegalArgumentException, CloneNotSupportedException {
		super(con, length);
		this.dsm = (DifferentiableStatisticalModel) dsm.clone();
		this.start=start;
	}

	public SubDiffSM(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	public SubDiffSM clone() throws CloneNotSupportedException {
		SubDiffSM clone = (SubDiffSM) super.clone();
		clone.dsm = (DifferentiableStatisticalModel) dsm.clone();
		return clone;
	}

	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		return dsm.getSizeOfEventSpaceForRandomVariablesOfParameter(index);
	}

	@Override
	public double getLogNormalizationConstant() {
		return dsm.getLogNormalizationConstant();
	}

	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex) throws Exception {
		return dsm.getLogPartialNormalizationConstant(parameterIndex);
	}

	@Override
	public double getLogPriorTerm() {
		return dsm.getLogPriorTerm();
	}

	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start) throws Exception {
		dsm.addGradientOfLogPriorTerm(grad, start);
	}

	@Override
	public double getESS() {
		return dsm.getESS();
	}

	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		DataSet[] sub = new DataSet[data.length];
		for( int i = 0; i < sub.length; i++ ) {
			sub[i] = data[i].getInfixDataSet(start, dsm.getLength());
		}
		dsm.initializeFunction(index, freeParams, sub, weights);
	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		dsm.initializeFunctionRandomly(freeParams);
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		return dsm.getLogScoreAndPartialDerivation(seq, this.start+start, indices, partialDer);
	}

	@Override
	public int getNumberOfParameters() {
		return dsm.getNumberOfParameters();
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		return dsm.getCurrentParameterValues();
	}

	@Override
	public void setParameters(double[] params, int start) {
		dsm.setParameters(params, start);
	}

	@Override
	public String getInstanceName() {
		return dsm.getInstanceName() + " offset="+start;
	}

	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		return dsm.getLogScoreFor(seq, this.start+start);
	}

	@Override
	public boolean isInitialized() {
		return dsm.isInitialized();
	}

	@Override
	public String toString(NumberFormat nf) {
		return dsm.toString(nf) + " offset="+start;
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml =new StringBuffer();
		XMLParser.appendObjectWithTags(xml, alphabets, "alphabets");
		XMLParser.appendObjectWithTags(xml, length, "length");
		XMLParser.appendObjectWithTags(xml, start, "start");
		XMLParser.appendObjectWithTags(xml, dsm, "internal");
		XMLParser.addTags(xml, "SubDiffSM");
		return xml;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, "SubDiffSM");
		
		alphabets =(AlphabetContainer) XMLParser.extractObjectForTags(xml, "alphabets");
		length = (Integer) XMLParser.extractObjectForTags(xml, "length");
		start = (Integer) XMLParser.extractObjectForTags(xml, "start");
		dsm = (DifferentiableStatisticalModel) XMLParser.extractObjectForTags(xml, "internal");
		
	}
}