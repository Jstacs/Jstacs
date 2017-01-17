package projects.dream2016.mix;

import java.text.NumberFormat;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.differentiable.AbstractDifferentiableSequenceScore;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

public class SubDiffSS extends AbstractDifferentiableSequenceScore {

	DifferentiableSequenceScore dsm;
	int start;
	
	public SubDiffSS(AlphabetContainer con, int length, DifferentiableSequenceScore dsm, int start) throws IllegalArgumentException, CloneNotSupportedException {
		super(con, length);
		this.dsm = (DifferentiableSequenceScore) dsm.clone();
		this.start=start;
	}

	public SubDiffSS(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	public SubDiffSS clone() throws CloneNotSupportedException {
		SubDiffSS clone = (SubDiffSS) super.clone();
		clone.dsm = (DifferentiableSequenceScore) dsm.clone();
		return clone;
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
		dsm = (DifferentiableSequenceScore) XMLParser.extractObjectForTags(xml, "internal");
		
	}
}