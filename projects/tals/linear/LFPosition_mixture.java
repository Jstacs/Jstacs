package projects.tals.linear;

import java.text.NumberFormat;

import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.differentiable.AbstractDifferentiableSequenceScore;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

public class LFPosition_mixture extends AbstractDifferentiableSequenceScore {

	private double beta1, beta2, beta3, a1, a2, b1, b2;
	
	private double[][] precomputed;
	
	public LFPosition_mixture() throws IllegalArgumentException {
		super(DNAAlphabetContainer.SINGLETON, 0);
		precomputed = new double[50][];
	}

	public LFPosition_mixture(StringBuffer xml) throws NonParsableException {
		super(xml);
		precomputed = new double[50][];
	}

	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		initializeFunctionRandomly(freeParams);

	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		beta1 = Math.random()*2.0-1.0;
		beta2 = Math.random()*2.0-1.0;
		beta3 = Math.random()*2.0-1.0;
		a1 = -Math.random()*3.0;
		a2 = Math.random()*3.0;
		b1 = -(Math.random()*0.5+0.25);
		b2 = -(Math.random()*0.5+0.25);
		precomputed = new double[50][];
	}


	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		double x = start/(double)seq.getLength();
		double y = (seq.getLength()-start)/(double)seq.getLength();
		
		double p1 = Math.exp(beta1);
		double p2 = Math.exp(beta2);
		double p3 = Math.exp(beta3);
		double sum = p1+p2+p3;
		p1 /= sum; p2 /= sum; p3 /= sum;
		
		double f1 = 1.0/(1.0 + Math.exp( -a1*(x+b1) ));
		double f2 = 1.0/(1.0 + Math.exp( -a2*(y+b2) ));
		
		double score =  p1*f1 + p2*f2 + p3*1.0;
		
		indices.add(0);
		partialDer.add(p1*(1.0-p1)*f1 - p1*p2*f2 - p1*p3);
		
		indices.add(1);
		partialDer.add(-p1*p2*f1 + p2*(1.0-p2)*f2 - p2*p3);
		
		indices.add(2);
		partialDer.add(-p1*p3*f1 - p2*p3*f2 + p3*(1.0-p3));
		
		indices.add(3);
		partialDer.add(p1*f1*(1.0-f1)*(x+b1));
		
		indices.add(4);
		partialDer.add(p2*f2*(1.0-f2)*(y+b2));
		
		indices.add(5);
		partialDer.add(p1*f1*(1.0-f1)*a1);
		
		indices.add(6);
		partialDer.add(p2*f2*(1.0-f2)*a2);
		
		return score;
	}

	@Override
	public int getNumberOfParameters() {
		return 7;
	}

	@Override
	public double[] getCurrentParameterValues() {
		return new double[]{beta1,beta2,beta3,a1,a2,b1,b2};
	}

	@Override
	public void setParameters(double[] params, int start) {
		beta1 = params[start];
		beta2 = params[start+1];
		beta3 = params[start+2];
		a1 = params[start+3];
		a2 = params[start+4];
		b1 = params[start+5];
		b2 = params[start+6];
		precomputed = new double[50][];
	}

	@Override
	public String getInstanceName() {
		return getClass().getSimpleName();
	}

	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		if(precomputed[seq.getLength()] == null){
			precomputed[seq.getLength()] = new double[seq.getLength()];
			for(int start2 = 0; start2<seq.getLength();start2++){
				double x = start2/(double)seq.getLength();
				double y = (seq.getLength()-start2)/(double)seq.getLength();
				
				double p1 = Math.exp(beta1);
				double p2 = Math.exp(beta2);
				double p3 = Math.exp(beta3);
				double sum = p1+p2+p3;
				p1 /= sum; p2 /= sum; p3 /= sum;
				
				double f1 = 1.0/(1.0 + Math.exp( -a1*(x+b1) ));
				double f2 = 1.0/(1.0 + Math.exp( -a2*(y+b2) ));
				
				precomputed[seq.getLength()][start2] = p1*f1 + p2*f2 + p3*1.0;
			}
		}
		return precomputed[seq.getLength()][start];
	}
	

	public double getLogScoreFor(int len, int start) {
		if(precomputed[len] == null){
			precomputed[len] = new double[len];
			for(int start2 = 0; start2<len;start2++){
				double x = start2/(double)len;
				double y = (len-start2)/(double)len;
				
				double p1 = Math.exp(beta1);
				double p2 = Math.exp(beta2);
				double p3 = Math.exp(beta3);
				double sum = p1+p2+p3;
				p1 /= sum; p2 /= sum; p3 /= sum;
				
				double f1 = 1.0/(1.0 + Math.exp( -a1*(x+b1) ));
				double f2 = 1.0/(1.0 + Math.exp( -a2*(y+b2) ));
				
				precomputed[len][start2] = p1*f1 + p2*f2 + p3*1.0;
			}
		}
		return precomputed[len][start];
	}

	@Override
	public boolean isInitialized() {
		return true;
	}

	@Override
	public String toString(NumberFormat nf) {
		double p1 = Math.exp(beta1);
		double p2 = Math.exp(beta2);
		double p3 = Math.exp(beta3);
		double sum = p1+p2+p3;
		p1 /= sum; p2 /= sum; p3 /= sum;
		
		return p1+"*1/(1+exp("+(-a1)+"*(x + "+b1+"))) + "+p2+"*1/(1+exp("+(-a2)+"*(y + "+b2+"))) + "+p3;
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, getCurrentParameterValues(), "params");
		XMLParser.addTags(xml, "mixture");
		return xml;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, "mixture");
		setParameters((double[])XMLParser.extractObjectForTags(xml, "params"), 0);
		this.alphabets = DNAAlphabetContainer.SINGLETON;
		this.length = 0;
	}

}
