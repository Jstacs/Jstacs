package projects.tals.linear;

import java.text.NumberFormat;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.QuickScanningSequenceScore;
import de.jstacs.sequenceScores.differentiable.AbstractDifferentiableSequenceScore;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import projects.tals.RVDSequence;

public class LFModularConditional9C extends AbstractDifferentiableSequenceScore {

	protected LFPosition_mixture position;
	protected LFSpecificity_parallel_cond9C specificity;
	protected LF0Conditional lf0;
	protected double[] a, b;

	
	public LFModularConditional9C(LF0Conditional lf0, LFSpecificity_parallel_cond9C specificity, LFPosition_mixture lfPosition, int numberOfGroups) throws IllegalArgumentException {
		super(DNAAlphabetContainer.SINGLETON, 0);
		this.lf0 = lf0;
		this.specificity = specificity;
		this.position = lfPosition;
		this.a = new double[numberOfGroups];
		this.b = new double[numberOfGroups];

	}
	
	public DifferentiableSequenceScore getSpecModel(){
		return specificity;
	}
	
	
	public LFModularConditional9C(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	

	



	public LFModularConditional9C clone() throws CloneNotSupportedException{
		LFModularConditional9C clone = (LFModularConditional9C) super.clone();
		clone.lf0 = lf0.clone();
		clone.specificity = specificity.clone();
		clone.a = a.clone();
		clone.b = b.clone();
		if(position != null){
			clone.position = (LFPosition_mixture)position.clone();
		}

		return clone;
	}

	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		initializeFunctionRandomly(freeParams);
	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		lf0.initializeFunctionRandomly(freeParams);
		specificity.initializeFunctionRandomly(freeParams);
		if(position != null){ position.initializeFunctionRandomly(freeParams); }
		for(int i=0;i<a.length;i++){
			a[i] = Math.log(Math.random()+0.5);
			b[i] = Math.random()-0.5;
		}
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		
		String gs = seq.getSequenceAnnotationByType("intgroup", 0).getIdentifier();
		
		
		String mask = null;
		SequenceAnnotation mann = seq.getSequenceAnnotationByType("mask", 0);
		if(mann != null){
			mask = mann.getIdentifier();
		}
		
		int group = Integer.parseInt( gs );
		
		double ag = Math.exp( a[group] );
		
		int off = indices.length();
		
		double score = 0.0;
		if(mask ==null || mask.charAt(start) == 'O'){
			score += lf0.getLogScoreAndPartialDerivation(seq, start, indices, partialDer);
		}
		
		partialDer.multiply(off, partialDer.length(), ag);
		indices.addToValues(off, indices.length(), a.length+b.length);
		
		off = indices.length();
		
		for(int i=start+1;i<seq.getLength();i++){
			if(mask == null || mask.charAt(i) == 'O'){
				
				int cs = partialDer.length();
				double spec = specificity.getLogScoreAndPartialDerivation(seq, i, indices, partialDer);
				int ci = partialDer.length();
				double pos = 1.0;
				int cp = partialDer.length();
				if(position != null){
					pos = position.getLogScoreAndPartialDerivation(seq, i, indices, partialDer);	
				}
				partialDer.multiply(cs, ci, pos*ag);
				indices.addToValues(cs, ci, a.length+b.length+lf0.getNumberOfParameters());
				
				partialDer.multiply(ci, cp, spec*pos*ag);
				indices.addToValues(ci, cp, a.length+b.length+lf0.getNumberOfParameters()+specificity.getNumberOfParameters());
				
				partialDer.multiply(cp, partialDer.length(), spec*ag);
				indices.addToValues(cp, indices.length(), a.length+b.length+lf0.getNumberOfParameters()+specificity.getNumberOfParameters());
				
				score += spec*pos;
				
			}
		}
		
		//partialDer.multiply(off, partialDer.length(), ag);
		//indices.addToValues(off, indices.length(), a.length+b.length+lf0.getNumberOfParameters());
		
		
		indices.add(group);
		partialDer.add(score*ag);
		
		indices.add(a.length + group);
		partialDer.add(1);
		
		score *= ag;
		score += b[group];
		
		return score;
		
	}

	@Override
	public int getNumberOfParameters() {
		return a.length + b.length + lf0.getNumberOfParameters() + specificity.getNumberOfParameters() 
			+ (position == null ? 0 : position.getNumberOfParameters());
	}
	

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		double[] params = new double[getNumberOfParameters()];
		System.arraycopy(a, 0, params, 0, a.length);
		System.arraycopy(b, 0, params, a.length, b.length);
		int off = a.length + b.length;
		System.arraycopy(lf0.getCurrentParameterValues(), 0, params, off, lf0.getNumberOfParameters());
		off += lf0.getNumberOfParameters();
		System.arraycopy(specificity.getCurrentParameterValues(), 0, params, off, specificity.getNumberOfParameters());
		off += specificity.getNumberOfParameters();
		
		if(position != null){
			
			System.arraycopy(position.getCurrentParameterValues(), 0, params, off, position.getNumberOfParameters());
		}
		return params;
	}

	@Override
	public void setParameters(double[] params, int start) {
		System.arraycopy(params, start, a, 0, a.length);
		start += a.length;
		System.arraycopy(params, start, b, 0, b.length);
		start += b.length;
		lf0.setParameters(params, start);
		start += lf0.getNumberOfParameters();
		specificity.setParameters(params, start);
		start += specificity.getNumberOfParameters();
		
		if(position != null){
		
			position.setParameters(params, start);
		}
	}

	@Override
	public String getInstanceName() {
		return "LFModular";
	}

	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		
		String mask = null;
		SequenceAnnotation mann = seq.getSequenceAnnotationByType("mask", 0);
		if(mann != null){
			mask = mann.getIdentifier();
		}
		SequenceAnnotation ann = seq.getSequenceAnnotationByType("intgroup", 0);
		int group = 0;
		if(ann != null){
			String gs = ann.getIdentifier();
			group = Integer.parseInt( gs );
		}
		
		double score = 0.0;
		if(mask == null || mask.charAt(start) == 'O'){
			score += lf0.getLogScoreFor(seq, start);
		}
		
		for(int i=start+1;i<seq.getLength();i++){
			if(mask == null || mask.charAt(i) == 'O'){
				
				double spec = specificity.getLogScoreFor(seq, i);
			
				double pos = (position == null ? 1.0 : position.getLogScoreFor(seq, i));
				
				score += spec*pos;
			}
		}
		
		score *= Math.exp( a[group] );
		score += b[group];
		
		return score;
		
	}
	
	public String getLogScorePoswiseFor(Sequence seq, int start) {
		String outScore="";
		String mask = null;
		SequenceAnnotation mann = seq.getSequenceAnnotationByType("mask", 0);
		if(mann != null){
			mask = mann.getIdentifier();
		}
		SequenceAnnotation ann = seq.getSequenceAnnotationByType("intgroup", 0);
		int group = 0;
		if(ann != null){
			String gs = ann.getIdentifier();
			group = Integer.parseInt( gs );
		}
		
		double score = 0.0;
		if(mask == null || mask.charAt(start) == 'O'){
			score += lf0.getLogScoreFor(seq, start);
			outScore="lf0: "+score;
		}
		
		for(int i=start+1;i<seq.getLength();i++){
			if(mask == null || mask.charAt(i) == 'O'){
				
				double spec = specificity.getLogScoreFor(seq, i);
			
				double pos = (position == null ? 1.0 : position.getLogScoreFor(seq, i));
				
				score += spec*pos;
				outScore+=", spec_pos"+i+": "+spec;
			}
		}
		
		score *= Math.exp( a[group] );
		score += b[group];
		
		return outScore;
	}

	@Override
	public boolean isInitialized() {
		return true;
	}

	@Override
	public String toString(NumberFormat nf) {
		return Arrays.toString(a)+"\n"+Arrays.toString(b)+"\n\n"+lf0.toString()+"\n\n"+specificity.toString()+"\n\n"+position;
	}

	public double[][] toPWM(RVDSequence rvds){
		double[][] pwm = new double[rvds.getLength()+1][];
		pwm[0] = lf0.getSpecs(rvds);
		
		for(int i=1;i<pwm.length;i++){
			pwm[i] = specificity.getSpecs(rvds,i);
			
			double pos = (position == null ? 1.0 : position.getLogScoreFor(rvds.getLength()+1, i));
			
			for(int j=0;j<pwm[i].length;j++){
				pwm[i][j] *= pos;
			}
			
		}
		
		for(int i=0;i<pwm.length;i++){
			for(int j=0;j<pwm[i].length;j++){
				pwm[i][j] /= pwm.length;
			}
		}
		return pwm;
		
	}
	
	
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, specificity, "specificity");
		
		XMLParser.appendObjectWithTags(xml, position, "position");
		XMLParser.appendObjectWithTags(xml, a, "a");
		XMLParser.appendObjectWithTags(xml, b, "b");
		XMLParser.appendObjectWithTags(xml, lf0, "lf0");
		XMLParser.addTags(xml, "LFMod");
		return xml;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, "LFMod");
		specificity = (LFSpecificity_parallel_cond9C) XMLParser.extractObjectForTags(xml, "specificity");
		
		position = (LFPosition_mixture) XMLParser.extractObjectForTags(xml, "position");
		this.a = (double[]) XMLParser.extractObjectForTags(xml, "a");
		this.b = (double[]) XMLParser.extractObjectForTags(xml, "b");
		this.lf0 = (LF0Conditional) XMLParser.extractObjectForTags(xml, "lf0");
		this.alphabets = DNAAlphabetContainer.SINGLETON;
		this.length = 0;
	}


}
