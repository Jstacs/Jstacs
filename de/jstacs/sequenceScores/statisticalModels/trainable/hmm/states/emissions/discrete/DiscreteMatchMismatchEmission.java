package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete;

import java.text.NumberFormat;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.alphabets.IUPACDNAAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sampling.SamplingFromStatistic;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.Emission;

/**
 * @author Jens Keilwagen
 */
public class DiscreteMatchMismatchEmission extends ReferenceSequenceDiscreteEmission {

	private static IUPACDNAAlphabet iupac = IUPACDNAAlphabet.SINGLETON;///XXX?
	
	private static AlphabetContainer dummyRefCon;
	static {
		try {
			dummyRefCon = new AlphabetContainer( new DiscreteAlphabet(true,  ".") );
		} catch (IllegalArgumentException | DoubleSymbolException e) {
			e.printStackTrace();
		}
	}
	
	private DiscreteMatchMismatchEmission ref;
	private DifferentiableEmission[] e;
	
	//create the first emission
	public DiscreteMatchMismatchEmission( AlphabetContainer con, int refIdx, double ess, DifferentiableEmission[] e ) throws IllegalArgumentException, CloneNotSupportedException, WrongAlphabetException {
		this( con, refIdx, getHyperParams(ess, 1, 2), null, null, e );
	}
	
	//create the second, third, ... emission
	public DiscreteMatchMismatchEmission( AlphabetContainer con, int refIdx, double ess, DiscreteMatchMismatchEmission ref ) throws IllegalArgumentException, CloneNotSupportedException, WrongAlphabetException {
		this( con, refIdx, getHyperParams(ess, 1, 2), null, ref, null );
	}
		
	public DiscreteMatchMismatchEmission clone() throws CloneNotSupportedException {
		DiscreteMatchMismatchEmission clone = (DiscreteMatchMismatchEmission) super.clone();
		if( refIdx == 0 ) {
			clone.e = ArrayHandler.clone(e);
			ref = clone;
			
			//System.out.println("create a new clone " + clone.toString() );
		} else {
			clone.setReference( ref.ref );
			//System.out.println(refIdx + "\trefering to the clone " + clone.ref.toString() );
		}
		return clone;
	}
	
	private DiscreteMatchMismatchEmission( AlphabetContainer con, int refIdx, double[][] hyperParams, double[][] initHyperParams, DiscreteMatchMismatchEmission ref, DifferentiableEmission[] e ) throws IllegalArgumentException, CloneNotSupportedException, WrongAlphabetException {
		super( con, dummyRefCon, refIdx, hyperParams, initHyperParams== null ? hyperParams : initHyperParams );
		if( hyperParams.length!= 1 || hyperParams[0].length!= 2 ) {
			throw new IllegalArgumentException("hyperParams has to be double[1][2]");
		}
		if( ref == null ) {
			this.e = ArrayHandler.clone(e);
		} else {
			setReference(ref);
		}
	}
	
	public void setReference( DiscreteMatchMismatchEmission ref ) throws CloneNotSupportedException {
		this.ref = ref;
		this.e = ArrayHandler.clone(ref.e);
	}
	
	public DiscreteMatchMismatchEmission(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	protected void appendFurtherInformation( StringBuffer xml ) {
		super.appendFurtherInformation(xml);
		if( refIdx == 0 ) {
			XMLParser.appendObjectWithTags(xml, e, "emission");
		}
	}
	
	//TODO protected void extractFurtherInformation(StringBuffer xml) throws NonParsableException
	
	protected int getConditionIndex( boolean forward, int seqPos, Sequence seq ) {
		return 0;
	}
	
	protected int getIndex( int seqPos, Sequence seq ) {
		return 
				//iupac.isPart(seq.discreteVal(seqPos), getReferenceSequence( seq ).discreteVal( refIdx ))
				seq.discreteVal(seqPos) == getReferenceSequence( seq ).discreteVal( refIdx )
				? 0 : 1;
	}
	
	@Override
	public String toString( NumberFormat nf ) {
		StringBuffer sb = new StringBuffer();
		sb.append( "P(R_" + refIdx + "=Match) = " + nf.format(probs[0][0]) + "\tP(R_" + refIdx + "=Mismatch) = " + nf.format(probs[0][1]) + "\n\n" );
		for(int i = 0; i < e.length; i++ ) {
			if( e[i] != null ) {
				sb.append( (i==0?"Match":"Mismatch") + ":\n" + e[i].toString(nf) + "\n" );
			}
		}
		return sb.toString();
	}
	
	public void initializeFunctionRandomly() {
		super.initializeFunctionRandomly();
		if(refIdx == 0){
			for(int i=0;i<e.length;i++){
				if( e[i] != null ) {
					e[i].initializeFunctionRandomly();
				}
			}
		} else {
			for(int i=0;i<e.length;i++){
				if( e[i] != null ) {
					e[i].setParameters(ref.e[i]);
				}
			}
		}
	}
	
	public void addToStatistic( boolean forward, int startPos, int endPos, double weight, Sequence seq ) throws OperationNotSupportedException {
		int s, e;
		Sequence current;
		
		if( forward ) {
			current = seq;
			s = startPos;
			e = endPos;
		} else {
			current = seq.reverseComplement();
			int len = current.getLength();
			s = len - endPos -1;
			e = len - startPos -1;
		}
		
		//System.out.println(current);
		while( s <= e ) {
			int condIdx = getConditionIndex( forward, s, seq);//TODO error? seq. vs. current
			int index = getIndex(s, current);
			statistic[condIdx][index] += weight;
			if( this.e[index] != null ) {
				//System.out.println(s + "\t" + current.discreteVal(s) + "\t" + index + "\t" + weight);
				this.e[index].addToStatistic(true, s, s, weight, current);
			}
			s++;
		}
	}
	
	public void estimateFromStatistic() {
		super.estimateFromStatistic();
		if(refIdx == 0){
			try {
				for(int i=0;i<e.length;i++){
					if( e[i] != null ) {
						//TODO A=T and C=G
						AbstractConditionalDiscreteEmission a = (AbstractConditionalDiscreteEmission) e[i];
						a.statistic[0][0] = a.statistic[0][3] = (a.statistic[0][0]+a.statistic[0][3])/2d;
						a.statistic[0][1] = a.statistic[0][2] = (a.statistic[0][1]+a.statistic[0][2])/2d;
						
						e[i].estimateFromStatistic();
					}
				}
			} catch (Exception e1) {
				throw new RuntimeException(e1);
			}
		} else {
			for(int i=0;i<e.length;i++){
				if( e[i] != null ) {
					e[i].setParameters(ref.e[i]);
				}
			}
		}
	}
	
	public void drawParametersFromStatistic() {
		super.drawParametersFromStatistic();
		if(refIdx == 0){
			try {
				for(int i=0;i<e.length;i++){
					if( e[i] != null ) {
						((SamplingFromStatistic) e[i]).drawParametersFromStatistic();
					}
				}
			} catch (Exception e1) {
				throw new RuntimeException(e1);
			}
		} else {
			for(int i=0;i<e.length;i++){
				if( e[i] != null ) {
					e[i].setParameters(ref.e[i]);
				}
			}
		}
	}
	
	public double getLogPosteriorFromStatistic() {
		double res = super.getLogPosteriorFromStatistic();
		if( refIdx == 0 ) {
			for(int i=0;i<e.length;i++){
				if( e[i] != null ) {
					((SamplingFromStatistic) e[i]).getLogPosteriorFromStatistic();
				}
			}
		}
		return res;
	}

	public void resetStatistic() {
		super.resetStatistic();
		for(int i=0;i<e.length;i++){
			if( e[i] != null ) {
				e[i].resetStatistic();
			}
		}
	}
	
	public void joinStatistics(Emission... emissions) {
		super.joinStatistics(emissions);
		Emission[] sub = new Emission[emissions.length+(refIdx == 0?0:1)];
		for(int i = 0; i < e.length; i++ ) {
			if( e[i] != null ) {
				for( int j = 0; j < emissions.length; j++ ) {
					sub[j] = ((DiscreteMatchMismatchEmission) emissions[j]).e[i];
				}
				if( refIdx != 0 ) {
					sub[emissions.length] = ref.e[i];
				}
				//System.out.println(refIdx + "\t" + i + "\t" + Arrays.toString( ((AbstractConditionalDiscreteEmission) e[i]).statistic[0]) );
				e[i].joinStatistics(sub);
				
				//System.out.println(refIdx + "\t" + i + "\t" + Arrays.toString( ((AbstractConditionalDiscreteEmission) (refIdx==0?this:ref).e[i]).statistic[0]) );
			}
		}
		//join stats for different refIdx
	}
	
	public void setParameters( Emission t ) throws IllegalArgumentException {
		super.setParameters(t);
		for(int i = 0; i < e.length; i++ ) {
			if( e[i] != null ) {
				e[i].setParameters( ((DiscreteMatchMismatchEmission) t).e[i] );
			}
		}
	}
	
	public double getLogPriorTerm() {
		double res = super.getLogPriorTerm();
		if(refIdx == 0){
			for(int i=0;i<e.length;i++){
				if( e[i] != null ) {
					res += e[i].getLogPriorTerm();
				}
			}
		}
		return res;
	}
	
	
	
	
	public int getNumberOfParameters() {
		int n = super.getNumberOfParameters();
		if(refIdx == 0){
			for( int i = 0; i < e.length; i++ ) {
				if( e[i] != null ) {
					n += e[i].getNumberOfParameters();
				}
			}
		}
		return n;
	}
	
	public void setParameter(double[] params, int offset) {
		super.setParameter(params, offset);
		if(refIdx == 0){
			offset += super.getNumberOfParameters();
			for(int i=0;i<e.length;i++){
				if( e[i] != null ) {
					e[i].setParameter(params, offset);
					offset += e[i].getNumberOfParameters();
				}
			}
		} else {
			for(int i=0;i<e.length;i++){
				if( e[i] != null ) {
					e[i].setParameters( ref.e[i] );
				}
			}
		}
	}	
}
