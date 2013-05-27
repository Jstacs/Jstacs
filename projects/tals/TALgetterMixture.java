/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 *
 * For more information on Jstacs, visit http://www.jstacs.de
 */
package projects.tals;
import java.text.NumberFormat;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * Class for the mixture component of TALgetter.
 * 
 * @author Annett Wolf, Jan Grau
 *
 */
public class TALgetterMixture extends AbstractDifferentiableStatisticalModel{
	
	private double ess; 
	private boolean isInitialized=false;
	private double[] params;
	private double[] probs;
	private double HyperParams[];
	private double HyperSum[];
	private int p_anz;
	private boolean p_gesamte_seq;
	
	/**
	 * Creates a mixture component.
	 * @param alphabetsRVD the alphabet of RVDs
	 * @param length the expected length of RVD sequences
	 * @param ess the equivalent sample size
	 * @param priorImp the prior importances
	 * @throws Exception if something went wrong
	 */
	public TALgetterMixture(AlphabetContainer alphabetsRVD, int length, double ess, double[] priorImp ) throws Exception{
		super(alphabetsRVD, 1);
		if(alphabetsRVD.getAlphabetLengthAt(0)>0){

			this.ess=ess;
			this.HyperParams=new double[(int)alphabetsRVD.getAlphabetLengthAt(0)+1];
			this.HyperSum=new double[(int)alphabetsRVD.getAlphabetLengthAt(0)+1];

			Arrays.fill(HyperSum, ess*(double)length/(alphabetsRVD.getAlphabetLengthAt(0)));
			for(int i=0;i<priorImp.length;i++){
				HyperParams[i] = HyperSum[i]*priorImp[i];
			}
			HyperParams[HyperParams.length-1] = 0.5*HyperSum[HyperParams.length-1];
			
			this.params=new double[(int)alphabetsRVD.getAlphabetLengthAt(0)+1];
			this.probs=new double[(int)alphabetsRVD.getAlphabetLengthAt(0)+1];
			
		}else{throw new Exception("Alphabet wrong");}
	}
	
	
	
	/**
	 * Creates a new {@link TALgetterMixture} from its XML description.
	 * @param xml the XML description
	 * @throws NonParsableException if the description could not be parsed.
	 */
	public TALgetterMixture( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}



	@Override
	public TALgetterMixture clone() throws CloneNotSupportedException {
		TALgetterMixture clone = (TALgetterMixture) super.clone();
		
		clone.params = params.clone();
		clone.probs = probs.clone();
		clone.HyperSum=HyperSum.clone();
		clone.HyperParams=HyperParams.clone();
		return clone;
	}

	
	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start)
			throws Exception {
        
		for (int j = 0; j < getNumberOfParameters()-1; j++, start++) {
			grad[start] += HyperParams[j] - HyperSum[j] * probs[j];
		}
		
	}

	@Override
	public double getESS() {
		return ess;
	}

	@Override
	public double getLogNormalizationConstant() {
		return 0;
	}

	@Override
	public double getLogPartialNormalizationConstant(int parameterIndex)
			throws Exception {
		return Double.NEGATIVE_INFINITY;
	}

	@Override
	public double getLogPriorTerm() {
		double[] logBeta=new double[params.length];
		double logPrior=0;
		for(int c=0;c<params.length-1;c++){
			logBeta[c]=( Gamma.logOfGamma( HyperParams[c] ) + Gamma.logOfGamma( HyperSum[c] - HyperParams[c] ) ) - Gamma.logOfGamma( HyperSum[c] );
			logPrior+= HyperParams[c]*Math.log( probs[c] ) + (HyperSum[c] - HyperParams[c])*Math.log1p( -probs[c] );
			logPrior -= logBeta[c];
		}
		return logPrior;
	}

	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		return 0;
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		return params.clone();
	}

	@Override
	public String getInstanceName() {
		return "TALgetterMixture";
	}

	/**
	 * Returns the importance of the RVD at position <code>pos</code>
	 * in <code>rvds</code>.
	 * @param rvds the RVD sequence
	 * @param pos the position
	 * @return the importance
	 */
	public double getImportance(Sequence rvds, int pos){
		return probs[rvds.discreteVal( pos )];
	}
	
	@Override
	public double getLogScoreFor(Sequence seq, int start) {

		ReferenceSequenceAnnotation data_anno =(ReferenceSequenceAnnotation)seq.getSequenceAnnotationByType("reference", 0);
		Sequence rvd_seq=data_anno.getReferenceSequence();
		
		double erg = 0;
		int p_ab=seq.getLength()-p_anz;
		
		
		if((start>=p_ab)||p_gesamte_seq){
			if(p_gesamte_seq){
				erg=Math.log(probs[rvd_seq.discreteVal(start-1)])+(start-1)*Math.log(probs[probs.length-1]);
			}
			else{
				erg=Math.log(probs[rvd_seq.discreteVal(start-1)])+(start-(seq.getLength()-1-p_anz))*Math.log(probs[probs.length-1]);
			}
		}else{
			erg=Math.log(probs[rvd_seq.discreteVal(start-1)]);
		}
		return erg;
	}
	

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start,
			IntList indices, DoubleList partialDer) {
		ReferenceSequenceAnnotation data_anno =(ReferenceSequenceAnnotation)seq.getSequenceAnnotationByType("reference", 0);
		Sequence rvd_seq=data_anno.getReferenceSequence();
		
		
		double logScore = 0;
		int index_rvd=rvd_seq.discreteVal(start-1);
		
		int p_ab=seq.getLength()-p_anz;;
		
		if((start>=p_ab)||p_gesamte_seq){
			if(p_gesamte_seq){
				logScore=Math.log(probs[index_rvd])+(start-1)*Math.log(probs[probs.length-1]);
				indices.add(probs.length-1);
				partialDer.add((start-1)*(1-probs[probs.length-1]));
			}else{
				logScore=Math.log(probs[index_rvd])+(start-(seq.getLength()-1-p_anz))*Math.log(probs[probs.length-1]);
				indices.add(probs.length-1);
				partialDer.add((start-(seq.getLength()-1-p_anz))*(1-probs[probs.length-1]));
			}
		}else{
			logScore=Math.log(probs[index_rvd]);
		}
		indices.add(index_rvd);
		partialDer.add(1-probs[index_rvd]);
		
		
		
		return logScore;
	}

	@Override
	public int getNumberOfParameters() {
		return params.length;
	}

	@Override
	public void initializeFunction(int index, boolean freeParams,
			DataSet[] data, double[][] weights) throws Exception {
		initializeFunctionRandomly(freeParams);
		
	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		
		double temp;
		for(int c=0;c<params.length;c++){

			
			temp=Math.random();
			
			params[c] = Math.log(temp);
			probs[c]=temp/(temp+1);
			
		}
		isInitialized=true;
	}

	@Override
	public boolean isInitialized() {
		return isInitialized;
	}

	@Override
	public void setParameters(double[] params, int start) {
		for (int j = 0; j < getNumberOfParameters(); j++, start++) {
			this.params[j] = params[start];
			probs[j] = Math.exp(this.params[j])/(1+Math.exp(this.params[j]));
		}
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, alphabets, "alphabets" );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		XMLParser.appendObjectWithTags( xml, ess, "ess" );
		XMLParser.appendObjectWithTags( xml, HyperParams, "hyperParams" );
		XMLParser.appendObjectWithTags( xml, HyperSum, "hyperSum" );
		XMLParser.appendObjectWithTags( xml, isInitialized, "isInitialized" );
		XMLParser.appendObjectWithTags( xml, p_anz, "panz" );
		XMLParser.appendObjectWithTags( xml, p_gesamte_seq,"pgesseq" );
		XMLParser.appendObjectWithTags( xml, params, "params" );
		XMLParser.appendObjectWithTags( xml, probs, "probs" );
		XMLParser.addTags( xml, "TALMSF" );
		return xml;
	}
	
	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, "TALMSF" );
		alphabets = (AlphabetContainer)XMLParser.extractObjectForTags( xml, "alphabets" );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
		HyperParams = (double[])XMLParser.extractObjectForTags( xml, "hyperParams" );
		HyperSum = (double[])XMLParser.extractObjectForTags( xml, "hyperSum" );
		isInitialized = XMLParser.extractObjectForTags( xml, "isInitialized", boolean.class );
		p_anz = XMLParser.extractObjectForTags( xml, "panz", int.class );
		p_gesamte_seq = XMLParser.extractObjectForTags( xml, "pgesseq", boolean.class );
		params = (double[])XMLParser.extractObjectForTags( xml, "params" );
		probs = (double[])XMLParser.extractObjectForTags( xml, "probs" );
	}

	@Override
	public String toString(NumberFormat nf){
		StringBuffer sb = new StringBuffer();
		for(int i=0;i<probs.length-1;i++){
			sb.append( alphabets.getSymbol( 0, i )+"\t"+nf.format(probs[i])+"\n" );
		}
		return sb.toString();
	}
	
	
}
