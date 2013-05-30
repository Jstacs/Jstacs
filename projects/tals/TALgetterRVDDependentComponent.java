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

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.DiscreteSequenceEnumerator;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;

/**
 * Class for the RVD-dependent component of TALgetter.
 * @author Annett Wolf, Jan Grau
 *
 */
public class TALgetterRVDDependentComponent extends AbstractDifferentiableStatisticalModel{

	/**
	 * The RVD alphabet
	 */
	protected AlphabetContainer alphabetsRVD;
	private double ess;
	/**
	 * The models for the binding specificities of the RVDs
	 */
	protected HomogeneousMMDiffSM[] hmm_c;
	private boolean isInitialized=false;
	
	/**
	 * Creates a new TALgetter RVD-dependent component.
	 *  
	 * @param alphabets the alphabet of target sites, typically {@link DNAAlphabetContainer}.
	 * @param alphabetsRVD the alphabet of RVDs
	 * @param length the expected length of target sites
	 * @param ess the equivalent sample size
	 * @param priorImp the prior importances, in same order as RVDs in <code>alphabetsRVD</code>
	 * @param priorPrefs the prior binding preferences, in same order as RVDs in <code>alphabetsRVD</code>
	 * @throws Exception if something went wrong
	 */
	public TALgetterRVDDependentComponent(AlphabetContainer alphabets,AlphabetContainer alphabetsRVD, int length, double ess, double[] priorImp, double[][] priorPrefs) throws Exception{
		super(alphabets,1);
		if(alphabets.getAlphabetLengthAt(0)>0&alphabetsRVD.getAlphabetLengthAt(0)>0){
			this.alphabetsRVD=alphabetsRVD;
			this.ess=ess;
			double norm = 0;
			for(int i=0;i<priorImp.length;i++){
				norm += priorImp[i];
			}
			this.hmm_c=new HomogeneousMMDiffSM[getNumberOfSymbols( alphabetsRVD )];
			double[][][] hypi = new double[hmm_c.length][1][(int)alphabets.getAlphabetLengthAt( 0 )];
			double[] priorImSum = new double[hmm_c.length];
			for(int c=0;c<alphabetsRVD.getAlphabetLengthAt( 0 );c++){
				for(int i=0;i<hypi[getMappedIndex( alphabetsRVD, c )][0].length;i++){
					hypi[getMappedIndex( alphabetsRVD, c )][0][i] += priorPrefs[c][i] * ess*priorImp[c]/norm*length;
				}
				priorImSum[getMappedIndex( alphabetsRVD, c )] += priorImp[c];
			}
			for(int c=0;c<hmm_c.length;c++){
				hmm_c[c]=new HomogeneousMMDiffSM(alphabets,0,ess*priorImSum[c]/norm,hypi[c],true,true,1);
			}
		}else{throw new Exception("Alphabet ist nicht OK!");}
	}
	
	/**
	 * Creates a new {@link TALgetterRVDDependentComponent} from its XML description.
	 * @param xml the XML description
	 * @throws NonParsableException if the description could not be parsed.
	 */
	public TALgetterRVDDependentComponent(StringBuffer xml) throws NonParsableException{
		super(xml);
	}
	
	@Override
	public TALgetterRVDDependentComponent clone() throws CloneNotSupportedException {
		TALgetterRVDDependentComponent clone=(TALgetterRVDDependentComponent) super.clone();
		clone.hmm_c=ArrayHandler.clone( hmm_c );
		return clone;
	}
	
	/**
	 * Returns the number of RVDs
	 * @param con the alphabet of RVDs
	 * @return the number of symbols
	 */
	protected int getNumberOfSymbols(AlphabetContainer con){
		return (int)con.getAlphabetLengthAt( 0 );
	}
	
	/**
	 * Returns the mapped index of the symbol with code <code>original</code>.
	 * @param con the RVD alphabet
	 * @param original the code of the original symbol
	 * @return the mapped code
	 * @see TAL_A_NSFMap
	 */
	protected int getMappedIndex(AlphabetContainer con, int original){
		return original;
	}
	
	/**
	 * Returns the index of the symbol at position <code>pos</code> of <code>seq</code>.
	 * @param seq the sequence
	 * @param pos the position
	 * @return the code
	 * @see TAL_A_NSFMap
	 */
	protected int getIndex(Sequence seq, int pos){
		return getMappedIndex(seq.getAlphabetContainer(), seq.discreteVal(pos));
	}
	
	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start)
			throws Exception {
		for(int c=0;c<hmm_c.length;c++){
			hmm_c[c].addGradientOfLogPriorTerm(grad, start);
			start+=hmm_c[c].getNumberOfParameters();
		}
	}

	@Override
	public double getESS() {
		return ess;
	}

	@Override
	public double getInitialClassParam(double classProb) {
		// TODO Auto-generated method stub
		return 0;
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
		double logPrior = 0;
		for(int c=0;c<hmm_c.length;c++){
			logPrior+=hmm_c[c].getLogPriorTerm();
		}
		return logPrior;
	}

	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		return 0;
	}

	@Override
	public boolean isNormalized() {
		return true;
	}


	@Override
	public double[] getCurrentParameterValues() throws Exception {
		double[] params=new double[getNumberOfParameters()];
		int off = 0;
		for(int c=0;c<hmm_c.length;c++){
			System.arraycopy( hmm_c[c].getCurrentParameterValues(), 0, params, off, hmm_c[c].getNumberOfParameters() );
			off += hmm_c[c].getNumberOfParameters();
		}
		return params;
	}

	@Override
	public String getInstanceName() {
		return "TALgetterRVDDependentComponent";
	}

	@Override
	public int getLength() {
		return 1;
	}


	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		ReferenceSequenceAnnotation data_anno =(ReferenceSequenceAnnotation)seq.getSequenceAnnotationByType("reference", 0);
		Sequence rvd_seq=data_anno.getReferenceSequence();
		
		double erg = 0;
		
		
		erg=hmm_c[getIndex( rvd_seq, start-1 )].getLogScoreFor(seq, start, start);
		
		return erg;
	}
	
	/**
	 * Returns the specificities for the RVD at position <code>pos</code> of RVD sequence <code>rvds</code>.
	 * @param rvds the RVD sequence
	 * @param pos the position
	 * @return the specificities
	 */
	public double[] getSpecificities(Sequence rvds, int pos){
		int rvd = getIndex( rvds, pos );
		double[] spec = new double[(int)hmm_c[rvd].getAlphabetContainer().getAlphabetLengthAt( 0 )];
		DiscreteSequenceEnumerator dse = new DiscreteSequenceEnumerator( hmm_c[rvd].getAlphabetContainer(), 1, false );
		int i=0;
		while(dse.hasMoreElements()){
			spec[i] = hmm_c[rvd].getLogScoreFor( dse.nextElement() );
			i++;
		}
		Normalisation.logSumNormalisation( spec );
		
		return spec;
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start,
			IntList indices, DoubleList partialDer) {
		ReferenceSequenceAnnotation data_anno =(ReferenceSequenceAnnotation)seq.getSequenceAnnotationByType("reference", 0);
		Sequence rvd_seq=data_anno.getReferenceSequence();
		
		int alph_length=(int)alphabets.getAlphabetLengthAt(0);
		
		double logScore = 0;
		int index_rvd=getIndex( rvd_seq, start-1 );
		IntList temp = new IntList();
			
		logScore = hmm_c[index_rvd].getLogScoreAndPartialDerivation(seq, start, start, temp, partialDer);
		for(int i=0;i<temp.length();i++){
			indices.add( temp.get( i ) + alph_length*index_rvd );
		}
		
		return logScore;
	}

	@Override
	public int getNumberOfParameters() {
		return hmm_c[0].getNumberOfParameters() * hmm_c.length;
	}

	@Override
	public int getNumberOfRecommendedStarts() {
		// TODO Auto-generated method stub
		return 1;
	}

	@Override
	public void initializeFunction(int index, boolean freeParams,
			DataSet[] data, double[][] weights) throws Exception {
		initializeFunctionRandomly(freeParams);
	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		for(int c=0;c<hmm_c.length;c++){
			hmm_c[c].initializeFunctionRandomly(freeParams);
		}
		isInitialized=true;
	}

	@Override
	public boolean isInitialized() {
		return isInitialized;
	}

	@Override
	public void setParameters(double[] params, int start) {
		for(int c=0;c<hmm_c.length;c++){
			hmm_c[c].setParameters(params, start);
			start += hmm_c[c].getNumberOfParameters();
		}
	}

	
	
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, alphabetsRVD, "alphabetsRVD" );
		XMLParser.appendObjectWithTags( xml, alphabets, "alphabets" );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		XMLParser.appendObjectWithTags( xml, ess, "ess" );
		XMLParser.appendObjectWithTags( xml, hmm_c, "hmmc" );
		XMLParser.appendObjectWithTags( xml, isInitialized, "isInitialized" );
		XMLParser.addTags( xml, "TALANSF" );
		return xml;
	}
	
	@Override
	public String toString(NumberFormat nf){
		StringBuffer sb = new StringBuffer();
		for(int i=0;i<hmm_c.length;i++){
			sb.append( alphabetsRVD.getSymbol( 0, i )+"\t" );
			sb.append( hmm_c[i].toString( nf )+"\n" );
		}
		return sb.toString();
	}

	@Override
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, "TALANSF" );
		alphabetsRVD = (AlphabetContainer)XMLParser.extractObjectForTags( xml, "alphabetsRVD" );
		alphabets = (AlphabetContainer)XMLParser.extractObjectForTags( xml, "alphabets" );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
		hmm_c = (HomogeneousMMDiffSM[])XMLParser.extractObjectForTags( xml, "hmmc" );
		isInitialized = XMLParser.extractObjectForTags( xml, "isInitialized", boolean.class );
	}
}
