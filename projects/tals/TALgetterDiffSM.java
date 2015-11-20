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

import java.io.PrintWriter;
import java.text.NumberFormat;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.DiscreteSequenceEnumerator;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.DifferentiableStatisticalModelWrapperTrainSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.ToolBox;



/**
 * Class for the statistical model of TALgetter.
 * 
 * @author Annett Wolf, Jan Grau
 *
 */
public class TALgetterDiffSM extends AbstractDifferentiableStatisticalModel{
	
	private TALgetterRVDDependentComponent tal_A_NSF;
	private HomogeneousMMDiffSM tal_U_NSF;
	private TALgetterMixture tal_M_NSF;
	private HomogeneousMMDiffSM FirstPosHMM;
	private boolean isInitialized=true;
	private int Ordnung_tal_U;
	private double Ess;

	private boolean isFixed;
	private double[][][][] condProbs;
	private double[] firstPosProbs;
	private int[] firstPosMatches;
	private int[][][][] matches;
	private int[] pows;
	
	/**
	 * Creates a new TALgetter model.
	 *  
	 * @param alphabets the alphabet of target sites, typically {@link DNAAlphabetContainer}.
	 * @param alphabetsRVD the alphabet of RVDs
	 * @param midLength the expected length of target sites
	 * @param ess the equivalent sample size
	 * @param order_talU the order of the RVD-independent component
	 * @param priorFP the prior probabilities for position 0
	 * @param priorImp the prior importances, in same order as RVDs in <code>alphabetsRVD</code>
	 * @param priorPrefs the prior binding preferences, in same order as RVDs in <code>alphabetsRVD</code>
	 * @throws Exception if something went wrong
	 */
	public TALgetterDiffSM(AlphabetContainer alphabets,AlphabetContainer alphabetsRVD, double midLength, double ess, int order_talU, double[] priorFP, double[] priorImp, double[][] priorPrefs) throws Exception{
		super(alphabets, 0);
		this.Ess=ess;
		priorFP = priorFP.clone();
		for(int i=0;i<priorFP.length;i++){
			priorFP[i] *= ess;
		}
		double part = 0;
		for(int i=0;i<priorImp.length;i++){
			part += priorImp[i];
		}
		part /= priorImp.length;

		this.FirstPosHMM=new HomogeneousMMDiffSM(alphabets,0,ess,new double[][]{priorFP},true,true,1);
		this.tal_A_NSF=getTalANsf(alphabets, alphabetsRVD, (int)midLength, ess, part, priorImp,priorPrefs);

		this.tal_U_NSF=new HomogeneousMMDiffSM(alphabets, order_talU, ess*(1-part), (int)midLength);
		this.tal_M_NSF=new TALgetterMixture(alphabetsRVD, (int)midLength, ess,priorImp);
		this.Ordnung_tal_U=order_talU;
		this.pows = new int[order_talU];
		for(int i=0;i<pows.length;i++){
			pows[i] = (int)Math.pow( (int)alphabets.getAlphabetLengthAt( 0 ), i );
		}
	} 

	/**
	 * Returns the RVD dependent component.
	 * 
	 * @param alphabets the alphabet of target sites, typically {@link DNAAlphabetContainer}.
	 * @param alphabetsRVD the alphabet of RVDs
	 * @param midLength the expected lenght of target sites
	 * @param ess the equivalent sample size
	 * @param part the part of the ess for the RVD-dependent component
	 * @param priorImp the prior importances, in same order as RVDs in <code>alphabetsRVD</code>
	 * @param priorPrefs the prior binding preferences, in same order as RVDs in <code>alphabetsRVD</code>
	 * @return the RVD-dependent component
	 * @throws Exception if something went wrong
	 */
	protected TALgetterRVDDependentComponent getTalANsf(AlphabetContainer alphabets, AlphabetContainer alphabetsRVD, int midLength, double ess, double part, double[] priorImp, double[][] priorPrefs) throws Exception{
		return new TALgetterRVDDependentComponent(alphabets, alphabetsRVD, (int)midLength, ess*part, priorImp,priorPrefs);
	}
	
	/**
	 * Creates a new {@link TALgetterDiffSM} from its XML description.
	 * @param xml the XML description
	 * @throws NonParsableException if the description could not be parsed.
	 */
	public TALgetterDiffSM( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}


	@Override
	public TALgetterDiffSM clone() throws CloneNotSupportedException {
		TALgetterDiffSM clone=(TALgetterDiffSM) super.clone();
		clone.FirstPosHMM=FirstPosHMM.clone();
		clone.tal_A_NSF=tal_A_NSF.clone();
		clone.tal_M_NSF=tal_M_NSF.clone();
		clone.tal_U_NSF=tal_U_NSF.clone();
		clone.pows = pows.clone();
		if(isFixed){
			clone.condProbs = ArrayHandler.clone( condProbs );
			clone.firstPosMatches = firstPosMatches.clone();
			clone.firstPosProbs = firstPosProbs.clone();
			clone.matches = ArrayHandler.clone( matches );
		}
		return clone;
	}

	/**
	 * Returns the RVD alphabet of this {@link TALgetterDiffSM}.
	 * @return the RVD alphabet
	 */
	public AlphabetContainer getRVDAlphabet(){
		return tal_M_NSF.getAlphabetContainer();
	}
	
	

	@Override
	public int getNumberOfRecommendedStarts() {
		return 1;
	}


	@Override
	public void addGradientOfLogPriorTerm(double[] grad, int start)
			throws Exception {
			FirstPosHMM.addGradientOfLogPriorTerm(grad, start);
			start+=FirstPosHMM.getNumberOfParameters();
			tal_A_NSF.addGradientOfLogPriorTerm(grad, start);
			start+=tal_A_NSF.getNumberOfParameters();
			tal_M_NSF.addGradientOfLogPriorTerm(grad, start);
			start+=tal_M_NSF.getNumberOfParameters();
	}

	@Override
	public double getESS() {
		return Ess;
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
		logPrior+=FirstPosHMM.getLogPriorTerm();
		logPrior+=tal_A_NSF.getLogPriorTerm();
		logPrior+=tal_M_NSF.getLogPriorTerm();
		return logPrior;
	}

	@Override
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		return 0;
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		double[] params=new double[getNumberOfParameters()];
		int off = 0;
		
		System.arraycopy(FirstPosHMM.getCurrentParameterValues(), 0, params, off, FirstPosHMM.getNumberOfParameters());
		off+=FirstPosHMM.getNumberOfParameters();
		System.arraycopy( tal_A_NSF.getCurrentParameterValues(), 0, params, off, tal_A_NSF.getNumberOfParameters() );
		off += tal_A_NSF.getNumberOfParameters();
		System.arraycopy( tal_M_NSF.getCurrentParameterValues(), 0, params, off, tal_M_NSF.getNumberOfParameters() );
		return params;
	}

	@Override
	public String getInstanceName() {
		return "TAL_Finder_NSF";
	}

	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		return getLogScoreFor(seq,start,seq.getLength()-1);
	}
	
	public double getBestPossibleScore(Sequence tal, double[] scs){
		return getBestPossibleScore( tal, scs, null );
	}
	
	/**
	 * Returns the score of the best possible target site for RVDs <code>tal</code>
	 * and fills the partial scores for each position into <code>scs</code>.
	 * 
	 * @param tal the RVD sequence
	 * @param scs the array for partial scores, may be null
	 * @return the best possible score
	 */
	public double getBestPossibleScore(Sequence tal, double[] scs, int[] bestSeq){
		if(isFixed){
			double sum = 0;
			
			//int[] bestSeq = new int[tal.getLength()+1];
			
			double temp = ToolBox.max( firstPosProbs );
			int idx = ToolBox.getMaxIndex( firstPosProbs );
			if(bestSeq != null){
				bestSeq[0] = idx;
			}
			sum += temp;
			if(scs != null){
				scs[0] = temp;
			}
			
			int context = 0;//TODO
			for(int i=1;i<=tal.getLength();i++){
				
				int rvd = tal.discreteVal( i-1 );
				int order = Math.min( i-1, Ordnung_tal_U );
				
				double max = Double.NEGATIVE_INFINITY;
				for(int j=0;j<condProbs[order][rvd].length;j++){
					temp = ToolBox.max( condProbs[order][rvd][j] );
					if(temp > max){
						max = temp;
						idx = ToolBox.getMaxIndex( condProbs[order][rvd][j] );
						if(bestSeq != null){
							bestSeq[i] = idx;
						}
					}
				}
				sum += max;
				if(scs != null){
					scs[i] = max;
				}
				
			}
			
			return sum;
			
		}else{
			throw new RuntimeException();
		}
	}
	
	/**
	 * Returns the partial score for RVDs <code>rvds</code> and a sequence represented by ints.
	 * @param rvds the RVD sequence
	 * @param seq the sequence as ints
	 * @param start the start position
	 * @param off the offset to start partial score from
	 * @param length the length of the sub-sequence for which partial score is computed
	 * @return the partial score
	 */
	public double getPartialLogScoreFor(Sequence rvds, int[] seq, int start, int off, int length) {
		if(isFixed){
			
			int initStart = start;
			int end = (start)+off+length-1;
			
			double logScore = 0;
			if(off == 0){
				logScore += firstPosProbs[seq[start]];
				start++;
			}

			int context = 0;
			for(int i=start+off;i<=end;i++){
				int rvd = rvds.discreteVal( i-1-initStart );
				int order = Math.min( i-start + (off == 0 ? 0 : -1), Ordnung_tal_U );
			//	System.out.println(order);
				if(order > 0){
					if(order == 1){
						context = seq[i-1];
					}else if(order < Ordnung_tal_U || i-(start+off-1)<order){
						context = getContextIndexFor( pows, seq, i-order, order );
					}else{
						context = getContextIndexFor( pows, seq, i-order, order, context);
					}
				}
				//System.out.println(i+" "+context);
				int curr = seq[i];
				//System.out.println(((DiscreteAlphabet)rvds.getAlphabetContainer().getAlphabetAt( 0 )).getSymbolAt( rvd )+" "+curr);
				logScore += condProbs[order][rvd][context][curr];
			}
			
			return logScore;
		}else{
			throw new RuntimeException();
		}
	}
	
	/**
	 * Returns the partial score for RVDs <code>rvds</code> and a sequence.
	 * @param rvds the RVD sequence
	 * @param seq the sequence
	 * @param start the start position
	 * @param off the offset to start partial score from
	 * @param length the length of the sub-sequence for which partial score is computed
	 * @return the partial score
	 */
	public double getPartialLogScoreFor(Sequence rvds, Sequence seq, int start, int off, int length) {
		if(isFixed){
			
			int initStart = start;
			int end = (start)+off+length-1;
			
			double logScore = 0;
			if(off == 0){
				logScore += firstPosProbs[seq.discreteVal( start )];
				start++;
			}

			int context = 0;
			for(int i=start+off;i<=end;i++){
				int rvd = rvds.discreteVal( i-1-initStart );
				int order = Math.min( i-start + (off == 0 ? 0 : -1), Ordnung_tal_U );
			//	System.out.println(order);
				if(order > 0){
					if(order == 1){
						context = seq.discreteVal( i-1 );
					}else if(order < Ordnung_tal_U || i-(start+off-1)<order){
						context = getContextIndexFor( pows, seq, i-order, order );
					}else{
						context = getContextIndexFor( pows, seq, i-order, order, context);
					}
				}
				//System.out.println(i+" "+context);
				int curr = seq.discreteVal( i );
				//System.out.println(((DiscreteAlphabet)rvds.getAlphabetContainer().getAlphabetAt( 0 )).getSymbolAt( rvd )+" "+curr);
				logScore += condProbs[order][rvd][context][curr];
			}
			
			return logScore;
		}else{
			throw new RuntimeException();
		}
	}
	
	@Override
	public double getLogScoreFor(Sequence seq, int start, int end) {
		double logScore=0;
		
		if(isFixed){
			
			logScore += firstPosProbs[seq.discreteVal( start )];
			//System.out.println(start+" "+logScore);
			start++;
			Sequence rvds = ((ReferenceSequenceAnnotation) seq.getSequenceAnnotationByType( ReferenceSequenceAnnotation.TYPE, 0 )).getReferenceSequence();
			
			int context = 0;//TODO
			for(int i=start;i<=end;i++){
				int rvd = rvds.discreteVal( i-1 );
				int order = Math.min( i-start, Ordnung_tal_U );
				//int context = 0;//TODO
				if(order > 0){
					if(order == 1){
						context = seq.discreteVal( i-1 );
					}else if(order < Ordnung_tal_U){
						context = getContextIndexFor( pows, seq, i-order, order );
					}else{
						context = getContextIndexFor( pows, seq, i-order, order, context);
					}
				}
				int curr = seq.discreteVal( i );
				//System.out.println(i+" "+condProbs[order][rvd][context][curr]+" "+order+" "+context+" "+curr);
				logScore += condProbs[order][rvd][context][curr];
			}
			
		}else{

			logScore=FirstPosHMM.getLogScoreFor(seq,start,start);
			start=start+1;


			double logScore_tal_M;
			double[] temp=new double[2];



			for(int i=start;i<=end;i++){

				logScore_tal_M = tal_M_NSF.getLogScoreFor( seq, i );


				temp[0] = logScore_tal_M + tal_A_NSF.getLogScoreFor( seq, i );
				//temp[1] = Math.log1p( -Math.exp(logScore_tal_M) );
				temp[1] = Math.log( 1.0-Math.exp(logScore_tal_M) );
				if(Ordnung_tal_U == 0 || i-start == 0){
					temp[1] += tal_U_NSF.getLogScoreFor( seq, i, i );
				}else{
					int l = Math.min( i-start, Ordnung_tal_U ) + 1;
					temp[1] += tal_U_NSF.getLogScoreFor( seq, i-l+1, i ) - tal_U_NSF.getLogScoreFor( seq, i-l+1, i-1 );
					//System.out.println((i-l+1)+" "+(i-1));
				}
				logScore += Normalisation.getLogSum( temp );
			}
		}
		return logScore;
	}
	
	

	/**
	 * Fixes the model, i.e., fixes all parameters and pre-computes several lookup tables
	 * to accelerate computation of scores.
	 * 
	 * @throws WrongAlphabetException if the alphabet did not fit
	 * @throws WrongSequenceTypeException if the alphabet did not fit an {@link IntSequence}
	 */
	public void fix() throws WrongAlphabetException, WrongSequenceTypeException{
		DiscreteAlphabet rvds = (DiscreteAlphabet)tal_M_NSF.getAlphabetContainer().getAlphabetAt( 0 );
		DiscreteAlphabet nucs = (DiscreteAlphabet)alphabets.getAlphabetAt( 0 );
		
		firstPosProbs = new double[(int)nucs.length()];
		firstPosMatches = new int[firstPosProbs.length];
		for(int i=0;i<nucs.length();i++){
			firstPosProbs[i] = FirstPosHMM.getLogScoreFor( new IntSequence( alphabets, i ) );
			if(firstPosProbs[i] > -Math.log( nucs.length() )){
				firstPosMatches[i] = 1;
			}
		}
		
		SequenceAnnotation[] rvdAnns = new SequenceAnnotation[(int)rvds.length()];
		for(int i=0;i<rvds.length();i++){
			rvdAnns[i] = new ReferenceSequenceAnnotation( "rvds", new IntSequence( getRVDAlphabet(), i ) );
		}
		
		int order = Ordnung_tal_U;
		
		double[] parts = new double[2];
		
		
		
		condProbs = new double[order+1][][][];
		matches = new int[order+1][][][];
		for(int i=0;i<order+1;i++){
			condProbs[i] = new double[(int)rvds.length()][(int)Math.pow( nucs.length(), i )][(int)nucs.length()];
			matches[i] = new int[(int)rvds.length()][(int)Math.pow( nucs.length(), i )][(int)nucs.length()];
			for(int j=0;j<rvdAnns.length;j++){
				for(int k=0;k<condProbs[i][j].length;k++){
					for(int l=0;l<condProbs[i][j][k].length;l++){
						Sequence temp = getSequenceFor(k,l,i);
						//System.out.println(i+" "+k+" "+l+" "+temp);
						temp = temp.annotate( false, rvdAnns[j] );
						Sequence temp2 = temp;
						if(temp.getLength()>2){
							temp2 = temp.getSubSequence( temp.getLength()-2 );
							temp2 = temp2.annotate( false, temp.getAnnotation() );
						}
						double logScore_tal_M = tal_M_NSF.getLogScoreFor( temp, 1 );
						
						parts[0] = logScore_tal_M + tal_A_NSF.getLogScoreFor( temp2, 1 );
						parts[1] = Math.log( 1.0-Math.exp(logScore_tal_M) );
						
						if(i == 0){
							parts[1] += tal_U_NSF.getLogScoreFor( temp, temp.getLength()-1, temp.getLength()-1 );
						}else{
							parts[1] += tal_U_NSF.getLogScoreFor( temp, 1, temp.getLength()-1 ) - tal_U_NSF.getLogScoreFor( temp, 1, temp.getLength()-2 );
						}
						if(parts[0] > parts[1] && parts[0] - logScore_tal_M > -Math.log( nucs.length() )){
							matches[i][j][k][l] = 1;
						}
						condProbs[i][j][k][l] = Normalisation.getLogSum( parts );
					}
				}
			}
		}
		isFixed = true;
	}
	
	/**
	 * Returns the best-matching RVDs for a given input sequence, where <code>allowed</code> indicates,
	 * which RVDs of the complete alphabet are allowed.
	 * @param target the input sequence
	 * @param allowed the allowed RVDs
	 * @return the best RVD sequence
	 * @throws WrongAlphabetException if the alphabet did not fit
	 * @throws WrongSequenceTypeException if the alphabet did not fit
	 */
	public Sequence getBestRVDsFor(Sequence target, int[] allowed) throws WrongAlphabetException, WrongSequenceTypeException{
		if(allowed == null){
			allowed = new int[(int)getRVDAlphabet().getAlphabetLengthAt( 0 )];
			for(int i=0;i<allowed.length;i++){
				allowed[i] = i;
			}
		}
		int[] tal = new int[target.getLength()-1];
		if(isFixed){
			
			int start = 1;
			int end = target.getLength()-1;
			
			int context = 0;//TODO
			for(int i=start;i<=end;i++){
				
				int order = Math.min( i-start, Ordnung_tal_U );
				//int context = 0;//TODO
				if(order > 0){
					if(order == 1){
						context = target.discreteVal( i-1 );
					}else if(order < Ordnung_tal_U){
						context = getContextIndexFor( pows, target, i-order, order );
					}else{
						context = getContextIndexFor( pows, target, i-order, order, context);
					}
				}
				int curr = target.discreteVal( i );
				//System.out.println(i+" "+condProbs[order][rvd][context][curr]+" "+order+" "+context+" "+curr);
				double max = Double.NEGATIVE_INFINITY;
				for(int j=0;j<allowed.length;j++){
					if(condProbs[order][allowed[j]][context][curr] > max){
						max = condProbs[order][allowed[j]][context][curr];
						tal[i-1] = allowed[j];
					}
				}
			}
			return new IntSequence( getRVDAlphabet(), tal );
			
		}else{
			throw new RuntimeException();
		}
	}
	
	private static final int getContextIndexFor( int[] pows, Sequence seq, int start, int order, int context ) {
		if(order > 1){
			//System.out.println(context + " "+pows[order-1]+" "+pows[1]);
			//context %= pows[order-1];
			//context *= pows[1];
			context = (context % pows[order-1])*pows[1] + seq.discreteVal( start+order-1 );
		}else{
			context = seq.discreteVal( start+order-1 );
		}
		return context;
	}
	
	private static final int getContextIndexFor( int[] pows, int[] seq, int start, int order, int context ) {
		if(order > 1){
			//System.out.println(context + " "+pows[order-1]+" "+pows[1]);
			//context %= pows[order-1];
			//context *= pows[1];
			context = (context % pows[order-1])*pows[1] + seq[ start+order-1 ];
		}else{
			context = seq[ start+order-1 ];
		}
		return context;
	}
	
	private static final int getContextIndexFor(int[] pows, Sequence seq, int start, int order){
		int index = 0;
		for(int i=0;i<order;i++){
			index += pows[order-i-1]*seq.discreteVal( start+i );
		}
		return index;
	}
	
	private static final int getContextIndexFor(int[] pows, int[] seq, int start, int order){
		int index = 0;
		for(int i=0;i<order;i++){
			index += pows[order-i-1]*seq[ start+i ];
		}
		return index;
	}
	
	private Sequence getSequenceFor(int context, int curr, int order) throws WrongAlphabetException, WrongSequenceTypeException{
		int[] cont = new int[order+2];
		cont[cont.length-1] = curr;
		for(int i=0;i<order;i++){
			cont[i+1] = context / pows[order-i-1];
			context = context % pows[order-i-1];
		}
		return new IntSequence( alphabets, cont );
	}
	
	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start,
			IntList indices, DoubleList partialDer) {
		return getLogScoreAndPartialDerivation( seq, start, seq.getLength()-1, indices, partialDer );
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, int end,
			IntList indices, DoubleList partialDer) {

		if(isFixed){
			throw new RuntimeException( "Model is fixed" );
		}
		
		double logScore=0;
		
		logScore=FirstPosHMM.getLogScoreAndPartialDerivation( seq.getSubSequence(start, 1), 0,indices, partialDer);
		start=start+1;
		
		
		double logScore_tal_M, logScore_tal_A, logScore_tal_U;
		double[] temp=new double[2];
		double localLogScore;
		DoubleList partialA = new DoubleList();
		DoubleList partialM = new DoubleList();
		DoubleList partialUPos = new DoubleList();
		DoubleList partialUNeg = new DoubleList();
		IntList tempIndices = new IntList();

		int offA=FirstPosHMM.getNumberOfParameters();
		int offM = offA+tal_A_NSF.getNumberOfParameters();
		int offU = offM + tal_M_NSF.getNumberOfParameters();
		
		
			
		for(int i=start;i<=end;i++){
			partialA.clear(); partialM.clear(); partialUPos.clear(); partialUNeg.clear();
			logScore_tal_A = tal_A_NSF.getLogScoreAndPartialDerivation( seq, i,tempIndices, partialA );
			
			for(int j=0;j<tempIndices.length();j++){
				indices.add( tempIndices.get( j )+offA );
			}
			tempIndices.clear();
			
			
		
				
			logScore_tal_M = tal_M_NSF.getLogScoreAndPartialDerivation( seq, i,tempIndices, partialM );
			
			
			temp[0] = logScore_tal_M + logScore_tal_A;
			
			for(int j=0;j<tempIndices.length();j++){
				indices.add( tempIndices.get( j )+offM );
			}
			tempIndices.clear();
			
			logScore_tal_U = 0;
			if(Ordnung_tal_U == 0 || i-start == 0){
				logScore_tal_U += tal_U_NSF.getLogScoreFor( seq, i, i );
			}else{
				int l = Math.min( i-start, Ordnung_tal_U ) + 1;
				logScore_tal_U += tal_U_NSF.getLogScoreFor( seq, i-l+1, i ) - tal_U_NSF.getLogScoreFor( seq, i-l+1, i-1 );
			}
			temp[1] = Math.log1p( -Math.exp(logScore_tal_M) ) + logScore_tal_U;
			
			
			localLogScore = Normalisation.getLogSum( temp );
			logScore += localLogScore;
			
			for(int j=0;j<partialA.length();j++){
				partialDer.add( partialA.get( j ) * Math.exp( logScore_tal_A + logScore_tal_M - localLogScore ) );
			}
			
			for(int j=0;j<partialM.length();j++){
				partialDer.add( partialM.get( j ) * Math.exp( logScore_tal_M - localLogScore )*(Math.exp( logScore_tal_A ) - Math.exp( logScore_tal_U )) );
			}
			
			for(int j=0;j<partialUPos.length();j++){
				partialDer.add( partialUPos.get( j ) * Math.exp( logScore_tal_U - localLogScore )*(1.0-Math.exp( logScore_tal_M )) );
			}
		}
		
		return logScore;
	}

	@Override
	public int getNumberOfParameters() {
		int number=FirstPosHMM.getNumberOfParameters();
		number+=tal_A_NSF.getNumberOfParameters();
		number+=tal_M_NSF.getNumberOfParameters();
		return number;
	}

	@Override
	public void initializeFunction(int index, boolean freeParams,
			DataSet[] data, double[][] weights) throws Exception {
		initializeFunctionRandomly(freeParams);
		
	}
	
	/**
	 * Trains the RVD-independent component
	 * @param trainDs the training data
	 * @throws Exception if the training went wrong
	 */
	public void trainIndependent(DataSet trainDs) throws Exception{
		DifferentiableStatisticalModelWrapperTrainSM temp = new DifferentiableStatisticalModelWrapperTrainSM( tal_U_NSF, 1, Optimizer.QUASI_NEWTON_BFGS, new SmallDifferenceOfFunctionEvaluationsCondition( 1E-12 ), 1E-12, 1E-4 );
		temp.setOutputStream( SafeOutputStream.getSafeOutputStream( null ) );
		temp.train( trainDs );
		tal_U_NSF = (HomogeneousMMDiffSM)temp.getFunction();
	}

	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		FirstPosHMM.initializeFunctionRandomly(freeParams);
		
		tal_A_NSF.initializeFunctionRandomly(freeParams);
			
		tal_M_NSF.initializeFunctionRandomly(freeParams);

		isInitialized=true;
		
	}

	@Override
	public boolean isInitialized() {
		return isInitialized;
	}

	@Override
	public void setParameters(double[] params, int start) {
		if(isFixed){
			throw new RuntimeException( "Model is fixed" );
		}
		FirstPosHMM.setParameters(params, start);
		start+=FirstPosHMM.getNumberOfParameters();
		tal_A_NSF.setParameters(params, start);
		start+=tal_A_NSF.getNumberOfParameters();
		tal_M_NSF.setParameters(params, start);
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, alphabets, "alphabets" );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		XMLParser.appendObjectWithTags( xml, Ess, "ess" );
		XMLParser.appendObjectWithTags( xml, FirstPosHMM, "firstPosHMM" );
		XMLParser.appendObjectWithTags( xml, isInitialized, "isInitialized" );
		XMLParser.appendObjectWithTags( xml, Ordnung_tal_U, "order" );
		XMLParser.appendObjectWithTags( xml, tal_A_NSF, "talansf" );
		XMLParser.appendObjectWithTags( xml, tal_M_NSF, "talmnsf" );
		XMLParser.appendObjectWithTags( xml, tal_U_NSF, "talunsf" );
		XMLParser.appendObjectWithTags( xml, isFixed, "isFixed" );
		XMLParser.addTags( xml, "TALFinder" );
		return xml;
		
	}
	
	@Override
	protected void fromXML(StringBuffer xml) throws de.jstacs.io.NonParsableException {
		xml = XMLParser.extractForTag( xml, "TALFinder" );
		alphabets = (AlphabetContainer)XMLParser.extractObjectForTags( xml, "alphabets" );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		Ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
		FirstPosHMM = (HomogeneousMMDiffSM)XMLParser.extractObjectForTags( xml, "firstPosHMM" );
		isInitialized = XMLParser.extractObjectForTags( xml, "isInitialized", boolean.class );
		Ordnung_tal_U = XMLParser.extractObjectForTags( xml, "order", int.class );
		tal_A_NSF = (TALgetterRVDDependentComponent)XMLParser.extractObjectForTags( xml, "talansf" );
		tal_M_NSF = (TALgetterMixture)XMLParser.extractObjectForTags( xml, "talmnsf" );
		tal_U_NSF = (HomogeneousMMDiffSM)XMLParser.extractObjectForTags( xml, "talunsf" );
		isFixed = XMLParser.extractObjectForTags( xml, "isFixed", boolean.class );
		this.pows = new int[Ordnung_tal_U];
		for(int i=0;i<pows.length;i++){
			pows[i] = (int)Math.pow( (int) alphabets.getAlphabetLengthAt( 0 ), i );
		}
		if(isFixed){
			try {
				fix();
			} catch (Exception e ) {
				NonParsableException ex = new NonParsableException( e.getMessage() );
				ex.setStackTrace( e.getStackTrace() );
				throw ex;
			}
		}
		
	}
	
	/**
	 * Returns for each position of an RVD sequence the corresponding specificities and importances.
	 * @param rvds the sequence of RVDs
	 * @return the specificities and importances
	 */
	public Pair<double[][],double[]> getSpecificitiesAndImportances(Sequence rvds){
		double[][] specs = new double[rvds.getLength()+1][];
		
		specs[0] = new double[(int)FirstPosHMM.getAlphabetContainer().getAlphabetLengthAt( 0 )];
		DiscreteSequenceEnumerator dse = new DiscreteSequenceEnumerator( FirstPosHMM.getAlphabetContainer(), 1, false );
		int i=0;
		while(dse.hasMoreElements()){
			specs[0][i] = FirstPosHMM.getLogScoreFor( dse.nextElement() );
			i++;
		}
		Normalisation.logSumNormalisation( specs[0] );
		for(i=0;i<rvds.getLength();i++){
			specs[i+1] = tal_A_NSF.getSpecificities( rvds, i );
		}
		
		double[] imp = new double[rvds.getLength()];
		for(i=0;i<rvds.getLength();i++){
			imp[i] = tal_M_NSF.getImportance( rvds, i );
		}
		
		return new Pair<double[][],double[]>(specs,imp);
	}
	
	/**
	 * Sets the initial probabilities of the RVD-independent component to the stationary distribution.
	 */
	public void setIndependentToStationary(){
		tal_U_NSF.setStartParamsToConditionalStationaryDistributions();
	}
	
	public void addAndSet(String[] rvds, double[][] specs, double[] fpSpec) throws IllegalArgumentException, WrongAlphabetException, DoubleSymbolException{
		int oldLen = (int)getRVDAlphabet().getAlphabetAt( 0 ).length();
		AlphabetContainer con = tal_A_NSF.addAndSet( rvds, specs );
		int newLen = (int)con.getAlphabetAt( 0 ).length();
		tal_M_NSF.addAndSet(con,rvds);
		if(fpSpec != null){
			double[] temp = fpSpec.clone();
			for(int i=0;i<temp.length;i++){
				temp[i] = Math.log( temp[i] );
			}
			FirstPosHMM.setParameters( temp, 0 );
		}
	}
	
	@Override
	public String toString(NumberFormat nf){
		StringBuffer buf = new StringBuffer();
		buf.append( "******************************************************" );
		buf.append( "\nModel position 0:\n" );
		buf.append( FirstPosHMM.toString(nf) );
		buf.append( "******************************************************" );
		buf.append( "\nModel mixture:\n" );
		buf.append( tal_M_NSF.toString(nf) );
		buf.append( "******************************************************" );
		buf.append( "\nModel dependent:\n" );
		buf.append( tal_A_NSF.toString(nf) );
		buf.append( "******************************************************" );
		buf.append( "\nModel context:" );
		buf.append( tal_U_NSF.toString(nf) );
		buf.append( "******************************************************" );
		return buf.toString();
	}

	/**
	 * Returns a string representing matches between <code>seq</code> and an annotated RVD sequence
	 * @param seq the sequence
	 * @return the matches
	 */
	public String getMatchString( Sequence seq ){

		ReferenceSequenceAnnotation data_anno =(ReferenceSequenceAnnotation)seq.getSequenceAnnotationByType("reference", 0);
		Sequence rvds=data_anno.getReferenceSequence();
		
		return getMatchString( rvds, seq );
	}
	
	/**
	 * Returns a string representing matches between <code>seq</code> and <code>rvds</code>
	 * @param rvds the RVD sequence
	 * @param seq the sequence
	 * @return the matches
	 */
	public String getMatchString( Sequence rvds, Sequence seq ){
		
		seq = seq.annotate( true, new ReferenceSequenceAnnotation( "seq", rvds ) );
		
		StringBuffer sb = new StringBuffer();
	
		if( FirstPosHMM.getLogProbFor( seq, 0, 0 ) > - Math.log( seq.getAlphabetContainer().getAlphabetLengthAt( 0 ) ) ){
			sb.append( "M" );
		}else{
			sb.append( "m" );
		}


		double logScore_tal_M;
		double[] temp=new double[2];

		int start = 1;
		int end = rvds.getLength();

		for(int i=start;i<=end;i++){

			logScore_tal_M = tal_M_NSF.getLogScoreFor( seq, i );


			temp[0] = logScore_tal_M + tal_A_NSF.getLogScoreFor( seq, i );
			//temp[1] = Math.log1p( -Math.exp(logScore_tal_M) );
			temp[1] = Math.log( 1.0-Math.exp(logScore_tal_M) );
			if(Ordnung_tal_U == 0 || i-start == 0){
				temp[1] += tal_U_NSF.getLogScoreFor( seq, i, i );
			}else{
				int l = Math.min( i-start, Ordnung_tal_U ) + 1;
				temp[1] += tal_U_NSF.getLogScoreFor( seq, i-l+1, i ) - tal_U_NSF.getLogScoreFor( seq, i-l+1, i-1 );
			}
			if(logScore_tal_M > Math.log( 0.7 ) && temp[0] > temp[1] && temp[0] - logScore_tal_M > Math.log( 0.7 )){
				sb.append( "|" );
			}else if( temp[0] - logScore_tal_M > -Math.log( seq.getAlphabetContainer().getAlphabetLengthAt( i ) ) ){
				sb.append( ":" );
			}else{
				sb.append( "x" );
			}
		}
		
		return sb.toString();
	}
	
	/**
	 * Counts the number of matches between <code>seq</code> and <code>rvds</code>.
	 * @param seq the sequence
	 * @param rvds the RVD sequence
	 * @return the number of matches
	 */
	public int getNumberOfMatches( Sequence seq, Sequence rvds ) {
		if(firstPosMatches == null){
			throw new RuntimeException("Only after fix");
		}
	
		int numMatches = firstPosMatches[seq.discreteVal( 0 )];
		
		for(int i=1;i<seq.getLength();i++){
			int rvd = rvds.discreteVal( i-1 );
			int order = Math.min( i, Ordnung_tal_U );
			int context = 0;
			if(order > 0){
				context = getContextIndexFor( pows, seq, i-order, order );
			}
			int curr = seq.discreteVal( i );
			numMatches += matches[order][rvd][context][curr];
		}
		return numMatches;
	}

	/**
	 * Returns the order of the RVD-independent component.
	 * @return the order
	 */
	public int getOrder() {
		return Ordnung_tal_U;
	}

}
