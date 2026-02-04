/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Locale;
import java.util.Random;

import javax.naming.OperationNotSupportedException;

import de.jstacs.algorithms.optimization.NumericalDifferentiableFunction;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.algorithms.optimization.termination.CombinedCondition;
import de.jstacs.algorithms.optimization.termination.IterationCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.algorithms.optimization.termination.TimeCondition;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LogGenDisMixFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.StorableResult;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.SamplingDifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.DifferentiableStatisticalModelWrapperTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.DifferentiableState;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.SimpleDifferentiableState;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableCombinedWrapperEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.DifferentiableSMWrapperEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.UniformEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.filter.Filter;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.BaumWelchParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.HMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.MaxHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.NumericalHMMTrainingParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.NumericalHMMTrainingParameterSet.TrainingType;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.ViterbiParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.DifferentiableTransition;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements.TransitionElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;


/**
 * This class combines an {@link HigherOrderHMM} and a {@link de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel} by implementing some of the declared methods. 
 * 
 * @author Jens Keilwagen
 */
public class DifferentiableHigherOrderHMM extends HigherOrderHMM implements SamplingDifferentiableStatisticalModel {	
	
	/**
	 * The number of parameters of this HMM
	 */
	protected int numberOfParameters;
	
	/**
	 * The equivalent sample size used for the prior
	 */
	protected double ess;
	
	/**
	 * Index array used for computing the gradient
	 */
	protected int[][] index;
	
	/**
	 * Help array for the gradient. (Dimensions: first = layer mod 2, second = context, third = parameter)
	 */
	protected double[][][] gradient, gradient2;
	protected double[][] slice, slice2;
	
	protected double[] logScore;
	protected double[] prop;
	
	/**
	 * Help array for the indexes of the parameters of the states
	 */
	protected IntList[] indicesState; 
	/**
	 * Help array for the indexes of the parameters of the transition
	 */
	protected IntList[] indicesTransition;
	/**
	 * Help array for the derivatives of the parameters of the states
	 */
	protected DoubleList[] partDerState; 
	/**
	 * Help array for the derivatives of the parameters of the transition
	 */
	protected DoubleList[] partDerTransition;
	
	private NumericalHMMTrainingParameterSet.TrainingType training;
	private Type score = Type.LIKELIHOOD;
	private boolean train;

	private double[] forwardIntermediate;
	private IntList childrenBW, childrenFW;
	/**
	 * This is the main constructor.
	 * 
	 * @param trainingParameterSet the {@link de.jstacs.parameters.ParameterSet} that determines the training algorithm and contains the necessary {@link de.jstacs.parameters.Parameter}s
	 * @param name the names of the states
	 * @param emissionIdx the indices of the emissions that should be used for each state, if <code>null</code> state <code>i</code> will use emission <code>i</code>
	 * @param forward a boolean array that indicates whether the symbol on the forward or the reverse complementary strand should be used,
	 * 				  if <code>null</code> all states use the forward strand
	 * @param emission the emissions
	 * @param ess the ess of the model
	 * @param te the {@link TransitionElement}s used for creating a {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.Transition}
	 * 
	 * @throws Exception if 
	 * 	<ul>
	 *  <li>some component could not be cloned</li> 
	 *  <li>some the length of <code>name, emissionIdx,</code> or <code>forward</code> is not equal to the number of states</li>
	 *  <li>not all emissions use the same {@link de.jstacs.data.AlphabetContainer}</li>
	 *  <li>the states can not be handled by the transition
	 *  </ul>
	 */
	public DifferentiableHigherOrderHMM( MaxHMMTrainingParameterSet trainingParameterSet, String[] name, int[] emissionIdx, boolean[] forward,
			DifferentiableEmission[] emission, double ess, TransitionElement... te ) throws Exception {
		this(null, null, trainingParameterSet, name, null, emissionIdx, forward, emission, ess, null, te);
	}
	
	public DifferentiableHigherOrderHMM( String type, int[][] statesGroups, MaxHMMTrainingParameterSet trainingParameterSet, String[] name, Filter[] filter, int[] emissionIdx, boolean[] forward,
			DifferentiableEmission[] emission, double ess, int[] transIndex, TransitionElement... te ) throws Exception {
		super( type, statesGroups, trainingParameterSet, name, filter, emissionIdx, forward, emission, transIndex, te );
		getOffsets();
		if( ess < 0 ) {
			throw new IllegalArgumentException();
		}
		this.ess = ess;
		childrenBW = new IntList();
		childrenFW = new IntList();
		forwardIntermediate  = new double[getNumberOfStates()];
		provideSlice();
	}
	
	protected void provideSlice() {
		int dim = transition.getNumberOfIndexes( transition.getMaximalMarkovOrder() );		
		slice = new double[2][dim];
		slice2 = new double[2][dim];
	}
	
	public void setTrainingParameter( HMMTrainingParameterSet trainingParameterSet ) throws CloneNotSupportedException {
		super.setTrainingParameter(trainingParameterSet);
		if( trainingParameter instanceof NumericalHMMTrainingParameterSet ) {
			training = ((NumericalHMMTrainingParameterSet) trainingParameter).getTrainingType();
		} else if( trainingParameter instanceof ViterbiParameterSet ){
			training = TrainingType.VITERBI;
		} else if( trainingParameter instanceof BaumWelchParameterSet ){
			training = TrainingType.LIKELIHOOD;
		} else {
			if( !(trainingParameterSet instanceof MaxHMMTrainingParameterSet) ) {
				throw new IllegalArgumentException( "Training parameter set must be a MaxHMMTrainingParameterSet.");
			}
		}
		setTrain(false);
	}	
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link DifferentiableHigherOrderHMM} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DifferentiableHigherOrderHMM} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public DifferentiableHigherOrderHMM( StringBuffer xml ) throws NonParsableException {
		super( xml );
		getOffsets();
		
		childrenBW = new IntList();
		childrenFW = new IntList();
		forwardIntermediate  = new double[getNumberOfStates()];
		provideSlice();
	}
	
	@Override
	protected void appendFurtherInformation( StringBuffer xml ) {
		super.appendFurtherInformation( xml );
		XMLParser.appendObjectWithTags( xml, ess, "ess" );
	}

	@Override
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		super.extractFurtherInformation( xml );
		ess = XMLParser.extractObjectForTags( xml, "ess", double.class );
	}

	protected void createHelperVariables() {
		if( container == null ) {
			if( numberOfParameters == 0 ) {
				getOffsets();
			}
			int maxOrder = transition.getMaximalMarkovOrder(), anz = 0, i;
			for( i = 0; i <= maxOrder; i++ ) {
				anz = Math.max( anz, transition.getNumberOfIndexes( i ) );
			}
			if( gradient == null || gradient[0].length != anz || gradient[0][0].length != numberOfParameters ) {
				gradient = new double[2][anz][numberOfParameters];
				gradient2 = new double[2][anz][numberOfParameters];
				
				index = new int[4][anz];
			}
			if( indicesState == null ) {
				anz = transition.getMaximalNumberOfChildren();
				try {
					indicesState = ArrayHandler.createArrayOf( new IntList(), states.length );
					partDerState = ArrayHandler.createArrayOf( new DoubleList(), states.length );
					
					indicesTransition = ArrayHandler.createArrayOf( new IntList(), anz );
					partDerTransition = ArrayHandler.createArrayOf( new DoubleList(), anz );
				} catch( CloneNotSupportedException cnse ) {
					throw getRunTimeException( cnse );
				}
			}
			this.logScore = new double[2];
			this.prop = new double[2];
		}
		super.createHelperVariables();
	}
	
	protected void createStates() {
		this.states = new SimpleDifferentiableState[emissionIdx.length];
		for( int i = 0; i < emissionIdx.length; i++ ) {
			this.states[i] = new SimpleDifferentiableState( (DifferentiableEmission) emission[emissionIdx[i]], name[i], forward[i] );
		}
	}

	public DifferentiableHigherOrderHMM clone() throws CloneNotSupportedException {
		//prepare for clone
		double[][][] grad = gradient;
		gradient = null;
		IntList[] ind = indicesState;
		indicesState = null;
		
		//clone
		DifferentiableHigherOrderHMM clone = (DifferentiableHigherOrderHMM) super.clone();
		clone.forwardIntermediate = forwardIntermediate.clone();
		clone.childrenBW = childrenBW.clone();
		clone.childrenFW = childrenFW.clone();
		clone.provideSlice();
				
		//reverse
		gradient = grad;
		indicesState = ind;
		return clone;
	}
	
	public double getESS() {
		return ess;
	}

	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception {
		for( int e = 0; e < emission.length; e++ ) {
			((DifferentiableEmission)emission[e]).addGradientOfLogPriorTerm( grad, start );
		}
		((DifferentiableTransition) transition).addGradientForLogPriorTerm( grad, start );
	}
	
	private void getOffsets() {
		numberOfParameters = 0;
		for( int e = 0; e < emission.length; e++ ) {
			numberOfParameters = ((DifferentiableEmission)emission[e]).setParameterOffset( numberOfParameters );
			if( numberOfParameters == UNKNOWN ) {
				return;
			}
		}
		numberOfParameters = ((DifferentiableTransition)transition).setParameterOffset( numberOfParameters );
		if( numberOfParameters == UNKNOWN ) {
			return;
		}
		createHelperVariables();
	}
	
	public int getNumberOfParameters() {
		return numberOfParameters;
	}

	public int getNumberOfRecommendedStarts() {
		return trainingParameter.getNumberOfStarts();
	}

	public double[] getCurrentParameterValues() throws Exception {
		int n = getNumberOfParameters();
		
		if( n != UNKNOWN ) {
			double[] params = new double[n];
			for( int e = 0; e < emission.length; e++ ) {
				((DifferentiableEmission)emission[e]).fillCurrentParameter( params );
			}
			((DifferentiableTransition)transition).fillParameters( params );
			return params;
		} else {
			throw new IllegalArgumentException();
		}
	}
	
	public boolean isInitialized() {
		return true;
	}
/*	
	public void setTrainingParameters(MaxHMMTrainingParameterSet params) throws CloneNotSupportedException{
		this.trainingParameter = (HMMTrainingParameterSet)params.clone();
	}
*/
	public void setParameters( double[] params, int start ) {
		for( int e = 0; e < emission.length; e++ ) {
			((DifferentiableEmission)emission[e]).setParameter( params, start );
		}
		((DifferentiableTransition)transition).setParameters( params, start );
	}
	
	public void initializeFunctionRandomly( boolean freeParams ) throws Exception {
		if(skipInit){
			return;
		}
		initializeRandomly();
		getOffsets();
	}
	
	public void initializeTransitionRandomly() throws Exception {
		transition.initializeRandomly();
	}
	
	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception {
		if(skipInit){
			return;
		}
		//XXX
		if( trainingParameter instanceof NumericalHMMTrainingParameterSet ) {
			boolean[] diffEM = new boolean[emission.length]; 
			for( int i = 0; i < emission.length; i++ ) {
				diffEM[i] = emission[i] instanceof DifferentiableSMWrapperEmission || emission[i] instanceof DifferentiableCombinedWrapperEmission;
			}
						
			try {
				AbstractTerminationCondition tc = new CombinedCondition(3, 
						new IterationCondition(25),
						new SmallDifferenceOfFunctionEvaluationsCondition(1E-3),
						new TimeCondition(3600)
					);
				MaxHMMTrainingParameterSet maxPar = new ViterbiParameterSet(1, tc, ((NumericalHMMTrainingParameterSet) trainingParameter).getNumberOfThreads() );
				
				//create simple variant
				DifferentiableEmission[] dEmission = new DifferentiableEmission[emission.length];
				ArrayList<Sequence>[] seqs = new ArrayList[emission.length];
				DoubleList[] subW = new DoubleList[emission.length];
				int[] l = new int[emission.length], offset=new int[l.length];
				for( int i = 0; i < emission.length; i++ ) {
					if( diffEM[i] ) {
/*							if( getAlphabetContainer().isDiscrete() ) {
								dEmission[i] = new DiscreteEmission(getAlphabetContainer(),getESS());
							} else {/**/
							dEmission[i] = new UniformEmission(getAlphabetContainer());
//							}
						seqs[i] = new ArrayList<Sequence>();
						subW[i] = new DoubleList();
						if( emission[i] instanceof DifferentiableSMWrapperEmission ) {
							DifferentiableSMWrapperEmission help = (DifferentiableSMWrapperEmission) emission[i];
							l[i] = help.getLength();
							offset[i] = help.getOffset(); 
						} else {
							DifferentiableCombinedWrapperEmission help = (DifferentiableCombinedWrapperEmission) emission[i];
							l[i] = help.getLength();
							offset[i] = help.getOffset();
						}
					} else {
						dEmission[i] = (DifferentiableEmission) emission[i];
					}
				}
				DifferentiableHigherOrderHMM simple;
				//maxPar.getParameterForName( "termination condition" ).setValue(tc);
				if( this instanceof FastDifferentiableHigherOrderHMM ) {
					simple = new FastDifferentiableHigherOrderHMM(type, statesGroups, maxPar, name, filter, emissionIdx, dEmission, ess, transIndex, getTransitionElements() );
				} else {
					simple = new DifferentiableHigherOrderHMM(type, statesGroups, maxPar, name, filter, emissionIdx, forward, dEmission, ess, transIndex, getTransitionElements() );
				}
				
				//train generatively without motifs
				simple.train(data[index], weights==null? null : weights[index] );
				
				//train discriminatively without motifs XXX
				/*
				transition.setParameters(simple.transition);
				for( int i = 0; i < dEmission.length; i++ ) {
					if( !diffEM[i] ) {
						dEmission[i].setParameters(simple.emission[i]);
					}
				}
				if( this instanceof FastDifferentiableHigherOrderHMM ) {
					simple = new FastDifferentiableHigherOrderHMM(type, statesGroups, (NumericalHMMTrainingParameterSet) trainingParameter, name, filter, emissionIdx, dEmission, ess, transIndex, getTransitionElements() );
				} else {
					simple = new DifferentiableHigherOrderHMM(type, statesGroups, (NumericalHMMTrainingParameterSet) trainingParameter, name, filter, emissionIdx, forward, dEmission, ess, transIndex, getTransitionElements() );
				}
				simple.setSkiptInit(true);
				simple.train(data[index], weights==null? null : weights[index] );
				/**/
				
				//ML estimation of (previously not used) models
				if( type!=null ) {
					Random random = new Random();
					for( int i = 0; i < data[index].getNumberOfElements(); i++ ) {
						Sequence seq = data[index].getElementAt(i);
						double w = (weights==null || weights[index]==null) ? 1: weights[index][i];
						SequenceAnnotation sa = seq.getSequenceAnnotationByType(type, 0);
						if( sa !=null ) {
							//select one of all possible paths
							StorableResult res = (StorableResult) sa.getResultAt(random.nextInt(sa.getNumberOfResults()) );
							int[] allowedStatesGroup=((AllowedStatesGroups) res.getResultInstance()).groups;
							for( int j = 0; j < allowedStatesGroup.length; j++ ) {
								int[] em = statesGroupsEm[allowedStatesGroup[j]];
								int k = 
									//0; ///XXX simple solution
									random.nextInt(em.length);
								int e = em[k];
								int st = j+offset[e];
								if( seqs[e] != null 
									&& st>=0 && st+l[e]<seq.getLength() ) {
									seqs[e].add( seq.getSubSequence(st, l[e]) );
									subW[e].add(w);
								}
							}
						}
					}
				}
				/**/
				
				//set
				transition.setParameters(simple.transition);
				NumberFormat nf = NumberFormat.getInstance(Locale.US);
				nf.setMaximumFractionDigits(3);
				for( int i = 0; i < emission.length; i++ ) {
					if( diffEM[i] ) {
						if( seqs[i].size()==0 ) {
							emission[i].initializeFunctionRandomly();
							//((DifferentiableSMWrapperEmission)emission[i]).initializeUniformly();
						} else {
System.out.println("initialize emission " + i + " with data");
							DataSet partial = new DataSet("",seqs[i]);
System.out.println("#="+seqs[i].size() + "\tlength=" + partial.getElementLength() + "\tfirst=" + seqs[i].get(0));
							if( emission[i] instanceof DifferentiableSMWrapperEmission ) {
								((DifferentiableSMWrapperEmission)emission[i]).initializeFunction( 0, false, new DataSet[] { partial }, new double[][] {subW[i].toArray()} );
							} else {
								((DifferentiableCombinedWrapperEmission)emission[i]).initializeFunction( 0, false, new DataSet[] { partial }, new double[][] {subW[i].toArray()} );
							}
System.out.println(emission[i].toString(nf));
						}/**/
					} else {
						emission[i].setParameters(simple.emission[i]);
					}
				}
			} catch( Exception e ) {
				e.printStackTrace();
				sostream.writeln("Problem while initialization from data. " + e.getClass().getSimpleName() + ": " + e.getCause() );
				initializeFunctionRandomly( freeParams );
			}
		} else {			
			initializeFunctionRandomly( freeParams );
		}
	}
	
	public void train( DataSet data, double[] weights ) throws Exception {
		train( data, weights,((MaxHMMTrainingParameterSet) trainingParameter).getTerminationCondition() );
	}
	
	public void train( DataSet data, double[] weights, AbstractTerminationCondition tc ) throws Exception {
		if( trainingParameter instanceof NumericalHMMTrainingParameterSet ) {
			setTrain(true);
			System.out.println(training + " threads=" + threads );
			NumericalHMMTrainingParameterSet params = (NumericalHMMTrainingParameterSet) trainingParameter;
			LogPrior p = DoesNothingLogPrior.defaultInstance;
					//new SimpleGaussianSumLogPrior(1);//TODO
			DifferentiableStatisticalModelWrapperTrainSM model = new DifferentiableStatisticalModelWrapperTrainSM( this, params.getNumberOfThreads(), params.getAlgorithm(), tc, params.getLineEps(), params.getStartDistance(), params.randomInitialization(), p, true );
			model.setOutputStream( sostream );
			model.train( data, weights );
			
			DifferentiableHigherOrderHMM hmm = (DifferentiableHigherOrderHMM) model.getFunction();
			this.emission = hmm.emission;
			createStates();
			this.transition = hmm.transition;
			setTrain(false);
		} else {
			super.train( data, weights );
		}
	}
	
	public void testGrad( DataSet d ) throws Exception {
		double[][] wei = new double[1][d.getNumberOfElements()];
		Arrays.fill(wei[0], 1);
		LogGenDisMixFunction f = new LogGenDisMixFunction( threads, new DifferentiableSequenceScore[] {this}, new DataSet[]{d}, wei, DoesNothingLogPrior.defaultInstance, LearningPrinciple.getBeta(LearningPrinciple.ML), true, false );
		f.reset();
		double[] x = f.getParameters(KindOfParameter.PLUGIN);
		NumericalDifferentiableFunction.compare(f, x, 1e-3);
	}
	

//XXX is normalized? start
	public boolean isNormalized() {
		boolean isNormalized=((DifferentiableTransition)transition).isNormalized();
System.out.println("trans\t" + isNormalized);
		if( isNormalized ) {
			int i = 0;
			while( i < emission.length && (isNormalized=((DifferentiableEmission)emission[i]).isNormalized()) ) {
				System.out.println("em\t"+i+"/"+emission.length+"\t" + isNormalized);
				i++;
			}
		}
		return isNormalized;
	}
	
	public double getLogNormalizationConstant() {
		return 0;
	}

	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception {
		return Double.NEGATIVE_INFINITY;
	}
	
	public double getInitialClassParam(double classProb) {
		return Math.log( classProb );
	}
//end

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#getLogScoreFor(de.jstacs.data.Sequence)
	 */
	@Override
	public double getLogScoreFor( Sequence seq ) {
		return getLogScoreFor( seq, 0 );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#getLogScoreFor(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor( Sequence seq, int start ) {
		return getLogScoreFor( seq, start, seq.getLength()-1 );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#getLogScoreFor(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor( Sequence seq, int start, int end ) {
		return getLogScoreAndPartialDerivation(seq, start, end, null, null);
	}	
	
	public double getLogScoreAndPartialDerivation( Sequence seq, IntList indices, DoubleList partialDer ) {
		return getLogScoreAndPartialDerivation( seq, 0, indices, partialDer );
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int startPos, IntList indices, DoubleList partialDer ) {
		return getLogScoreAndPartialDerivation( seq, startPos, seq.getLength()-1, indices, partialDer);
	}
		
	public double getLogScoreAndPartialDerivation( Sequence seq, int startPos, int endPos, IntList indices, DoubleList partialDer ) {
		if( !train || type==null ) {
			return compute( seq, null, startPos, endPos, indices, partialDer, score, false );
		} else {
			SequenceAnnotation sa = seq.getSequenceAnnotationByType(type, 0);
			if( sa !=null ) {
				//determine best
				Type s = training.isViterbiLike()?Type.VITERBI:Type.LIKELIHOOD;
				double max = Double.NEGATIVE_INFINITY;
				int[] best = null;
				if( sa.getNumberOfResults()==1 ) {
					best = ((AllowedStatesGroups) ((StorableResult) sa.getResultAt(0)).getResultInstance()).groups;
				} else {
					int[] allowedStatesGroup=null;
					for( int i = 0; i <sa.getNumberOfResults(); i++ ) {
						StorableResult res = (StorableResult) sa.getResultAt(i);
						AllowedStatesGroups asg = (AllowedStatesGroups) res.getResultInstance();
						allowedStatesGroup=asg.groups;
						double current = compute( seq, allowedStatesGroup, startPos, endPos, null, null, s, false );
						if( current > max ) {
							max = current;
							best=allowedStatesGroup;
						}
					}
				}
				//compute score for best
				return compute( seq, best, startPos, endPos, indices, partialDer, s, training.isDiscrimnative() );
			} else {
				throw new IllegalArgumentException("Missing allowedStatesGroup");
			}
		}
	}
		
	/*
	//TODO To be removed
	private double[] toGrad( IntList ind, DoubleList partD ) {
		double[] grad = new double[numberOfParameters];
		for( int i = 0; i < ind.length(); i++ ) {
			grad[ind.get(i)] += partD.get(i);
		}
		return(grad);
	}
	private void gradComp( IntList oldIndices, DoubleList oldPartDer, IntList newIndices, DoubleList newPartDer, boolean showAll ) {
		System.out.println(oldIndices.length() + "\t" +oldPartDer.length() + "\t" + newIndices.length() + "\t" + newPartDer.length());
		double[] oldGrad = toGrad( oldIndices, oldPartDer );
		double[] newGrad = toGrad( newIndices, newPartDer );
		double max = 0;
		for( int i = 0; i < numberOfParameters; i++ ) {
			if( showAll ) System.out.println(i + "\t" +oldGrad[i] + "\t" + newGrad[i]);
			double delta = Math.abs( oldGrad[i] - newGrad[i] );
			max = Math.max(delta, max);
		}
		System.out.println("delta=" + max);
		System.out.println();
	}
	
	public void myTest( Sequence seq ) throws Exception {
		setTrain(true);
		Type s = Type.VITERBI;
		
		boolean gradOut=false;
		int start=0;
		int end = seq.getLength()-1;
		int[] asg = null;
		IntList newIndices = new IntList();
		DoubleList newPartDer = new DoubleList();
		IntList oldIndices = new IntList();
		DoubleList oldPartDer = new DoubleList();
		
		//generative
		fillBwdOrViterbiMatrix( s, start, end, 0, seq, asg, false );
		System.out.println(
				bwdMatrix[0][0]
				+ "\t"+ logScoreAndPartialDerivation(seq, start, end, oldIndices, oldPartDer, asg, s, 1)
				+ "\t" + compute( seq, asg, start, end, null, null, s, false )
				+ "\t" + compute( seq, asg, start, end, newIndices, newPartDer, s, false )
		);
		gradComp(oldIndices, oldPartDer, newIndices, newPartDer, gradOut);
				
		//generative
		SequenceAnnotation sa = seq.getSequenceAnnotationByType(type, 0);
		asg = ((AllowedStatesGroups) ((StorableResult) sa.getResultAt(0)).getResultInstance()).groups;
		newIndices.clear();
		newPartDer.clear();
		oldIndices.clear();
		oldPartDer.clear();
		fillBwdOrViterbiMatrix( Type.VITERBI, start, end, 0, seq, asg, false );
		System.out.println(
				bwdMatrix[0][0]
				+ "\t"+ logScoreAndPartialDerivation(seq, start, end, oldIndices, oldPartDer, asg,Type.VITERBI, 1)
				+ "\t" + compute( seq, asg, start, end, null, null, s, false )
				+ "\t" + compute( seq, asg, start, end, newIndices, newPartDer, s, false )
		);
		gradComp(oldIndices, oldPartDer, newIndices, newPartDer, gradOut);

		//Discriminative
		newIndices.clear();
		newPartDer.clear();
		oldIndices.clear();
		oldPartDer.clear();
		System.out.println(
				getLogScoreFor(seq, start, end)
				+ "\t" + getLogScoreAndPartialDerivation(seq, oldIndices, oldPartDer)
				+ "\t" + compute( seq, asg, 0, seq.getLength()-1, null, null, s, true )
				+ "\t" + compute( seq, asg, 0, seq.getLength()-1, newIndices, newPartDer, s, true )
		);
		gradComp(oldIndices, oldPartDer, newIndices, newPartDer, gradOut);
	}/**/

	//compute score (and if asked for the partial derivatives)
	private double compute( Sequence seq, int[] allowedStatesGroup, int startPos, int endPos, IntList indices, DoubleList partialDer, Type s, boolean discriminative ) {
		try {
			//slice => Viterbi/Likelihood in allowed states
			//slice2 => Likelihood
			Type t = discriminative ? Type.LIKELIHOOD : null;

			int maxOrder = transition.getMaximalMarkovOrder();
			boolean zero = maxOrder == 0;
			int l = endPos-startPos+1, stateID, context, n;
			
			DifferentiableTransition diffTransition = (DifferentiableTransition) transition;
			
			//init
			double val, emTrans;
			for( int i = 0; i < slice.length; i++ ) {
				Arrays.fill( slice[i], 0 );
				Arrays.fill( slice2[i], 0 );
			}
			boolean computePartialDer = indices != null && partialDer != null;
			double[][][] g = null, g2 = null;
			if( computePartialDer ) {
				g=gradient;
				g2=gradient2;
				for( int idx2 = 0; idx2 < gradient[1].length; idx2++ ) {
					Arrays.fill( gradient[0][idx2], 0 );
					Arrays.fill( gradient[1][idx2], 0 );
					Arrays.fill( gradient2[0][idx2], 0 );
					Arrays.fill( gradient2[1][idx2], 0 );
				}
				for( stateID = 0; stateID < states.length; stateID++ ) {
					indicesState[stateID].clear();
					partDerState[stateID].clear();
				}
			}
			fillFilter(endPos, seq);
			int[] all = getAllowedContext(endPos+1, startPos, null, maxOrder);
			int[] allowed = getAllowedContext(endPos+1, startPos, allowedStatesGroup, maxOrder);
			int a=allowed.length-1;
	
			for( int x = all.length-1; x>=0; x-- ) {
				context = all[x];
				n = transition.getNumberOfChildren( l, context );
				boolean allowedState = a >=0 ? context == allowed[a] : false;
	
				//the stateIds which are filled
				childrenBW.clear();
				childrenFW.clear();
				if( zero || finalState[transition.getLastContextState( l, context )] ) {
					val = 0;
				} else {
					val = Double.NEGATIVE_INFINITY;
				}
				//for all different children states
				for( stateID = 0; stateID < n; stateID++ ) {
					transition.fillTransitionInformation( l, context, stateID, container );
					if( filterRes[container[0]] ) {
						if( states[container[0]].isSilent() ) {
							emTrans=0;//there is no emission (silent state)
							if( computePartialDer ) {
								indicesTransition[stateID].clear();
								partDerTransition[stateID].clear();
								emTrans += diffTransition.getLogScoreAndPartialDerivation( l, context, stateID, indicesTransition[stateID], partDerTransition[stateID], seq, endPos );
							} else {
								emTrans += diffTransition.getLogScoreFor( l, context, stateID, seq, endPos );
							}
							
							if( allowedState ) {
								backwardIntermediate[stateID] = emTrans
									+ slice[l%2][container[1]]; //viterbi score until next position
								childrenBW.add(stateID);
							}
							if( discriminative ) {
								forwardIntermediate[stateID] = emTrans
									+ slice2[l%2][container[1]]; //score until next position
								childrenFW.add(stateID);
							}
							
							index[0][stateID] = container[0];
							index[1][stateID] = container[1];
							index[2][stateID] = container[2];
							index[3][stateID] = 0;
						}
					}
				}
				slice[l%2][context] = merge( backwardIntermediate, childrenBW, l, context, val, s, g );
				slice2[l%2][context] = merge( forwardIntermediate, childrenFW, l, context, val, t, g2 );
				
				if( allowedState ) {
					a--;
				}
			}
			//System.out.println( seq.toString(endPos, endPos+1) + "\t" + l + "\t" + Arrays.toString( slice[l%2] ) );
			
			//compute scores for all positions backward
			while( --l >= 0 ) {
				fillLogEmissionAndPartialDer( endPos, seq, computePartialDer );
				fillFilter(endPos, seq);
				
				//for all different contexts
				all = getAllowedContext(endPos, startPos, null, maxOrder);
				allowed = getAllowedContext(endPos, startPos, allowedStatesGroup, maxOrder);
				a=allowed.length-1;
				for( int x = all.length-1; x>=0; x-- ) {
					context = all[x];
					n = transition.getNumberOfChildren( l, context );
					boolean allowedState = a >=0 ? context == allowed[a] : false;
					
					//for all different children states
					childrenBW.clear();
					childrenFW.clear();
					for( stateID = 0; stateID < n; stateID++ ) {
						transition.fillTransitionInformation( l, context, stateID, container );
						if( filterRes[container[0]] ) {
							emTrans = logEmission[getIndex(container[0])]; //emission
							if( computePartialDer ) {
								indicesTransition[stateID].clear();
								partDerTransition[stateID].clear();
								emTrans += diffTransition.getLogScoreAndPartialDerivation( l, context, stateID, indicesTransition[stateID], partDerTransition[stateID], seq, endPos );
							} else {
								emTrans += diffTransition.getLogScoreFor( l, context, stateID, seq, endPos );
							}
	
							index[0][stateID] = container[0];
							index[1][stateID] = container[1];
							index[2][stateID] = container[2];
							index[3][stateID] = 0;
							
							if( allowedState ) {
								backwardIntermediate[stateID] = emTrans
									+ slice[(l+container[2])%2][container[1]]; //inside viterbi score until next position
								childrenBW.add(stateID);
							}
							
							if( discriminative ) {
								forwardIntermediate[stateID] = emTrans
										+ slice2[(l+container[2])%2][container[1]];
								childrenFW.add(stateID);	
							}							
						}
					}
					slice[l%2][context]=merge( backwardIntermediate, childrenBW, l, context, Double.NEGATIVE_INFINITY, s, g );
					slice2[l%2][context]=merge( forwardIntermediate, childrenFW, l, context, Double.NEGATIVE_INFINITY, t, g2 );
	
					if( allowedState ) {
						a--;
					}
				}
				
	
				//System.out.println("X\t" + l + "\t" + (l==0?" ":seq.toString(endPos, endPos+1)) + "\t" + allowedStatesGroup[endPos] + "\t" + Arrays.toString( slice[l%2] ) );
				//System.out.println("Y\t" + l + "\t" + (l==0?" ":seq.toString(endPos, endPos+1)) + "\t" + allowedStatesGroup[endPos] + "\t" + Arrays.toString( slice2[l%2] ) );
				
				endPos--;
			}
			logScore[0] = slice[0][0];
			double ls = 0;
			if( discriminative ) {
				ls = slice2[0][0];//Normalisation.logSumNormalisation(slice2[0], 0, slice2[0].length, slice2[1], 0);
				//System.out.println(slice[0][0] + "\t" + ls);
			}
			if( computePartialDer ) {
				/*if( discriminative ) {
					Arrays.fill(gradient2[1][0], 0);
					for( int k = 0; k < gradient2[1].length; k++ ) {
						if( slice2[1][k]!= 0 ) {
							for( int p = 0; p < numberOfParameters; p++ ) {
								gradient2[1][0][p] += slice2[1][k] * gradient2[0][k][p];
							}
						}
					}
				}/**/
				
				for( int p = 0; p < numberOfParameters; p++ ) {
					double v = gradient[0][0][p] - (discriminative ? gradient2[0][0][p] : 0);
					if( v!= 0 ) {
						indices.add( p );
						partialDer.add( v );
					}
				}
			}
//System.out.println(slice[0][0] + "\t" + ls);
			return slice[0][0] - ls;
		} catch( Exception e ) {
			throw getRunTimeException( e );
		}
	}
	
	public double test( Sequence seq, int[] allowedStatesGroup, int startPos, int endPos ) {
		return compute( seq, allowedStatesGroup, startPos, endPos, null, null, Type.VITERBI, true );
	}
	
	//this is slow if emissons are used several times
	protected void fillLogEmissionAndPartialDer( int endPos, Sequence seq, boolean grad ) throws OperationNotSupportedException, WrongLengthException {
		if( grad ) {
			for( int stateID = 0; stateID < states.length; stateID++ ) {
				indicesState[stateID].clear();
				partDerState[stateID].clear();
				logEmission[stateID] = ((DifferentiableState) states[stateID]).getLogScoreAndPartialDerivation( endPos, endPos, indicesState[stateID], partDerState[stateID], seq );
			}
		} else {
			for( int stateID = 0; stateID < states.length; stateID++ ) {
				logEmission[stateID] = states[stateID].getLogScoreFor(endPos, endPos, seq);
			}
		}
	}

	protected double max( double[] intermediate, IntList children, int layer, int context, double val, double[][][] gradient, double[][][] altGradient ) {
		if( children.length() == 0 )  {
			if( gradient != null ) resetGradient( gradient, layer, context, 0 );
			return val;
		} else {
			int idx = ToolBox.getMaxIndex( children, intermediate );
			if( gradient != null ) {
				int h = layer % 2;
				
				double[] old;
				int x = (layer+index[2][idx]) % 2;
				if( index[3][idx]==0 ) {
					old = gradient[x][index[1][idx]];
				} else {
					old = altGradient[x][index[1][idx]];
				}
				for( int p = 0; p < numberOfParameters; p++ ) {
					gradient[h][context][p] = old[p];
				}
				
				miniMerge(idx, 1, h, context, gradient);
			}
			return intermediate[idx];
		}
	}
	
	protected double merge( double[] intermediate, IntList children, int layer, int context, double val, Type score, double[][][] gradient ) {
		if( score==null ) {
			return Double.NEGATIVE_INFINITY;
		}
		if( children.length() == 0 )  {
			if( gradient != null ) resetGradient( gradient, layer, context, 0 );
			return val;
		} else {
			int h = layer % 2;
			if( score == Type.VITERBI ) {
				//determine best 
				int idx = ToolBox.getMaxIndex( children, intermediate );
				//sum gradient
				if( gradient != null ) {
					System.arraycopy( gradient[(layer+index[2][idx]) % 2][index[1][idx]], 0, gradient[h][context], 0, numberOfParameters );
					miniMerge( idx, 1, h, context, gradient );				
				}
				//set partial log probability for (viterbi) path
				return intermediate[idx];
			} else { //LIKELIHOOD
				double res = Normalisation.logSumNormalisation( intermediate, children, intermediate );
				
				if( gradient != null ) {
					// S = \sum_i u_i v_i
					// \frac{\partial \log(S)}{\partial \lambda}
					// = \sum_i 
					//		(u_i v_i / S) \frac{\partial \log u_i}{\partial \lambda}   (AAA)
					//	 	+ (u_i v_i / S) \frac{\partial \log v_i}{\partial \lambda} (BBB)
					
					//help[0][0][i] = u_i v_i / S
					
					Arrays.fill( gradient[h][context], 0 );
					
					//slow version
					/*
					// old = (AAA)
					for( int p = 0; p < numberOfParameters; p++ ) {
						for( int i = 0; i < anz; i++ ) {
							int x = (layer+index[2][i]) % 2;
							gradient[h][context][p] += help[0][0][i] * gradient[x][index[1][i]][p];
						}
					}	
					// transition & emission = (BBB)
					for( int i = 0; i < anz; i++ ) {
						miniMerge( i, help[0][0][i], h, context );
					}
					*/
					
					//fast version
					for( int j = 0; j < children.length(); j++ ) {
						int i = children.get(j);
						// old = (AAA)
						int x = (layer+index[2][i]) % 2;
						for( int p = 0; p < numberOfParameters; p++ ) {
							gradient[h][context][p] += intermediate[i] * gradient[x][index[1][i]][p];
						}
		
						// transition & emission = (BBB)
						miniMerge( i, intermediate[i], h, context, gradient );
					}
				}
				return res;
			}
		}
	}
	
	//add the partial derivation with given weight
	private void miniMerge( int i, double weight, int h, int context, double[][][] gradient ) {
		for( int p = 0; p < indicesTransition[i].length(); p++ ) {
			gradient[h][context][indicesTransition[i].get(p)] += weight * partDerTransition[i].get(p); 
		}
		int j = getIndex(index[0][i]);
		for( int p = 0; p < indicesState[j].length(); p++ ) {
			gradient[h][context][indicesState[j].get(p)] += weight * partDerState[j].get(p); 
		}
	}
	
	private void resetGradient( double[][][] gradient, int layer, int context, double val ) {
		Arrays.fill( gradient[layer % 2][context], val );
	}
	
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		int off = 0;
		for(int i=0;i<emission.length;i++){
			int num = ((DifferentiableEmission)emission[i]).getNumberOfParameters();
			if( num > 0){
				if(index >= off && index < off + num){
					return ((DifferentiableEmission)emission[i]).getSizeOfEventSpace();
				}
			}
			off += num;
		}
		return ((DifferentiableTransition)transition).getSizeOfEventSpace(index);
	}
	
	@Override
	public int[][] getSamplingGroups( int parameterOffset ) {
		LinkedList<int[]> list = new LinkedList<int[]>();
		for(int i=0;i<emission.length;i++){
			((DifferentiableEmission)emission[i]).fillSamplingGroups(parameterOffset, list);
		}
		((DifferentiableTransition)transition).fillSamplingGroups(parameterOffset, list);
		return list.toArray( new int[0][0] );
	}
	
	public String getInstanceName() {
		return "differentiable HMM(" + transition.getMaximalMarkovOrder() + ", " + training + ")";
	}
	
	public void setTrain( boolean train ) {
		this.train=train;
		if( train ) {
			if( training==TrainingType.LIKELIHOOD ) {
				score = Type.LIKELIHOOD;
			} else {
				score = Type.VITERBI;
			}
		} else {
			score = Type.LIKELIHOOD;
		}
	}
	
	//TODO to be removed
	public void check( DataSet data ) throws Exception {
		double[] params = new double[1+getNumberOfParameters()];
		System.arraycopy( getCurrentParameterValues(), 0, params, 1, params.length-1);
		setTrain(true);
		double[][] w = new double[1][data.getNumberOfElements()];
		Arrays.fill(w[0], 1);
		LogGenDisMixFunction log = new LogGenDisMixFunction( ((NumericalHMMTrainingParameterSet) trainingParameter).getNumberOfThreads(), new DifferentiableHigherOrderHMM[] { this }, new DataSet[] { data }, w, null, LearningPrinciple.getBeta(ess==0?LearningPrinciple.ML:LearningPrinciple.MAP), true, false);
		log.reset();
		double[] grad = log.evaluateGradientOfFunction(params);
		int min = ToolBox.getMinIndex(grad);
		int max = ToolBox.getMaxIndex(grad);
		System.out.println( "grad\t" + grad.length + "\t" + grad[min] + "\t" + grad[max] + "\t"  + dist(grad,null) + "\t" + Arrays.toString(grad)	);
		setTrain(false);
	}
	
	private static double dist( double[] a, double[] b ) {
		double dist=0;
		for( int i = 0; i < a.length; i++ ) {
			double diff = a[i] - (b==null?0:b[i]);
			dist += diff*diff;
		}
		return Math.sqrt(dist)/a.length;
	}
}