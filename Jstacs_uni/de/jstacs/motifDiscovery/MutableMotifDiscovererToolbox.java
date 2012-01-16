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

package de.jstacs.motifDiscovery;

import java.io.OutputStream;
import java.util.Arrays;

import de.jstacs.algorithms.optimization.ConstantStartDistance;
import de.jstacs.algorithms.optimization.DifferentiableFunction;
import de.jstacs.algorithms.optimization.NegativeDifferentiableFunction;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.StartDistanceForecaster;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.algorithms.optimization.termination.CombinedCondition;
import de.jstacs.algorithms.optimization.termination.IterationCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.classifier.differentiableSequenceScoreBased.DiffSSBasedOptimizableFunction;
import de.jstacs.classifier.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.data.DataSet;
import de.jstacs.data.RecyclableSequenceEnumerator;
import de.jstacs.data.Sequence;
import de.jstacs.data.WrongLengthException;
import de.jstacs.io.ArrayHandler;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder.RandomSeqType;
import de.jstacs.motifDiscovery.history.History;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.SafeOutputStream;

/**
 * This class contains some important methods for the initiation and optimization of {@link MutableMotifDiscoverer}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public final class MutableMotifDiscovererToolbox extends MotifDiscovererToolBox {
	
	/**
	 * This method allows to enumerate all possible seeds for a motif in the {@link MutableMotifDiscoverer} of a specific class.
	 * 
	 * @param funs the {@link DifferentiableSequenceScore}s
	 * @param classIndex the index of the class 
	 * @param motifIndex the index of the motif in the {@link MutableMotifDiscoverer}
	 * @param rse an {@link RecyclableSequenceEnumerator} that contains {@link Sequence} objects tested for initialization of the motif <code>motifIndex</code>
	 * @param weight the weight of the seed {@link Sequence}
	 * @param opt the objective function
	 * @param out a stream that allows to write some output if necessary
	 * 
	 * @return the best {@link Sequence} with respect to the {@link DiffSSBasedOptimizableFunction} 
	 * 
	 * @throws Exception if something went wrong 
	 * 
	 * @see MutableMotifDiscovererToolbox#enumerate(DifferentiableSequenceScore[], int[], int[], RecyclableSequenceEnumerator[], double, DiffSSBasedOptimizableFunction, OutputStream)
	 */
	public static Sequence enumerate( DifferentiableSequenceScore[] funs, int classIndex, int motifIndex, RecyclableSequenceEnumerator rse, double weight, DiffSSBasedOptimizableFunction opt, OutputStream out ) throws Exception
	{
		return enumerate( funs, new int[]{classIndex}, new int[]{motifIndex}, new RecyclableSequenceEnumerator[]{rse}, weight, opt, out )[0];
	}
	
	/**
	 * This method allows to enumerate all possible seeds for a number of motifs in the {@link MutableMotifDiscoverer}s of a specific classes.
	 * 
	 * @param funs the {@link DifferentiableSequenceScore}s
	 * @param classIndex the indices of the classes 
	 * @param motifIndex the indices of the motif in the {@link MutableMotifDiscoverer}s
	 * @param rse an array of {@link RecyclableSequenceEnumerator} that contains {@link Sequence} objects tested for initialization of the corresponding motif
	 * @param weight the weight of the seed {@link Sequence}
	 * @param opt the objective function
	 * @param out a stream that allows to write some output if necessary
	 * 
	 * @return an array containing the best {@link Sequence}s with respect to the {@link DiffSSBasedOptimizableFunction} 
	 * 
	 * @throws Exception if something went wrong 
	 */
	public static Sequence[] enumerate( DifferentiableSequenceScore[] funs, int[] classIndex, int[] motifIndex, RecyclableSequenceEnumerator[] rse, double weight, DiffSSBasedOptimizableFunction opt, OutputStream out ) throws Exception
	{
		DataSet[] data = opt.getData();
		double[][] dataWeights = opt.getSequenceWeights();
			
		int num = 0, idx, i;
		Sequence[] seq = new Sequence[classIndex.length], bestSeq = new Sequence[classIndex.length];
		DataSet[] s = new DataSet[classIndex.length];
		
		boolean[] adjust = new boolean[funs.length];
		Arrays.fill( adjust, false );
		for( i = 0; i < classIndex.length; i++ ) {
			adjust[classIndex[i]] = true;
		}

		MutableMotifDiscoverer[] mmd = new MutableMotifDiscoverer[funs.length];
		for( i = 0; i < funs.length; i++ ) {
			if( adjust[i] ) {
				mmd[i] = (MutableMotifDiscoverer) funs[i];
			} else {
				mmd[i] = null;
			}
		}
		
		int[] len = new int[classIndex.length];
		for( i = 0; i < classIndex.length; i++ ) {
			len[i] = mmd[classIndex[i]].getMotifLength( motifIndex[i] );
			rse[i].reset();
			seq[i] = rse[i].nextElement();
			s[i] = new DataSet( "sample " + i, seq[i] );
		}
		idx = classIndex.length-1;
		
		double curr, best = Double.NEGATIVE_INFINITY;
		double[] pars;
		double[][] weights = new double[motifIndex.length][];
		Arrays.fill(  weights, new double[]{weight} );
		while( true )
		{
			initMotif( idx, classIndex, motifIndex, s, weights, adjust, mmd, len, data, dataWeights );
			
			opt.reset();
			pars = opt.getParameters( KindOfParameter.PLUGIN );
			
			curr = opt.evaluateFunction(pars);
			out.write( (num++ + "\t" + Arrays.toString( seq ) + "\t" + curr + "\n").getBytes() );

			if( curr > best ) {
				best = curr;
				System.arraycopy( seq, 0, bestSeq, 0, seq.length );
			}
			
			idx = 0;
			while( idx < rse.length && !rse[idx].hasMoreElements() ) {
				rse[idx].reset();
				seq[idx] = rse[idx].nextElement();
				s[idx] = new DataSet( "sample " + idx, seq[idx] );
				idx++;
			}
			if( idx < rse.length ){
				seq[idx] = rse[idx].nextElement();
				s[idx] = new DataSet( "sample " + idx, seq[idx] );
			} else {
				break;
			}
		}
		out.write( ( "best: " + Arrays.toString( bestSeq ) + " " + best + "\n" ).getBytes() );
		for( i = 0; i < classIndex.length; i++ ) {
			s[i] = new DataSet( "sample " + i, bestSeq[i] );
		}
		initMotif( classIndex.length-1, classIndex, motifIndex, s, weights, adjust, mmd, len, data, dataWeights );
		
		return bestSeq;
	}
	
	/**
	 * This method allows to initialize a number of motifs.
	 * 
	 * @param idx the index indicates how many motifs are initialized
	 * @param classIndex the indices of the classes of each motif
	 * @param motifIndex the indices of each motif within the {@link MutableMotifDiscoverer}
	 * @param s the {@link DataSet}s to be used for the initialization
	 * @param seqWeights the weights corresponding to the {@link DataSet}s
	 * @param adjust an array of switches indicating whether to adjust hidden parameters of not
	 * @param mmd the array of {@link MutableMotifDiscoverer}s to be initialized
	 * @param len the length of each motif
	 * @param data the complete data sets
	 * @param dataWeights the weights corresponding to the complete data sets
	 * 
	 * @throws Exception if something went wrong
	 */
	public static void initMotif( int idx, int[] classIndex, int[] motifIndex, DataSet[] s, double[][] seqWeights, boolean[] adjust, MutableMotifDiscoverer[] mmd, int[] len, DataSet data[], double[][] dataWeights ) throws Exception {
		int i, sl;
		for( i = 0; i <= idx; i++ ) {
			sl = s[i].getElementLength();
			if( sl > len[i] ){
				throw new WrongLengthException( sl );
			} else if( sl < len[i] ){
				mmd[classIndex[i]].modifyMotif( motifIndex[i], 0, sl-len[i] );
			}
			
			//System.out.println( i + "\t" + classIndex[i] + "\t" + motifIndex[i] + "\t" + s[i].getElementAt(0) );
						
			mmd[classIndex[i]].initializeMotif( motifIndex[i], s[i], seqWeights[i] );
			if( sl < len[i] ){
				mmd[classIndex[i]].modifyMotif( motifIndex[i], -(int)Math.floor( (len[i]-sl)/2.0 ), (int)Math.ceil( (len[i]-sl)/2.0 ) );
			}
		}
		for( i = 0; i < adjust.length; i++ ) {
			if( adjust[i] ) {
				mmd[classIndex[i]].adjustHiddenParameters( classIndex[i], data, dataWeights );
			}
		}
	}
	
	/**
	 * This enum defines some constants for the method {@link MutableMotifDiscovererToolbox#getSortedInitialParameters(DifferentiableSequenceScore[], InitMethodForScoringFunction[], DiffSSBasedOptimizableFunction, int, OutputStream, int)}.
	 * These constants define how to initialize the {@link DifferentiableSequenceScore}s. 
	 * 
	 * @author Jens Keilwagen
	 */
	public static enum InitMethodForScoringFunction {
		/**
		 * This constants indicates that a {@link DifferentiableSequenceScore} should be initialized using {@link DifferentiableSequenceScore#initializeFunction(int, boolean, DataSet[], double[][])}.
		 */
		PLUG_IN,
		/**
		 * This constants indicates that the motifs of a {@link MutableMotifDiscoverer} should be initialized using {@link MutableMotifDiscoverer#initializeMotifRandomly(int)}.
		 */
		MOTIF_RANDOMLY,
		/**
		 * This constants indicates that a {@link DifferentiableSequenceScore} should be initialized using {@link DifferentiableSequenceScore#initializeFunctionRandomly(boolean)}.
		 */
		RANDOMLY,
		/**
		 * This constants indicates that a {@link DifferentiableSequenceScore} should not be initialized, i.e. the instance is not changed and uses the current parameters.
		 */
		NOTHING;
	}
	
	/**
	 * This method allows to initialize the {@link DifferentiableSequenceScore} using different {@link InitMethodForScoringFunction}. It returns an array of {@link ComparableElement}s that contain the parameters and the
	 * 
	 * @param funs the {@link DifferentiableSequenceScore}s
	 * @param init the specific {@link InitMethodForScoringFunction}, the entries correspond one to one to those of <code>fun</code>
	 * @param opt the objective function
	 * @param n the number of initializations
	 * @param stream a stream that allows to write some output if necessary
	 * @param optimizationSteps the number of initial steps that should be performed before evaluating the function
	 * 
	 * @return a sorted array containing {@link ComparableElement}s of parameter arrays and corresponding values of the {@link DiffSSBasedOptimizableFunction}
	 * 
	 * @throws Exception if something went wrong
	 */
	@SuppressWarnings("unchecked")
	public static ComparableElement<double[],Double>[] getSortedInitialParameters( DifferentiableSequenceScore[] funs, InitMethodForScoringFunction[] init, DiffSSBasedOptimizableFunction opt, int n, OutputStream stream, int optimizationSteps ) throws Exception
	{
		SafeOutputStream info = SafeOutputStream.getSafeOutputStream(stream);
		DataSet[] data = opt.getData();
		double[][] oldParams = new double[funs.length][];
		for( int j = 0; j < funs.length; j++ ) {
			if( init[j] == InitMethodForScoringFunction.NOTHING ) {
				oldParams[j] = funs[j].getCurrentParameterValues();
			}			
		}
		
		ComparableElement<double[],Double>[] erg = new ComparableElement[n];
		double[] params;
		double c;
		ConstantStartDistance cs = new ConstantStartDistance( 1 );
		SafeOutputStream out = SafeOutputStream.getSafeOutputStream( null );
		NegativeDifferentiableFunction nOpt = new NegativeDifferentiableFunction( opt );
		TerminationCondition condition = new IterationCondition( optimizationSteps );
		for( int j, i = 0; i < n; i++ )
		{
			for( j = 0; j < funs.length; j++ )
			{
				switch( init[j] )
				{
					case PLUG_IN: 
						funs[j].initializeFunction( j, false, data, null );
						break;
					case MOTIF_RANDOMLY:
						if( funs[j] instanceof MutableMotifDiscoverer ) {
							MutableMotifDiscoverer mmd = (MutableMotifDiscoverer) funs[j];
							int m = mmd.getNumberOfMotifs();
							for( int k = 0; k < m; k++ ) {
								mmd.initializeMotifRandomly( k );
							}
						}
						break;
					case RANDOMLY: 
						funs[j].initializeFunctionRandomly( false );
						break;
					case NOTHING:
						//leave the parameters as they are
						funs[j].setParameters( oldParams[j], 0 );
						break;
				}
				
			}
			opt.reset();
			params = opt.getParameters( KindOfParameter.PLUGIN );
			/*c = opt.evaluateFunction( params );
			stream.write( i + "\t" + c + "\t" ); //TODO*/
			if( optimizationSteps > 0 ) {
				Optimizer.optimize( (byte) 10, nOpt, params, condition, 1E-10, cs, out );
			}
			c = opt.evaluateFunction( params );
			info.writeln( i + "\t" + c );
			erg[i] = new ComparableElement<double[],Double>( params, c );			
		}
		Arrays.sort( erg );
		info.writeln( "[" + erg[0].getWeight() + " .. " + erg[n-1].getWeight() + "]" );
		return erg;
	}
	
	/**
	 * This method creates a minimalNewLength-array that can be used in an optimization.
	 * 
	 * @param funs the ScoringFunctions used in an optimization
	 * 
	 * @return an minimalNewLength-array for the given ScoringFunctions
	 */
	public static int[][] createMinimalNewLengthArray( DifferentiableSequenceScore[] funs ) {
		int[][] minimalNewLength = new int[funs.length][];
		MutableMotifDiscoverer disc;
		for( int i, j = 0; j < funs.length; j++ )
		{
			if( funs[j] instanceof MutableMotifDiscoverer )
			{
				disc = (MutableMotifDiscoverer)funs[j];
				minimalNewLength[j] = new int[disc.getNumberOfMotifs()];
				for( i = 0; i < minimalNewLength[j].length; i++ )
				{
					minimalNewLength[j][i] = disc.getMotifLength( i );
				}
			}
		}
		return minimalNewLength;
	}
	
	/**
	 * This method creates a History-array that can be used in an optimization.
	 * 
	 * @param funs the ScoringFunctions used in an optimization
	 * @param template the template history instance
	 * 
	 * @return an History-array for the given ScoringFunctions
	 * 
	 * @throws CloneNotSupportedException if the t<code>template</code> could not be cloned
	 */
	public static History[][] createHistoryArray( DifferentiableSequenceScore[] funs, History template ) throws CloneNotSupportedException{
		History[][] history = new History[funs.length][];
		for( int i, j = 0; j < funs.length; j++ )
		{
			if( funs[j] instanceof MutableMotifDiscoverer )
			{
				history[j] = new History[((MutableMotifDiscoverer)funs[j]).getNumberOfMotifs()];
				if( template != null ) {
					for( i = 0; i < history[j].length; i++ )
					{
						history[j][i] = template.clone();
					}
				}
			}
		}
		return history;
	}
	
	/**
	 * This method clears all elements of an History-array, so that it can be used again.
	 * 
	 * @param history the array
	 */
	public static void clearHistoryArray( History[][] history ) {
		for( int i, j = 0; j < history.length; j++ )
		{
			if( history[j] != null )
			{
				for( i = 0; i < history[j].length; i++ )
				{
					if( history[j][i] != null )
					{
						history[j][i].clear();
					}
				}
			}
		}
	}
	
	/**
	 * This method tries to optimize the problem at hand as good as possible. If the optimization uses {@link MutableMotifDiscoverer}s it tries to perform modify operations as long as they seem to be promising.
	 * 
	 * @param funs the {@link DifferentiableSequenceScore}s for scoring sequences
	 * @param opt the {@link DiffSSBasedOptimizableFunction}
	 * @param algorithm used for the optimization
	 * @param condition used for the optimization
	 * @param linEps used for the optimization
	 * @param startDistance used for the optimization
	 * @param out an stream that allows to obtain some information while optimization
	 * @param breakOnChanged a switch that decides whether a new optimization should be started after one successful modify or after all motifs have been tried to modify.
	 * @param template a history instance used to build an array with this instance
	 * @param plugIn a switch whether to take the internal parameters or not
	 * @param maxPos a switch whether to take the maximal shift position or not in the heuristic
	 * 
	 * @return the optimized value (res[0][0]) and the array for the class parameters (res[1])
	 * 
	 * @throws Exception if something went wrong while optimization
	 *
	 * @see MutableMotifDiscovererToolbox#clearHistoryArray(de.jstacs.motifDiscovery.history.History[][])
	 * @see MutableMotifDiscovererToolbox#optimize(DifferentiableSequenceScore[], DiffSSBasedOptimizableFunction, byte, AbstractTerminationCondition, double, StartDistanceForecaster, SafeOutputStream, boolean, History[][], int[][], de.jstacs.classifier.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter, boolean)
	 */
	public static double[][] optimize( DifferentiableSequenceScore[] funs, DiffSSBasedOptimizableFunction opt, byte algorithm, AbstractTerminationCondition condition, double linEps, StartDistanceForecaster startDistance, SafeOutputStream out, boolean breakOnChanged, History template, KindOfParameter plugIn, boolean maxPos ) throws Exception {
		return optimize( funs, opt, algorithm, condition, linEps, startDistance, out, breakOnChanged, createHistoryArray( funs, template ), createMinimalNewLengthArray( funs ), plugIn, maxPos );
	}
	
	/**
	 * This method tries to optimize the problem at hand as good as possible. If the optimization uses {@link MutableMotifDiscoverer}s it tries to perform modify operations as long as they seem to be promising.
	 * 
	 * @param funs the {@link DifferentiableSequenceScore}s for scoring sequences 
	 * @param opt the {@link DiffSSBasedOptimizableFunction}
	 * @param algorithm used for the optimization
	 * @param condition used for the optimization
	 * @param linEps used for the optimization
	 * @param startDistance used for the optimization
	 * @param out an stream that allows to obtain some information while optimization
	 * @param breakOnChanged a switch that decides whether a new optimization should be started after one successful modify or after all motifs have been tried to modify.
	 * @param hist an array that is used to check whether a modify-operation can be performed
	 * @param minimalNewLength the minimal new length for each motif in each class, that will be used in an expand if the motif was shortened before
	 * @param plugIn a switch whether to take the internal parameters or not
	 * @param maxPos a switch whether to take the maximal shift position or not in the heuristic
	 * 
	 * @return the optimized value (res[0][0]) and the array for the class parameters (res[1])
	 * 
	 * @throws Exception if something went wrong while optimization
	 */
	public static double[][] optimize( DifferentiableSequenceScore[] funs, DiffSSBasedOptimizableFunction opt, byte algorithm, AbstractTerminationCondition condition, double linEps, StartDistanceForecaster startDistance, SafeOutputStream out, boolean breakOnChanged, History[][] hist, int[][] minimalNewLength, KindOfParameter plugIn, boolean maxPos ) throws Exception {
		NegativeDifferentiableFunction neg = new NegativeDifferentiableFunction(opt);
		int k;
		DifferentiableSequenceScore[] best = null;
		double[] params, classParams = null;
		DataSet[] data = opt.getData();
		double[][] weights = opt.getSequenceWeights();
		double bestVal = Double.NEGATIVE_INFINITY, current;
		do{
			opt.reset();
			startDistance.reset();
			params = opt.getParameters( plugIn );
			plugIn = KindOfParameter.LAST;
			Optimizer.optimize( algorithm, neg, params, condition, linEps, startDistance, out );
			current = opt.evaluateFunction( params );
			if( current > bestVal )
			{
				best = null;
				System.gc();
				best = ArrayHandler.clone( funs );
				bestVal = current;
				classParams = opt.getClassParams(params);
			}
 			if( condition instanceof SmallDifferenceOfFunctionEvaluationsCondition ) {
 				condition = new CombinedCondition( 1, condition, steps );
 			}
		}while( doHeuristicSteps( funs, data, weights, opt, neg, algorithm, linEps, startDistance, out, breakOnChanged, hist, minimalNewLength, maxPos ) );
		for( k = 0; k < funs.length; k++ )
		{
			funs[k] = best[k];
		}
		return new double[][]{{bestVal}, classParams};
	}
		
	/**
	 * This method tries to make some heuristic step if at least one {@link DifferentiableSequenceScore} is a {@link MutableMotifDiscoverer}.
	 * These heuristic steps include shift, shrink, and expand as far as the user allows those operations by the {@link History} array.
	 * 
	 * @param funs the {@link DifferentiableSequenceScore}s  for scoring sequences
	 * @param data array of {@link DataSet} containing the data for each class
	 * @param weights the weights corresponding to the {@link Sequence}s in <code>data</code>
	 * @param opt the {@link DiffSSBasedOptimizableFunction}
	 * @param neg the {@link NegativeDifferentiableFunction} used in the optimization
	 * @param algorithm used for the optimization
	 * @param linEps used for the optimization
	 * @param startDistance used for the optimization
	 * @param out an stream that allows to obtain some information while optimization
	 * @param breakOnChanged a switch that decides whether a new optimization should be started after one successful modify or after all motifs have been tried to modify.
	 * @param hist an array that is used to check whether a modify-operation can be performed
	 * @param minimalNewLength the minimal new length for each motif in each class, that will be used in an expand if the motif was shortened before
	 * @param maxPos a switch whether to take the maximal shift position or not in the heuristic
	 * 
	 * @return <code>true</code> if some heuristic steps has been performed otherwise <code>false</code>
	 * 
	 * @throws Exception if something went wrong
	 */
	public static boolean doHeuristicSteps( DifferentiableSequenceScore[] funs, DataSet[] data, double[][] weights, DiffSSBasedOptimizableFunction opt,
			DifferentiableFunction neg, byte algorithm, double linEps, StartDistanceForecaster startDistance, 
			SafeOutputStream out, boolean breakOnChanged, History[][] hist, int[][] minimalNewLength, boolean maxPos ) throws Exception {
		boolean changed = false, changedThisOne;
		double normOld, normNew;
		DifferentiableStatisticalModel nsf;
		for( int k = 0; k < funs.length && !( changed && breakOnChanged ); k++ )
		{
			if( funs[k] instanceof MutableMotifDiscoverer  && ((MutableMotifDiscoverer)funs[k]).getNumberOfMotifs() > 0 )
			{
				out.writeln( "MutableMotifDiscoverer " + k + ":\n" + funs[k].toString() );
				if( funs[k] instanceof DifferentiableStatisticalModel ) {
					nsf = (DifferentiableStatisticalModel) funs[k];
					normOld = nsf.getLogNormalizationConstant();
				} else {
					nsf = null;
					normOld = 0;
				}
				MutableMotifDiscoverer currMD = (MutableMotifDiscoverer) funs[k];
				int numMotifs = currMD.getNumberOfMotifs();
				changedThisOne = false;
				for( int l = 0 ; l < numMotifs  && !( changedThisOne && breakOnChanged ); l++ ){
					if( hist [k][l] != null ) {
						changedThisOne |= 
							//modify( currMD, k, l, data, weights, hist[k][l], minimalNewLength[k][l], out );
							findModification( k, l, currMD, funs, data, weights, opt, neg, algorithm, linEps, startDistance, out, hist[k][l], minimalNewLength[k][l], maxPos );
					}
				}
				if( changedThisOne ) {
					changed = true;
					if( nsf != null ) {
						normNew = nsf.getLogNormalizationConstant();
						//set new class parameter
						opt.addTermToClassParameter( k, normOld - normNew );
					}
				}
			}
		}
		return changed;
	}
	
	private static double pVal = 1E-4, threshold = 0.8;
	private static AbstractTerminationCondition steps;
	static{
		try {
			steps = new IterationCondition( 10 );
		} catch ( Exception e ) {
			//will not happen
			throw new RuntimeException( e.getMessage() );
		}
	}
	
	/**
	 * This method tries to find a modification, i.e. shifting, shrinking, or expanding a motif, that is promising.
	 * The method returns <code>true</code> a modification was found and could be performed.
	 * 
	 * <br><br>
	 * 
	 * For finding a promising modification, the method test various shifts and computes the number of sequences predicted to be bound.
	 * 
	 * @param clazz the class index for which the Scoring function will be tested for modification
	 * @param motif the motif index for which the Scoring function will be tested for modification
	 * @param mmd the {@link MutableMotifDiscoverer} that will be tested
	 * @param score the {@link DifferentiableSequenceScore}s for scoring sequences
	 * @param data array of {@link DataSet} containing the data for each class
	 * @param weights array of <code>double[]</code> containing the weights for the data of each class
	 * @param opt the {@link DiffSSBasedOptimizableFunction}
	 * @param neg the {@link NegativeDifferentiableFunction} used in the optimization
	 * @param algo used for the optimization
	 * @param linEps used for the optimization
	 * @param startDistance used for the optimization
	 * @param out an stream that allows to obtain some information
	 * @param hist an instance to check whether a modify-operation can be performed
	 * @param minimalNewLength the minimal new length for each motif in each class, that will be used in an expand if the motif was shortened before
	 * @param maxPos a switch whether to take the maximal shift position or not in the heuristic
	 * 
	 * @return <code>true</code> if a modification has been performed
	 * 
	 * @throws Exception if something went wrong
	 * 
	 * @see SignificantMotifOccurrencesFinder#getNumberOfBoundSequences(DataSet, double[], int)
	 */
	public static boolean findModification( int clazz, int motif, MutableMotifDiscoverer mmd, DifferentiableSequenceScore[] score, DataSet[] data, double[][] weights, DiffSSBasedOptimizableFunction opt, DifferentiableFunction neg, byte algo, double linEps, StartDistanceForecaster startDistance, SafeOutputStream out, History hist, int minimalNewLength, boolean maxPos ) throws Exception {
		double[] params = opt.getParameters( KindOfParameter.LAST );
		int len = mmd.getMotifLength( motif ) / 2;
		DataSet[] my = {data[clazz], null};
		double[][] myWeights = { weights[clazz], null };
		if( data.length > 2 ) {
			boolean[] in = new boolean[data.length];
			Arrays.fill( in, true ) ;
			in[clazz] = false;
			my[1] = DataSet.union( data, in );
			int anz = 0;
			for( int i = 0; i < data.length; i++ ) {
				if( in[i] ) {
					anz += data[i].getNumberOfElements();
				}
			}
			myWeights[1] = new double[anz];
			anz = 0;
			for( int n, i = 0; i < data.length; i++ ) {
				if( in[i] ) {
				n = data[i].getNumberOfElements();
					if( weights[i] != null ) {
						System.arraycopy( weights[i], 0, myWeights[1], anz, n );
					} else {
						Arrays.fill( myWeights, anz, anz+n, 1d );
					}
					anz += n;
				}
			}
		} else {
			my[1] = data[1-clazz];
			myWeights[1] = weights[1-clazz];
		}
		SignificantMotifOccurrencesFinder smof;
		if( my[1] != null ) {
			smof = new SignificantMotifOccurrencesFinder( mmd, my[1], myWeights[1], pVal );
		} else {
			smof = new SignificantMotifOccurrencesFinder( mmd, RandomSeqType.PERMUTED, true, NUMBER_OF_PERMUTATIONS, pVal );
		}
		double current = smof.getNumberOfBoundSequences( my[0], myWeights[0], motif );
		out.writeln( "optimized predicted bound sequences: " + current );
		out.writeln( "====================================" );
	
		int[] notSignif = new int[2];
		
		out.writeln( "shift downstream" );
		notSignif[0] = heuristic( clazz, motif, len, 1, my, myWeights, opt, neg, algo, linEps, startDistance, params, score, current, out, maxPos );
		
		out.writeln();
		out.writeln( "shift upstream" );
		notSignif[1] = heuristic( clazz, motif, len, -1, my, myWeights, opt, neg, algo, linEps, startDistance, params, score, current, out, maxPos );
		
		out.writeln();
		out.writeln( "not significant: " + Arrays.toString( notSignif ) );
		boolean modified = modify( notSignif, mmd, clazz, motif, hist, minimalNewLength, out );
		if( modified ) {
			opt.reset();
		}
		return modified;
	}
	
	private static final SafeOutputStream DISCARD_OUT = SafeOutputStream.getSafeOutputStream( null );
	private static final int NUMBER_OF_PERMUTATIONS = 1000;
	
	private static int heuristic( int clazz, int motif, int len, int direction, DataSet[] data, double[][] weights, DiffSSBasedOptimizableFunction test, DifferentiableFunction neg, byte algo, double linEps, StartDistanceForecaster startDistance, double[] params, DifferentiableSequenceScore[] score, double pred, SafeOutputStream out, boolean maxPos ) throws Exception {
		test.setParams(params);
		MutableMotifDiscoverer mmd = (MutableMotifDiscoverer) score[clazz];
		//System.out.println( mmd );
		double[] params_copy;
		double normOld = Double.NEGATIVE_INFINITY;
		double[] val = new double[len];
		int dir1, dir2;
		dir1=dir2=direction;//TODO
		if( direction > 0 ) {
			dir1 = direction;
		} else {
			dir2 = direction;
		}
		int i;
		for( i=1;i<=len;i++){
			
			if( score[clazz] instanceof DifferentiableStatisticalModel ) {
				normOld = ((DifferentiableStatisticalModel) score[clazz]).getLogNormalizationConstant();
			}
			if( mmd.modifyMotif( motif, i*dir1, i*dir2 ) ){
				if( score[clazz] instanceof DifferentiableStatisticalModel ) {
					//set new class parameter
					test.addTermToClassParameter( 0, normOld - ((DifferentiableStatisticalModel) score[clazz]).getLogNormalizationConstant() );
				}
				
				test.reset();
				
				//optimize a little bit
				params_copy = test.getParameters( KindOfParameter.LAST );
				Optimizer.optimize( algo, neg, params_copy, steps, linEps, startDistance, DISCARD_OUT );
				test.setParams( params_copy );
				SignificantMotifOccurrencesFinder smof;
				if( data.length > 1 && data[1] != null ) {
					smof = new SignificantMotifOccurrencesFinder( mmd, data[1], weights[1], pVal );
				} else {
					smof = new SignificantMotifOccurrencesFinder( mmd, RandomSeqType.PERMUTED, true, NUMBER_OF_PERMUTATIONS, pVal );
				}
				val[i-1] = smof.getNumberOfBoundSequences( data[0], weights[0], motif );
				
				mmd.modifyMotif( motif, -i*dir1, -i*dir2 );
				test.reset();
				test.setParams( params );
			} else {
				val[i-1] = 0;
			}
			out.writeln( i + "\t" + val[i-1] );
		}
		if( !maxPos ) {
			i = 0;
			while( i < val.length && val[i] >= pred*threshold ) {
				i++;
			}
			i--;
		} else {
			i=-1;
			double max = pred*threshold;
			for(int j=0;j<val.length;j++){
				if(val[j] >= max){
					max = val[j];
					i=j;
				}
			}			
		}
		return i+1;
	}
	
	/*
	private static boolean modifyOld( MutableMotifDiscoverer md, int clazz, int motif, Sample[] data, double[][] weights, History hist, int minimalNewLength, SafeOutputStream out ) throws Exception
	{
		int[] notSignif = md.determineNotSignificantPositionsFor( motif, data, weights, clazz );
		return modify( notSignif, md, clazz, motif, hist, minimalNewLength, out );
	}
	*/
	
	private static boolean modify( int[] notSignif, MutableMotifDiscoverer md, int clazz, int motif, History hist, int minimalNewLength, SafeOutputStream out ) throws Exception
	{
		boolean modified = false;
		
		int sum = notSignif[0]+notSignif[1], ml = md.getMotifLength( motif );
		int[] modification = new int[2];
	
		if( sum > 0 ) { // at least one side is not significant
			if( sum < ml ) {
				//try to shift from the side with the most not significant positions
				if( notSignif[0] >= notSignif[1] ) {
					// (more) not significant positions at the left side
					modification[0] = notSignif[0];
				} else {
					// (more) not significant positions at the right side
					modification[0] = -notSignif[1];
				}
				modification[1] = modification[0];
				if( hist.operationAllowed( modification ) ){
					modified = md.modifyMotif( motif, modification[0], modification[1] );
				}
			}
			
			if( !modified ) {
				//shrink
				if( sum < ml ) {
					// shrink is possible
					modification[0] = notSignif[0];
					modification[1] = -notSignif[1];
				} else {
					// shrink to size 1
					modification[0] = 0;
					modification[1] = 1-md.getMotifLength( motif );
				}
				if( modification[0] != modification[1] // modification [0,0] forbidden
				        && hist.operationAllowed( modification ) ){
					modified = md.modifyMotif( motif, modification[0], modification[1] );
				}
			}
		} else {
			//expand
			if( md.getMotifLength( motif ) < minimalNewLength ) {
				// expand to minimal length on one side
				modification[0] = 0;
				modification[1] = minimalNewLength - md.getMotifLength( motif );
				if( !hist.operationAllowed( modification ) ){
					modification[0] = -(minimalNewLength - md.getMotifLength( motif ));
					modification[1] = 0;
				}
			} else {
				// expand 1 to the left and right
				modification[0] = -1;
				modification[1] = 1;
			}
			
			if( hist.operationAllowed( modification ) ){
				modified = md.modifyMotif( motif, modification[0], modification[1] );
			}
		}
		
		if( modified ){
			hist.operationPerfomed( modification );
			out.writeln( "class " + clazz + ": modified motif " + motif + " => " + Arrays.toString( modification ) );
		}
		return modified;
	}
}
