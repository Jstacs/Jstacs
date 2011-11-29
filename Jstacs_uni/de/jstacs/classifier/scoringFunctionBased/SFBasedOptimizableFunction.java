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

package de.jstacs.classifier.scoringFunctionBased;

import java.util.Arrays;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.classifier.scoringFunctionBased.logPrior.DoesNothingLogPrior;
import de.jstacs.classifier.scoringFunctionBased.logPrior.LogPrior;
import de.jstacs.data.Sample;
import de.jstacs.scoringFunctions.NormalizableScoringFunction;
import de.jstacs.scoringFunctions.ScoringFunction;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This abstract class is the basis of all multi-threaded {@link OptimizableFunction}s that are based on {@link ScoringFunction}s. 
 * 
 * @author Jens Keilwagen
 */
public abstract class SFBasedOptimizableFunction extends AbstractMultiThreadedOptimizableFunction {

	/**
	 * These shortcuts indicate the beginning of a new part in the parameter vector.
	 * Has to be set by {@link #reset()}.
	 */
	protected int[] shortcut;
	
	/**
	 * These {@link ScoringFunction}s are used during the parallel computation.
	 * <code>score[t]</code> contains all {@link ScoringFunction}s that are used by thread <code>t</code>.
	 */
	protected ScoringFunction[][] score;

	/**
	 * These {@link DoubleList}s are used during the parallel computation of the gradient.
	 * <code>dList[t]</code> contains all {@link DoubleList}s that are used by thread <code>t</code>.
	 */
	protected DoubleList[][] dList;

	/**
	 * These {@link IntList}s are used during the parallel computation of the gradient.
	 * <code>iList[t]</code> contains all {@link IntList}s that are used by thread <code>t</code>.
	 */
	protected IntList[][] iList;

	/**
	 * The prior that is used in this function.
	 */
	protected LogPrior prior;
	
	/**
	 * Creates an instance with the underlying infrastructure.
	 * Before using this instance, one has to invoke {@link #reset()}.
	 * 
	 * @param threads the number of threads used for evaluating the function and determining the gradient of the function
	 * @param score an array containing the {@link ScoringFunction}s that are used for determining the sequences scores
	 * @param data the array of {@link Sample}s containing the data that is needed to evaluate the function
	 * @param weights the weights for each {@link de.jstacs.data.Sequence} in each {@link Sample} of  <code>data</code> 
	 * @param prior the prior that is used for learning the parameters 
	 * @param norm
	 *            the switch for using the normalization (division by the number
	 *            of sequences)
	 * @param freeParams
	 *            the switch for using only the free parameters
	 * 
	 * @throws IllegalArgumentException
	 *             if the number of threads is not positive, the number of classes or the dimension of the weights is not correct
	 */
	public SFBasedOptimizableFunction( int threads, ScoringFunction[] score, Sample[] data, double[][] weights, LogPrior prior, boolean norm, boolean freeParams )																														throws IllegalArgumentException {
		super( threads, data, weights, norm, freeParams );
		shortcut = new int[cl + 1];
		if( freeParams ) {
			shortcut[0] = cl - 1;
		} else {
			shortcut[0] = cl;
		}
		this.prior = (prior == null) ? DoesNothingLogPrior.defaultInstance : prior;
		dList = new DoubleList[threads][cl];
		iList = new IntList[threads][cl];
		this.score = new ScoringFunction[threads][cl];
		int i = 0, j;
		for( ; i < cl; i++ )
		{
			this.score[0][i] = score[i];
			for( j = 0; j < threads; j++ )
			{
				dList[j][i] = new DoubleList();
				iList[j][i] = new IntList();
			}
		}
	}
	
	

	/**
	 * Returns from the complete vector of parameters those that are for the
	 * classes.
	 * 
	 * @param params
	 *            the current parameters
	 * 
	 * @return the parameters for the classes
	 */
	public final double[] getClassParams( double[] params ) {
		double[] res = new double[cl];
		System.arraycopy( params, 0, res, 0, shortcut[0] );
		if( freeParams ) {
			res[shortcut[0]] = 0;
		}
		return res;
	}
	
	/* (non-Javadoc)
	 * @see de.jstacs.algorithms.optimization.Function#getDimensionOfScope()
	 */
	public final int getDimensionOfScope() {
		return shortcut[cl];
	}

	protected void setThreadIndependentParameters() throws DimensionException {
		if( params == null || params.length != shortcut[cl] ) {
			if( params != null ) {
				throw new DimensionException( params.length, shortcut[cl] );
			} else {
				throw new DimensionException( 0, shortcut[cl] );
			}
		}
		for( int counter1 = 0; counter1 < shortcut[0]; counter1++ ) {
			logClazz[counter1] = params[counter1];
			clazz[counter1] = Math.exp( logClazz[counter1] );
		}
		if( freeParams ) {
			clazz[cl-1] = Math.exp( logClazz[cl-1] );
		}
	}
	
	public void getParameters( KindOfParameter kind, double[] erg ) throws Exception
	{
		switch( kind ) {
			case PLUGIN:
				double[] ess = new double[score[0].length];
				double discount = 0, e = 0;
				for( int i = 0; i < cl; i++ )
				{
					if( score[0][i] instanceof NormalizableScoringFunction ) {
						ess[i] = ((NormalizableScoringFunction)score[0][i]).getEss();
					}
					e += ess[i];
				}
				if( freeParams )
				{
					discount = score[0][cl-1].getInitialClassParam( (sum[cl - 1] + ess[cl-1] ) / ( sum[cl] + e ) );
				}
				for( int i = 0; i < shortcut[0]; i++ )
				{
					erg[i] = score[0][i].getInitialClassParam( (sum[i] + ess[i] ) / ( sum[cl] + e ) ) - discount;
				}
				break;
			case LAST:
				for( int i = 0; i < shortcut[0]; i++ )
				{
					erg[i] = logClazz[i];
				}
				break;
			case ZEROS:
				Arrays.fill( erg, 0 );
				break;
			default:
				throw new IllegalArgumentException( "Unknown kind of parameter" );
		}
		for( int i = 0; i < cl; i++ )
		{
			System.arraycopy( score[0][i].getCurrentParameterValues(), 0, erg, shortcut[i], score[0][i]
					.getNumberOfParameters() );
		}
	}
	
	protected void setParams( int index ) throws DimensionException
	{
		for( int counter1 = 0; counter1 < cl; counter1++ )
		{
			score[index][counter1].setParameters( params, shortcut[counter1] );
		}
	}

	/**
	 * This method adds the <code>term</code> to the class parameter of the
	 * class with index <code>classIndex</code>.
	 * 
	 * @param classIndex
	 *            the index of the class
	 * @param term
	 *            the term to be added to the class parameter
	 */
	public final void addTermToClassParameter( int classIndex, double term ) {
		if( classIndex < 0 || classIndex >= cl ) {
			throw new IndexOutOfBoundsException( "check the class index" );
		}
		if( freeParams && classIndex == cl - 1 ) {
			// this parameter is not free so we have to change all other class parameters
			for( int i = 0; i < shortcut[0]; i++ ) {
				logClazz[i] -= term;
				clazz[i] = Math.exp( logClazz[i] );
			}
		} else {
			logClazz[classIndex] += term;
			clazz[classIndex] = Math.exp( logClazz[classIndex] );
		}
	}
	
	/**
	 * This method allows to reset the internally used functions and the corresponding objects.
	 * 
	 * @param funs the new instances
	 * 
	 * @throws Exception if something went wrong
	 */
	public abstract void reset( ScoringFunction[] funs ) throws Exception;

	public final void reset() throws Exception {
		reset( score[0] );
		System.gc();
	}
}
