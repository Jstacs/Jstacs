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

package de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix;

import java.util.Arrays;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.classifiers.differentiableSequenceScoreBased.DiffSSBasedOptimizableFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.utils.Normalisation;

/**
 * This class implements the the following function
 * 
 * {@latex.ilb \\[f(\\underline{\\lambda}|C,D,\\underline{\\alpha},\\underline{\\beta})
 * := \\beta_0 \\log p(C|D,\\underline{\\lambda}) + \\beta_1 \\log p(C,D|\\underline{\\lambda}) + \\beta_2 \\log p(\\underline{\\lambda}|\\underline{\\alpha})
 * .\\]}
 * 
 * The weights {@latex.inline $\\beta_i$} have to sum to 1. For special weights the optimization turns out to be
 * well known
 * <ul>
 * <li> if the weights are (0,1,0), one obtains maximum likelihood,
 * <li> if the weights are (0,0.5,0.5), one obtains maximum a posteriori,
 * <li> if the weights are (1,0,0), one obtains maximum conditional likelihood,
 * <li> if the weights are (0.5,0,0.5), one obtains maximum supervised posterior,
 * <li> if the {@latex.inline $\\beta_2$}=0, one obtains the generative-discriminative trade-off,
 * <li> if the {@latex.inline $\\beta_2$}=0.5, one obtains the penalized generative-discriminative trade-off.
 * </ul>
 * 
 * Of course, there are also some very interesting cases with other weights.
 * 
 * <br>
 * <br>
 * 
 * It can be used to maximize the parameters.
 * 
 * <br>
 * <br>
 * 
 * This class enables the user to exploit all CPUs of the computer by using threads. The number of compute threads can be
 * determined in the constructor. 
 * 
 * <br>
 * <br>
 * 
 * It is very important for this class that the {@link DifferentiableSequenceScore#clone()} method works correctly, since each thread works on its own clones. 
 * 
 * @author Jens Keilwagen
 */
public class LogGenDisMixFunction extends DiffSSBasedOptimizableFunction
{
	/**
	 * General temporary array
	 */
	protected double[][] helpArray;
	
	/**
	 * Array for the gradient of the log-likelihood
	 */
	protected double[][] llGrad;
	/**
	 * Array for the gradient of the conditional log-likelihood
	 */
	protected double[][] cllGrad;

	/**
	 * The mixture parameters of the GenDisMix
	 */
	protected double[] beta;
	
	/**
	 * Array for the gradient of the prior
	 */
	protected double[] prGrad;

	/**
	 * The constructor for creating an instance that can be used in an {@link de.jstacs.algorithms.optimization.Optimizer}.
	 * 
	 * @param threads the number of threads used for evaluating the function and determining the gradient of the function
	 * @param score an array containing the {@link DifferentiableSequenceScore}s that are used for determining the sequences scores;
	 * 			if the weight <code>beta[LearningPrinciple.LIKELIHOOD_INDEX]</code> is positive all elements of <code>score</code> have to be {@link DifferentiableStatisticalModel}
	 * @param data the array of {@link DataSet}s containing the data that is needed to evaluate the function
	 * @param weights the weights for each {@link Sequence} in each {@link DataSet} of  <code>data</code> 
	 * @param prior the prior that is used for learning the parameters
	 * @param beta the beta-weights for the three terms of the learning principle 
	 * @param norm
	 *            the switch for using the normalization (division by the number
	 *            of sequences)
	 * @param freeParams
	 *            the switch for using only the free parameters
	 * 
	 * @throws IllegalArgumentException
	 *             if the number of threads is not positive, the number of classes or the dimension of the weights is not correct
	 */
	public LogGenDisMixFunction( int threads, DifferentiableSequenceScore[] score, DataSet[] data, double[][] weights,
			LogPrior prior, double[] beta, boolean norm, boolean freeParams ) throws IllegalArgumentException
	{
		super( threads, score, data, weights, prior, norm, freeParams );
		if( cl < 1 || ( beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] != 0 && cl < 2 ) ) 
		{
			throw new IllegalArgumentException(
					"The number of classes is not correct. You can use this class for the generative training of one class or the (in some kind) discriminative training of more than one class." );
		}
		this.beta = LearningPrinciple.checkWeights( beta );
		check();
		helpArray = new double[threads][Math.max(2,cl)];
	}

	protected double[] joinGradients() throws EvaluationException
	{	
		int index = 0, t;
		for( ; index < llGrad[0].length; index++ )
		{
			for( t = 1; t < llGrad.length; t++ )
			{
				llGrad[0][index] += llGrad[t][index];
				cllGrad[0][index] += cllGrad[t][index];
			}
		}

		if( beta[LearningPrinciple.LIKELIHOOD_INDEX] != 0 )
		{
			try
			{
				int counter1, counter3;
				double norm = Double.NEGATIVE_INFINITY;
				for( counter3 = 0; counter3 < cl; counter3++ )
				{
					norm = Normalisation.getLogSum( norm, logClazz[counter3] + ((DifferentiableStatisticalModel)score[0][counter3]).getLogNormalizationConstant() );
				}
				for( counter3 = 0; counter3 < cl; counter3++ )
				{
					if( counter3 < shortcut[0] )
					{
						llGrad[0][counter3] -= sum[cl] * Math.exp( logClazz[counter3] + ((DifferentiableStatisticalModel)score[0][counter3]).getLogNormalizationConstant() - norm );
					}
					for( counter1 = shortcut[counter3]; counter1 < shortcut[counter3 + 1]; counter1++ )
					{
						llGrad[0][counter1] -= sum[cl] * Math.exp( logClazz[counter3] + ((DifferentiableStatisticalModel)score[0][counter3]).getLogPartialNormalizationConstant( counter1 - shortcut[counter3] ) - norm );
					}
				}
			}
			catch( Exception e )
			{
				EvaluationException eva = new EvaluationException( e.getMessage() );
				eva.setStackTrace( e.getStackTrace() );
				throw eva;
			}
		}
		
		// prior
		Arrays.fill( prGrad, 0 );
		if( beta[LearningPrinciple.PRIOR_INDEX] != 0 )
		{
			prior.addGradientFor( params, prGrad );
		}
		// to avoid: 0*Infinity, ...
		if( beta[LearningPrinciple.LIKELIHOOD_INDEX] == 0 )
		{
			Arrays.fill( llGrad[0], 0 );
		}
		if( beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] == 0 )
		{
			Arrays.fill( cllGrad[0], 0 );
		}

		double[] grad = new double[shortcut[cl]];
		double weight;
		// normalization
		if( norm )
		{
			weight = sum[cl];
		}
		else
		{
			weight = 1d;
		}
		
		for( index = 0; index < grad.length; index++ )
		{
			grad[index] = (beta[LearningPrinciple.LIKELIHOOD_INDEX] * llGrad[0][index]
			                   + beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] * cllGrad[0][index]
			                   + beta[LearningPrinciple.PRIOR_INDEX] * prGrad[index])
			                   / weight;
		}
		return grad;
	}
	
	protected void evaluateGradientOfFunction( int index, int startClass, int startSeq, int endClass, int endSeq )
	{
		Arrays.fill( llGrad[index], 0 );
		Arrays.fill( cllGrad[index], 0 );

		double weight;
		int counter1, counter2, counter3 = startClass, counter4 = 0, start, end;

		Sequence s;
		for( ; counter3 <= endClass; counter3++ )
		{
			if( counter3 == startClass )
			{
				start = startSeq;
			}
			else
			{
				start = 0;
			}
			if( counter3 == endClass )
			{
				end = endSeq;
			}
			else
			{
				end = data[counter3].getNumberOfElements();
			}
			
			for( counter2 = start; counter2 < end; counter2++ )
			{
				s = data[counter3].getElementAt( counter2 );
				weight = weights[counter3][counter2];
				if( beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] != 0 )
				{
					for( counter1 = 0; counter1 < cl; counter1++ )
					{
						iList[index][counter1].clear();
						dList[index][counter1].clear();
						helpArray[index][counter1] = logClazz[counter1]
								+ score[index][counter1].getLogScoreAndPartialDerivation( s, 0, iList[index][counter1],
										dList[index][counter1] );
					}
				}
				else
				{
					iList[index][counter3].clear();
					dList[index][counter3].clear();
					helpArray[index][counter3] = logClazz[counter3]
							+ score[index][counter3].getLogScoreAndPartialDerivation( s, 0, iList[index][counter3], dList[index][counter3] );
				}
				if( beta[LearningPrinciple.LIKELIHOOD_INDEX] != 0 )
				{
					if( counter3 < shortcut[0] )
					{
						llGrad[index][counter3] += weight;
					}
					for( counter4 = 0; counter4 < iList[index][counter3].length(); counter4++ )
					{
						llGrad[index][shortcut[counter3] + iList[index][counter3].get( counter4 )] += weight
								* dList[index][counter3].get( counter4 );
					}
				}

				if( beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] != 0 )
				{
					Normalisation.logSumNormalisation( helpArray[index], 0, helpArray[index].length, helpArray[index], 0 );

					for( counter1 = 0; counter1 < shortcut[0]; counter1++ )
					{
						if( counter1 != counter3 )
						{
							cllGrad[index][counter1] -= weight * helpArray[index][counter1];
						}
						else
						{
							cllGrad[index][counter1] += weight * (1 - helpArray[index][counter1]);
						}
					}
					for( counter1 = 0; counter1 < cl; counter1++ )
					{
						if( counter1 != counter3 )
						{
							for( counter4 = 0; counter4 < iList[index][counter1].length(); counter4++ )
							{
								cllGrad[index][shortcut[counter1] + iList[index][counter1].get( counter4 )] -= weight
										* dList[index][counter1].get( counter4 ) * helpArray[index][counter1];
							}
						}
						else
						{
							for( counter4 = 0; counter4 < iList[index][counter1].length(); counter4++ )
							{
								cllGrad[index][shortcut[counter1] + iList[index][counter1].get( counter4 )] += weight
										* dList[index][counter1].get( counter4 ) * (1d - helpArray[index][counter1]);
							}
						}
					}
				}
			}
		}
	}

	protected double joinFunction() throws DimensionException, EvaluationException
	{		
		int i = 0;
		double cll = 0, ll = 0, lpr = 0, z = Double.NEGATIVE_INFINITY;
		for( ; i < helpArray.length; i++ )
		{
			ll += helpArray[i][0];
			cll += helpArray[i][1];
		}
		if( beta[LearningPrinciple.LIKELIHOOD_INDEX] != 0 )
		{
			for( i = 0; i < cl; i++ )
			{
				z = Normalisation.getLogSum( z, logClazz[i] + ((DifferentiableStatisticalModel)score[0][i]).getLogNormalizationConstant() );
			}
			ll -= (sum[cl] * z);
		}

		if( beta[LearningPrinciple.PRIOR_INDEX] != 0 )
		{
			lpr = prior.evaluateFunction( params );
		}
		// to avoid: 0*Infinity, ...
		if( beta[LearningPrinciple.LIKELIHOOD_INDEX] == 0 )
		{
			ll = 0;
		}
		if( beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] == 0 )
		{
			cll = 0;
		}


		double res = beta[LearningPrinciple.LIKELIHOOD_INDEX] * ll
			+ beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] * cll
			+ beta[LearningPrinciple.PRIOR_INDEX] * lpr;
	
		if( Double.isNaN( res ) || Double.isInfinite( res ) )
		{
			System.out.println( res + "\t= "
					+ beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] + " * " + cll
					+ " + " + beta[LearningPrinciple.LIKELIHOOD_INDEX] + " * " + ll					
					+ " + " + beta[LearningPrinciple.PRIOR_INDEX] + " * " + lpr );
			System.out.println( "params " + Arrays.toString( params ) );
			System.out.flush();
			throw new EvaluationException( "Evaluating the function gives: "
					+ beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] + " * " + cll
					+ " + " + beta[LearningPrinciple.LIKELIHOOD_INDEX] + " * " + ll
					+ " + " + beta[LearningPrinciple.PRIOR_INDEX] + " * " + lpr );
		}
		else if( norm )
		{
			// normalization
			return res / sum[cl];
		}
		else
		{
			return res;
		}
	}
	
	protected void evaluateFunction( int index, int startClass, int startSeq, int endClass, int endSeq ) throws EvaluationException
	{
		double cll = 0, ll = 0;
		int counter1, counter2, counter3 = startClass, start , end;

		Sequence s;
		for( ; counter3 <= endClass; counter3++ )
		{
			if( counter3 == startClass )
			{
				start = startSeq;
			}
			else
			{
				start = 0;
			}
			if( counter3 == endClass )
			{
				end = endSeq;
			}
			else
			{
				end = data[counter3].getNumberOfElements();
			}
			
			for( counter2 = start; counter2 < end; counter2++ )
			{
				s = data[counter3].getElementAt( counter2 );
				if( beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] != 0 )
				{
					for( counter1 = 0; counter1 < cl; counter1++ )
					{
						// class weight + class score
						helpArray[index][counter1] = logClazz[counter1] + score[index][counter1].getLogScoreFor( s, 0 );
					}
					cll += weights[counter3][counter2] * (helpArray[index][counter3] - Normalisation.getLogSum( helpArray[index] ));
				}
				else
				{
					helpArray[index][counter3] = logClazz[counter3] + score[index][counter3].getLogScoreFor( s, 0 );
				}
				ll += weights[counter3][counter2] * helpArray[index][counter3];
			}
		}

		helpArray[index][0] = ll;
		helpArray[index][1] = cll;
	}

	private void check() throws IllegalArgumentException {
		if( beta[LearningPrinciple.LIKELIHOOD_INDEX] != 0 )
		{
			for( int i =0; i < score[0].length; i++ ) {
				if( !( score[0][i] instanceof DifferentiableStatisticalModel ) ) {
					throw new IllegalArgumentException( "For evaluating the likelihood we the " );
				}
			}
		}
	}
	
	public void reset( DifferentiableSequenceScore[] funs ) throws Exception
	{
		for( int j, i = 0; i < cl; i++ )
		{
			score[0][i] = funs[i];
			for( j = 1; j < score.length; j++ )
			{
				score[j][i] = score[0][i].clone();
			}
			shortcut[i + 1] = shortcut[i] + score[0][i].getNumberOfParameters();
		}
		check();
		if( beta[LearningPrinciple.PRIOR_INDEX] > 0 && prior != null )
		{
			prior.set( freeParams, score[0] );
		}
		llGrad = new double[getNumberOfThreads()][shortcut[cl]];
		cllGrad = new double[llGrad.length][shortcut[cl]];
		prGrad = new double[shortcut[cl]];
	}
}
