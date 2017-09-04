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

import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.utils.Normalisation;

/**
 * This class implements the the following function
 * 
 * {@latex.ilb \\[f(\\underline{\\lambda}|C,D,\\underline{w},\\underline{\\alpha},\\underline{\\beta})
 * := \\left(\\sum_c \\sum_n w_{c,n} \\left(\\beta_0 \\log p(c|d_n,\\underline{\\lambda}) + \\beta_1 \\log p(c,d_n|\\underline{\\lambda}) \\right) \\right) + \\beta_2 \\log p(\\underline{\\lambda}|\\underline{\\alpha})
 * \\]}
 * where {@latex.inline $w_{c,n}$} is the weight for sequence {@latex.inline $d_n$} and class {@latex.inline $c$}.
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
 * It can be used to maximize the parameters. This can also be done with {@link LogGenDisMixFunction}.
 * However, we implemented this class to allow a faster function and gradient evaluation leading to a faster optimization.
 * This becomes especially interesting if the number of classes increases.
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
 * It is very important for this class that the {@link de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#clone()} method works correctly, since each thread works on its own clones. 
 * 
 * @author Jens Keilwagen
 */
public class OneDataSetLogGenDisMixFunction extends LogGenDisMixFunction
{
	/**
	 * The constructor for creating an instance that can be used in an {@link de.jstacs.algorithms.optimization.Optimizer}.
	 * 
	 * @param threads the number of threads used for evaluating the function and determining the gradient of the function
	 * @param score an array containing the {@link DifferentiableSequenceScore}s that are used for determining the sequences scores;
	 * 			if the weight <code>beta[LearningPrinciple.LIKELIHOOD_INDEX]</code> is positive all elements of <code>score</code> have to be {@link de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel}
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
	public OneDataSetLogGenDisMixFunction( int threads, DifferentiableSequenceScore[] score, DataSet data, double[][] weights,
			LogPrior prior, double[] beta, boolean norm, boolean freeParams ) throws IllegalArgumentException
	{
		super(threads, score, new DataSet[]{data}, weights, prior, beta, norm, freeParams);
	}
	
	public void setDataAndWeights( DataSet[] data, double[][] weights ) throws IllegalArgumentException {
		if( data.length != 1 || weights == null || weights.length != cl ) {
			throw new IllegalArgumentException( "The dimension of the data set or weights (array) is not correct."  );
		}
		this.data = data;
		this.weights = weights;
		sum[cl] = 0;
		int i = 0, j;
		for( ; i < cl; i++ ) {
			sum[i] = 0;
			if( data[0].getNumberOfElements() != weights[i].length ) {
				throw new IllegalArgumentException( "The dimension of the " + i + "-th weights (array) is not correct (" + data[0].getNumberOfElements() + " vs. " + weights[i].length + ")."  );
			}
			for( j = 0; j < weights[i].length; j++ ) {
				sum[i] += weights[i][j];
			}
			sum[cl] += sum[i];
		}
		if( worker != null ) {
			prepareThreads();
		}
	}
	
	public DataSet[] getData() {
		DataSet[] d = new DataSet[weights.length];
		Arrays.fill( d, data[0] );
		return d;
	}
	
	protected void evaluateGradientOfFunction( int index, int startClass, int startSeq, int endClass, int endSeq )
	{
		Arrays.fill( llGrad[index], 0 );
		Arrays.fill( cllGrad[index], 0 );

		double weight;
		int counter1, counter2, counter3, counter4 = 0;

		Sequence s;
		for( counter2 = startSeq; counter2 < endSeq; counter2++ )
		{
			s = data[0].getElementAt( counter2 );
			
			for( counter1 = 0; counter1 < cl; counter1++ )
			{
				iList[index][counter1].clear();
				dList[index][counter1].clear();
				helpArray[index][counter1] = logClazz[counter1]
						+ score[index][counter1].getLogScoreAndPartialDerivation( s, 0, iList[index][counter1],
								dList[index][counter1] );
			}
			
			Normalisation.logSumNormalisation( helpArray[index], 0, helpArray[index].length, helpArray[index], 0 );
			for( counter3 = 0; counter3 < cl; counter3++ ) {
				weight = weights[counter3][counter2];
			
			
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

	protected void evaluateFunction( int index, int startClass, int startSeq, int endClass, int endSeq ) throws EvaluationException
	{
		double cll = 0, ll = 0, offset = 0;
		int counter1, counter2, counter3;

		Sequence s;		
		for( counter2 = startSeq; counter2 < endSeq; counter2++ )
		{
			s = data[0].getElementAt( counter2 );
			for( counter1 = 0; counter1 < cl; counter1++ )
			{
				// class weight + class score
				helpArray[index][counter1] = logClazz[counter1] + score[index][counter1].getLogScoreFor( s, 0 );
			}
			if( beta[LearningPrinciple.CONDITIONAL_LIKELIHOOD_INDEX] != 0 ) {
				offset = Normalisation.getLogSum( helpArray[index] );
			}
			for( counter3 = 0; counter3 < cl; counter3++ ) {
				double part = weights[counter3][counter2] * (helpArray[index][counter3] - offset);
				cll+=part;
				
				/*if( Double.isNaN( cll ) || Double.isInfinite( cll ) ) {
					System.out.println("counters: " + index + ", " + counter3 + ", " + counter2);
					System.out.println("likelihoods: " + Arrays.toString(helpArray[index]));
					System.out.println("weight: " + weights[counter3][counter2] + "\tcll: " + (helpArray[index][counter3] - offset) );
					System.out.println("part of global cll: " + part);
					System.out.println(s);
					SequenceAnnotation[] annot = s.getAnnotation();
					for( int z = 0; z < annot.length; z++ ) {
						System.out.println(annot[z]);
					}
					System.out.println();
					System.out.println(score[index][0]);
					System.out.println();
					throw new EvaluationException("CLL will be: " + cll );
				}/**/
				
				ll += weights[counter3][counter2] * helpArray[index][counter3];
			}
		}

		helpArray[index][0] = ll;
		helpArray[index][1] = cll;
	}
}
