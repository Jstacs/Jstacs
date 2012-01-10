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

package de.jstacs.classifier.differentiableSequenceScoreBased.logPrior;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.differentiableStatisticalModels.DifferentiableSequenceScore;
import de.jstacs.differentiableStatisticalModels.DifferentiableStatisticalModel;
import de.jstacs.utils.Normalisation;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * This class implements a composite prior that can be used for DifferentiableStatisticalModel. The prior for each
 * DifferentiableStatisticalModel should be the (transformed) prior of the corresponding generative model. So e.g. for a
 * PWM one should use an product of Dirichlets. The prior is more or less implemented in the
 * DifferentiableStatisticalModel. For the class variables the prior uses a (transformed) Dirichlet with hyperparameters
 * equal to the ESS of the classes.
 * 
 * <br>
 * <br>
 * 
 * If this class uses only the free parameters the class implements a real prior the is normalized to 1. If it used all
 * parameters the function does not have (and is in general) not normalized to 1. Fortunately this is no problem, since
 * it can be shown the it makes no difference in the optimization.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see DifferentiableStatisticalModel
 * @see DifferentiableStatisticalModel#addGradientOfLogPriorTerm(double[], int)
 * @see DifferentiableStatisticalModel#getLogPriorTerm()
 */
public class CompositeLogPrior extends LogPrior
{

	private DifferentiableStatisticalModel[] function;

	private double fullEss, logGammaSum;

	private double[] ess, classPars, logPartNorm;

	private boolean freeParameters;

	/**
	 * The main constructor.
	 */
	public CompositeLogPrior(){}

	/**
	 * The constructor for the {@link de.jstacs.Storable} interface.
	 * 
	 * @param xml the StringBuffer
	 */
	public CompositeLogPrior( StringBuffer xml )
	{
		this();
	}

	public void set( boolean freeParameters, DifferentiableSequenceScore... funs ) throws Exception
	{
		function = new DifferentiableStatisticalModel[funs.length];
		ess = new double[funs.length];
		classPars = new double[funs.length];
		logPartNorm = new double[funs.length];
		fullEss = 0;
		logGammaSum = 0;
		for( int i = 0; i < funs.length; i++ )
		{
			if( !(funs[i] instanceof DifferentiableStatisticalModel) )
			{
				throw new Exception( "Only DifferentiableStatisticalModel allowed." );
			}
			else
			{
				function[i] = (DifferentiableStatisticalModel) funs[i];
			}
			ess[i] = function[i].getESS();
			if( ess[i] == 0 ) {
				throw new IllegalArgumentException( "The ess of the function " + i + " is zero, but should be positive." );
			}
			fullEss += ess[i];
			logGammaSum -= Gamma.logOfGamma(ess[i]);
		}
		logGammaSum += Gamma.logOfGamma(fullEss);
		this.freeParameters = freeParameters;
	}

	public void addGradientFor( double[] params, double[] grad ) throws EvaluationException
	{
		try
		{
			int start = 0, j = function.length - (freeParameters ? 1 : 0), k = 0, l;
			for( ; k < j; k++ )
			{
				classPars[k] = params[k];
				logPartNorm[k] = classPars[k] + function[k].getLogNormalizationConstant();
			}
			if( freeParameters )
			{
				classPars[j] = 0;
				logPartNorm[j] = function[j].getLogNormalizationConstant();
			}
			double fullNorm = Normalisation.logSumNormalisation( logPartNorm );
			for( start = 0; start < j; start++ )
			{
				grad[start] += ess[start] - fullEss * logPartNorm[start];
			}

			for( j = 0 ; j < function.length; j++ )
			{
				function[j].addGradientOfLogPriorTerm( grad, start );
				k = function[j].getNumberOfParameters();
				for( l = 0; l < k; l++, start++ )
				{
					grad[start] -= (fullEss * Math.exp( classPars[j] + function[j].getLogPartialNormalizationConstant( l ) - fullNorm ) );
				}
			}
		}
		catch( Exception e )
		{
			e.printStackTrace();
			throw new EvaluationException( e.getMessage() );
		}
	}

	public double evaluateFunction( double[] x ) throws DimensionException, EvaluationException
	{
		try
		{
			double logNorm = Double.NEGATIVE_INFINITY, logProductPart = 0;
			int i = 0, j =  function.length - (freeParameters ? 1 : 0);
			for( ; i < j; i++ )
			{
				//TODO set parameters!?
				logNorm = Normalisation.getLogSum( logNorm, x[i] + function[i].getLogNormalizationConstant() );
				logProductPart += x[i] * ess[i] + function[i].getLogPriorTerm();
			}
			if( freeParameters )
			{
				logNorm = Normalisation.getLogSum( logNorm, function[j].getLogNormalizationConstant() );
				logProductPart += function[j].getLogPriorTerm();
			}
			return logGammaSum - fullEss * logNorm + logProductPart;
		}
		catch( Exception e )
		{
			e.printStackTrace();
			EvaluationException eva = new EvaluationException( e.getCause().getMessage() );
			eva.setStackTrace( e.getStackTrace() );
			throw eva;
		}
	}

	public int getDimensionOfScope()
	{
		int current, all = function.length - (freeParameters?1:0);
		for( int i = 0; i < function.length; i++ )
		{
			current = function[i].getNumberOfParameters();
			if( current == UNKNOWN )
			{
				return UNKNOWN;
			}
			else
			{
				all += current;
			}
		}
		return all;
	}

	public CompositeLogPrior getNewInstance() throws CloneNotSupportedException
	{
		return new CompositeLogPrior();
	}

	public StringBuffer toXML()
	{
		return new StringBuffer( 1 );
	}

	public String getInstanceName()
	{
		return "Composite log prior";
	}
}
