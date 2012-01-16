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

package de.jstacs.sequenceScores.statisticalModels.differentiable.mix.motif;

import de.jstacs.NonParsableException;
import de.jstacs.data.DataSet;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This scoring function implements a uniform distribution for positions.
 * The class has no parameters, so the distribution does not change and it is possible to save parameters in an optimization.
 */
public final class UniformDurationDiffSM extends DurationDiffSM
{
	private double logP;
	
	/**
	 * This is the main constructor that creates an instance for the given interval.
	 * 
	 * @param min the minimal value
	 * @param max the maximal value
	 */
	public UniformDurationDiffSM( int min, int max )
	{
		this( min, max, 0 );
	}
	
	/**
	 * This is the main constructor that creates an instance for the given interval and given ESS.
	 * 
	 * @param min the minimal value
	 * @param max the maximal value
	 * @param ess the equivalent sample size (used for the class probability)
	 */
	public UniformDurationDiffSM( int min, int max, double ess )
	{
		super( min, max, ess );
		computeLogP();
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link UniformDurationDiffSM} out of a {@link StringBuffer}.
	 * 
	 * @param b
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public UniformDurationDiffSM( StringBuffer b ) throws NonParsableException
	{
		super( b );
		computeLogP();
	}
	
	private void computeLogP()
	{
		logP = -Math.log( delta + 1 );
	}

	public String getInstanceName()
	{
		return "uniform";
	}

	public int getNumberOfParameters()
	{
		return 0;
	}

	public void setParameters( double[] params, int start )
	{
	}

	public void initializeFunction( int index, boolean meila, DataSet[] data, double[][] weights )
	{
		// does nothing
	}
	
	protected String getRNotation( String distributionName )
	{
		return "l = " + min + ":" + max + "; n = length(l); " + distributionName + " = rep(1,n)/n;";
	}
	
	public double getLogPriorTerm()
	{
		// since the normalization constant does not depend on any parameter,
		// it is constant and therefore left out
		return 0;
	}

	public void addGradientOfLogPriorTerm( double[] grad, int start ){}

	public double getLogScore( int... values )
	{
		return logP;
	}

	public double getLogScoreAndPartialDerivation( IntList indices, DoubleList partialDer, int... values )
	{
		return logP;
	}

	/**
	 * This method draws from the distribution and returns the result in the given array.
	 * 
	 * @param positions an array for the result.
	 */
	public void drawPosition( int[] positions )
	{
		positions[0] = min + r.nextInt( delta + 1 );		
	}

	public double[] getCurrentParameterValues() throws Exception
	{
		return new double[0];
	}

	public boolean isInitialized()
	{
		return true;
	}
	
	public boolean isNormalized()
	{
		return true;
	}
	
	public void initializeFunctionRandomly(boolean freeParams) throws Exception
	{	
	}

	public void initializeUniformly(){}

	public void adjust( int[] length, double[] weight ){}
}
