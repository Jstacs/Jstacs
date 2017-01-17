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

package projects.dream2016.mix;

import java.util.Arrays;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.classifiers.differentiableSequenceScoreBased.AbstractMultiThreadedOptimizableFunction;
import de.jstacs.data.DataSet;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;


/**
 * This class implements maximum supervised posterior (MSP) based on an {@link OptimizableClassifier}. 
 * 
 * @author Jens
 */
public class MSPClassifierObjective extends AbstractMultiThreadedOptimizableFunction {

	private OptimizableClassifier[] optClassifiers;
	private double[][] grad;
	private double[] value;
	private IntList[] indices;
	private DoubleList[] partDer;
	
	public MSPClassifierObjective( int threads, OptimizableClassifier optCl, DataSet[] data, double[][] weights, boolean norm )
																										throws IllegalArgumentException {
		super( threads, data, weights, norm, false );
		
		optClassifiers = new OptimizableClassifier[threads];
		this.optClassifiers[0] = optCl;
		
		value = new double[threads];
		grad = new double[threads][];
		
		indices = new IntList[threads];
		partDer = new DoubleList[threads];
		for( int t = 0; t < indices.length; t++ ) {
			indices[t] = new IntList();
			partDer[t] = new DoubleList();
		}
	}

	@Override
	protected void evaluateFunction( int index, int startClass, int startSeq, int endClass, int endSeq ) throws EvaluationException {
		int seqIndex, dataSet = startClass, start, end;
		value[index] = 0;
		for( ; dataSet <= endClass; dataSet++ )
		{
			if( dataSet == startClass )
			{
				start = startSeq;
			}
			else
			{
				start = 0;
			}
			if( dataSet == endClass )
			{
				end = endSeq;
			}
			else
			{
				end = data[dataSet].getNumberOfElements();
			}
			
			for( seqIndex = start; seqIndex < end; seqIndex++ )
			{
				double v= weights[dataSet][seqIndex] * optClassifiers[index].getLogProb(dataSet, data[dataSet].getElementAt( seqIndex ));;
				if( Double.isNaN(v) ) {
					System.out.println("PROBLEM: NaN");
					System.out.println("Classifier:");
					System.out.println( optClassifiers[index] );
					System.out.println("Sequence:");
					System.out.println( data[dataSet].getElementAt(seqIndex) );
					try {
						System.out.println("Parameter:");
						System.out.println( Arrays.toString(optClassifiers[0].getCurrentParameterValues( KindOfParameter.LAST ) ) );
					} catch( Exception e ) {
						
					}
					throw new IllegalArgumentException( "NaN" );
				}
				value[index] += v;
			}
		}
	}

	@Override
	protected double joinFunction() throws EvaluationException, DimensionException {
		double res = 0;
		
		//sum conditional likelihood of a mixture classifier
		for( int t = 0; t < value.length; t++ ) {
			res += value[t];
		}
		
		//add prior
		res += optClassifiers[0].getLogPriorTerm();
		
		if( Double.isNaN( res ) ) {
			throw new EvaluationException( "Error in evaluation: " + res );
		}
		
		//normalize
		if( norm ) {
			res /= sum[cl];
		}
		//XXX System.out.println( res );
		return res;
	}
	
	@Override
	protected void evaluateGradientOfFunction( int index, int startClass, int startSeq, int endClass, int endSeq ) {
		int seqIndex, dataSet = startClass, start, end, h;
		Arrays.fill( grad[index], 0 );
		for( ; dataSet <= endClass; dataSet++ )
		{
			if( dataSet == startClass )
			{
				start = startSeq;
			}
			else
			{
				start = 0;
			}
			if( dataSet == endClass )
			{
				end = endSeq;
			}
			else
			{
				end = data[dataSet].getNumberOfElements();
			}
			
			for( seqIndex = start; seqIndex < end; seqIndex++ )
			{
				indices[index].clear();
				partDer[index].clear();
				optClassifiers[index].getLogProbAndPartialDerivations(dataSet, data[dataSet].getElementAt( seqIndex ), indices[index], partDer[index] );
				for( h = 0; h < indices[index].length(); h++ ) {
					grad[index][indices[index].get(h)] += weights[dataSet][seqIndex] * partDer[index].get(h);
				}
			}
		}
	}

	@Override
	protected double[] joinGradients() throws EvaluationException {
		double[] res = new double[grad[0].length];
		//sum gradients of the conditional likelihood of a mixture classifier
		for( int i = 0; i < res.length; i++ ) {
			for( int t = 0; t < grad.length; t++ ) {
				res[i] += grad[t][i];
			}
		}
		
		//add gradient of the prior
		optClassifiers[0].addGradient( res, 0 );
		
		//normalize
		if( norm ) {
			for( int i = 0; i < res.length; i++ ) {
				res[i] /= sum[cl];
			}
		}

		/*
		System.out.println( Arrays.toString( res ) );
		double l = 0;
		for( int i = 0; i < res.length; i++ ) {
			l += (res[i] * res[i]);
		}
		System.out.println("grad-length: " + l);
		/**/
		return res;
	}

	@Override
	public void getParameters( KindOfParameter kind, double[] erg ) throws Exception {
		double[] params = optClassifiers[0].getCurrentParameterValues( kind );
		System.arraycopy( params, 0, erg, 0, params.length );
	}

	@Override
	public void reset() throws Exception {
		optClassifiers[0].reset();
		
		for( int j = 1; j < optClassifiers.length; j++ ) {
			optClassifiers[j] = optClassifiers[0].clone();
			optClassifiers[j].reset();
		}
		
		for( int t = 0; t < grad.length; t++ ) {
			grad[t] = new double[getDimensionOfScope()];
		}
	}
	
	@Override
	protected void setThreadIndependentParameters() throws DimensionException {
	}
	
	@Override
	protected void setParams( int index ) throws DimensionException {
		try {
			optClassifiers[index].setParameters( params, 0 );
		} catch (Exception e) {
			DimensionException d = new DimensionException();
			d.setStackTrace( e.getStackTrace() );
			throw d;
		}
	}

	public int getDimensionOfScope() {
		return optClassifiers[0].getNumberOfParameters();
	}
}