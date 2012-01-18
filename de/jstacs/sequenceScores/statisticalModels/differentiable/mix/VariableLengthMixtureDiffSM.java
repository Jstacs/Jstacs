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

package de.jstacs.sequenceScores.statisticalModels.differentiable.mix;

import java.util.Arrays;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.sequenceScores.statisticalModels.differentiable.VariableLengthDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;


/**
 * This class implements a mixture of {@link VariableLengthDiffSM} by extending {@link MixtureDiffSM} and implementing the methods of {@link VariableLengthDiffSM}.
 * 
 * @author Jens Keilwagen
 */
public class VariableLengthMixtureDiffSM extends MixtureDiffSM implements VariableLengthDiffSM {

	/**
	 * This constructor creates a new {@link VariableLengthMixtureDiffSM}.
	 * 
	 * @param starts
	 *            the number of starts that should be done in an optimization
	 * @param plugIn
	 *            indicates whether the initial parameters for an optimization
	 *            should be related to the data or randomly drawn
	 * @param component
	 *            the {@link VariableLengthDiffSM}s
	 * 
	 * @throws CloneNotSupportedException
	 *             if an element of <code>component</code> could not be cloned
	 */
	public VariableLengthMixtureDiffSM( int starts, boolean plugIn, VariableLengthDiffSM... component ) throws CloneNotSupportedException {
		super( starts, plugIn, component );
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link VariableLengthMixtureDiffSM} out of a
	 * {@link StringBuffer}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public VariableLengthMixtureDiffSM( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.AbstractDifferentiableSequenceScore#getLogScoreFor(int, de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getLogScoreFor( Sequence seq, int start, int end ) {
		for( int i = 0; i < function.length; i++ ) {
			componentScore[i] = logHiddenPotential[i] + ((VariableLengthDiffSM)function[i]).getLogScoreFor( seq, start, end );
		}
		return Normalisation.getLogSum( componentScore );
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int start, int end, IntList indices, DoubleList partialDer ) {
		int i = 0, j = 0, k = paramRef.length - 1;
		k = paramRef[k] - paramRef[k - 1];
		for( ; i < function.length; i++ ) {
			iList[i].clear();
			dList[i].clear();
			componentScore[i] = logHiddenPotential[i] + ((VariableLengthDiffSM)function[i]).getLogScoreAndPartialDerivation( seq, start, end, iList[i], dList[i] );
		}
		double logScore = Normalisation.logSumNormalisation( componentScore, 0, function.length, componentScore, 0 );
		for( i = 0; i < function.length; i++ ) {
			for( j = 0; j < iList[i].length(); j++ ) {
				indices.add( paramRef[i] + iList[i].get( j ) );
				partialDer.add( componentScore[i] * dList[i].get( j ) );
			}
		}
		for( j = 0; j < k; j++ ) {
			indices.add( paramRef[i] + j );
			partialDer.add( componentScore[j] - ( isNormalized() ? hiddenPotential[j] : 0 ) );
		}
		return logScore;
	}

	public double getLogNormalizationConstant( int length ) {
		double n = Double.NEGATIVE_INFINITY;
		for( int i = 0; i < logHiddenPotential.length; i++ ) {
			n = Normalisation.getLogSum( n, logHiddenPotential[i] + ((VariableLengthDiffSM)function[i]).getLogNormalizationConstant( length ) );
		}
		return n;
	}

	public double getLogPartialNormalizationConstant( int parameterIndex, int length ) throws Exception {
		if( isNormalized() ) {
			return Double.NEGATIVE_INFINITY;
		} else {
			int[] ind = getIndices( parameterIndex );
			if( ind[0] == function.length ) {
				return logHiddenPotential[ind[1]] + ((VariableLengthDiffSM) function[ind[1]]).getLogNormalizationConstant( length );
			} else {
				return logHiddenPotential[ind[0]] + ((VariableLengthDiffSM) function[ind[0]]).getLogPartialNormalizationConstant( ind[1], length );
			}
		}
	}

	public void setStatisticForHyperparameters( int[] length, double[] weight ) throws Exception {
		double[] w = new double[getNumberOfComponents()];
		for( int i = 0; i < w.length; i++ ) {
			w[i] = getHyperparameterForHiddenParameter( i );
		}
		Normalisation.sumNormalisation( w );
		double[] nw = new double[weight.length];
		for( int i = 0; i < function.length; i++ ) {
			for( int j = 0; j < nw.length; j++ ) {
				nw[j] = weight[j] * w[i];
			}
			( (VariableLengthDiffSM)function[i] ).setStatisticForHyperparameters( length, nw );
		}
	}
}
