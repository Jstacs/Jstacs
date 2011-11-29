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

package de.jstacs.scoringFunctions.mix;

import de.jstacs.NonParsableException;
import de.jstacs.data.Sequence;
import de.jstacs.scoringFunctions.VariableLengthScoringFunction;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;


/**
 * This class implements a mixture of {@link VariableLengthScoringFunction} by extending {@link MixtureScoringFunction} and implementing the methods of {@link VariableLengthScoringFunction}.
 * 
 * @author Jens Keilwagen
 */
public class VariableLengthMixtureScoringFunction extends MixtureScoringFunction implements VariableLengthScoringFunction {

	/**
	 * This constructor creates a new {@link VariableLengthMixtureScoringFunction}.
	 * 
	 * @param starts
	 *            the number of starts that should be done in an optimization
	 * @param plugIn
	 *            indicates whether the initial parameters for an optimization
	 *            should be related to the data or randomly drawn
	 * @param component
	 *            the {@link VariableLengthScoringFunction}s
	 * 
	 * @throws CloneNotSupportedException
	 *             if an element of <code>component</code> could not be cloned
	 */
	public VariableLengthMixtureScoringFunction( int starts, boolean plugIn, VariableLengthScoringFunction... component ) throws CloneNotSupportedException {
		super( starts, plugIn, component );
	}

	/**
	 * This is the constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link VariableLengthMixtureScoringFunction} out of a
	 * {@link StringBuffer}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation could not be parsed
	 */
	public VariableLengthMixtureScoringFunction( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	public double getLogScore( Sequence seq, int start, int length ) {
		for( int i = 0; i < function.length; i++ ) {
			componentScore[i] = logHiddenPotential[i] + ((VariableLengthScoringFunction)function[i]).getLogScore( seq, start, length );
		}
		return Normalisation.getLogSum( componentScore );
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int start, int length, IntList indices, DoubleList partialDer ) {
		int i = 0, j = 0, k = paramRef.length - 1;
		k = paramRef[k] - paramRef[k - 1];
		for( ; i < function.length; i++ ) {
			iList[i].clear();
			dList[i].clear();
			componentScore[i] = logHiddenPotential[i] + ((VariableLengthScoringFunction)function[i]).getLogScoreAndPartialDerivation( seq, start, length, iList[i], dList[i] );
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
			n = Normalisation.getLogSum( n, logHiddenPotential[i] + ((VariableLengthScoringFunction)function[i]).getLogNormalizationConstant( length ) );
		}
		return n;
	}

	public double getLogPartialNormalizationConstant( int parameterIndex, int length ) throws Exception {
		if( isNormalized() ) {
			return Double.NEGATIVE_INFINITY;
		} else {
			int[] ind = getIndices( parameterIndex );
			if( ind[0] == function.length ) {
				return logHiddenPotential[ind[1]] + ((VariableLengthScoringFunction) function[ind[1]]).getLogNormalizationConstant( length );
			} else {
				return logHiddenPotential[ind[0]] + ((VariableLengthScoringFunction) function[ind[0]]).getLogPartialNormalizationConstant( ind[1], length );
			}
		}
	}

	public void setStatisticForHyperparameters( int[] length, double[] weight ) throws Exception {
		// TODO Auto-generated method stub
		throw new RuntimeException( "not implemented yet" );
	}
}
