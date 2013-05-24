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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete;

import java.text.NumberFormat;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;

/**
 * This class implements a simple discrete emission without any condition.
 * 
 * @author Jens Keilwagen, Michael Scharfe, Jan Grau
 */
public class DiscreteEmission extends AbstractConditionalDiscreteEmission {
	
	/**
	 * This is a simple constructor for a {@link DiscreteEmission} based on the equivalent sample size.
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param ess the equivalent sample size (ess) of this emission that is equally distributed over all parameters
	 * 
	 * @see #DiscreteEmission(AlphabetContainer, double[])
	 */
	public DiscreteEmission( AlphabetContainer con, double ess ) {
		super( con, 1, ess );
	}

	/**
	 * This is a simple constructor for a {@link DiscreteEmission} defining the individual hyper parameters.
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param hyperParams the individual hyper parameters for each parameter
	 */
	public DiscreteEmission( AlphabetContainer con, double[] hyperParams ) {
		super( con, new double[][]{hyperParams} );
	}

	/**
	 * Creates a {@link DiscreteEmission} from its XML representation.
	 * @param xml the XML representation.
	 * @throws NonParsableException if the XML representation could not be parsed
	 */
	public DiscreteEmission( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	protected int getConditionIndex( boolean forward, int seqPos, Sequence seq ) {
		return 0;
	}

	@Override
	public String toString( NumberFormat nf ) {
		String res = "";
		DiscreteAlphabet abc = (DiscreteAlphabet) con.getAlphabetAt( 0 );
		for( int i = 0; i < probs[0].length; i++ ) {
			res += "P(X=" + abc.getSymbolAt( i ) + ") = " + nf.format( probs[0][i] ) + "\t";
		}
		return res+"\n";
	}
}