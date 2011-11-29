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
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.models.hmm;

import java.text.NumberFormat;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.Sequence;
import de.jstacs.data.WrongLengthException;

/**
 * This interface declares the methods of any state used in a hidden Markov model.
 * 
 * @author Jan Grau, Jens Keilwagen, Michael Scharfe
 * 
 * @see AbstractHMM
 */
public interface State  {
	
	/**
	 * This method returns the logarithm of the score for a given sequence with given start and end position.
	 * The method provides the basics for quantifying emission probability/density.
	 * 
	 * 
	 * @param startPos the start position within the sequence(s)
	 * @param endPos the end position within the sequence(s)
	 * @param seq the sequence(s)
	 * 
	 * @return the logarithm of score for the given sequence(s)
	 * 
	 * @throws WrongLengthException if the length can not be modeled 
	 * @throws OperationNotSupportedException if the reverse complement of the sequence can not be computed 
	 */
	public double getLogScoreFor(int startPos, int endPos, Sequence seq) throws WrongLengthException, OperationNotSupportedException ;
	
	/**
	 * This method returns a {@link String} representation of the node options that
	 * can be used in <i>Graphviz</i> to create the node for this state.
	 * 
	 * @param weight for state
	 * @param maxWeight the maximal weight for the state
	 * @param nf the {@link NumberFormat} for the output, can be <code>null</code>
	 * 
	 * @return {@link String} representation of the state
	 */
	public String getGraphvizNodeOptions( double weight, double maxWeight, NumberFormat nf );
	
	/**
	 * This method returns whether a state is silent or not.
	 * 
	 * @return <code>true</code> if the state is silent, otherwise <code>false</code>
	 * 
	 * @see de.jstacs.models.hmm.states.emissions.SilentEmission
	 */
	public boolean isSilent();
	
	/**
	 * Returns a {@link String} indentifier for the type of this state.
	 * @return the identifier
	 */
	public String getEmissionType();
}
