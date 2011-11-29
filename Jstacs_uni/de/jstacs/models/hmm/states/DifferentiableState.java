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

package de.jstacs.models.hmm.states;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.Sequence;
import de.jstacs.data.WrongLengthException;
import de.jstacs.models.hmm.State;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This interface declares a method that allows to evaluate the gradient which is essential for numerical optimization. 
 * 
 * @author Jens Keilwagen
 */
public interface DifferentiableState extends State {

	/**
	 * This method allows to compute the logarithm of the score and the gradient for the given subsequences.
	 * 
	 * @param startPos the start position (inclusive)
	 * @param endPos the end position (inclusive)
	 * @param indices a list for the parameter indices
	 * @param partDer a list for the partial derivations
	 * @param seq the {@link Sequence}(s)
	 * 
	 * @return the logarithm of the score
	 * 
	 * @throws WrongLengthException if the length can not be modeled 
	 * @throws OperationNotSupportedException if the reverse complement of the sequence can not be computed
	 */
	public double getLogScoreAndPartialDerivation( int startPos, int endPos, IntList indices, DoubleList partDer, Sequence seq ) throws OperationNotSupportedException, WrongLengthException;
}
