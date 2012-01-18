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

package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.sequences.Sequence;

/**
 * This class implements method that allows to fill a statistic, which is used to estimate the parameters of a state during, for instance, the Baum-Welch training.
 * All other methods that allow, for instance, to reset the statistic, to estimate the parameters from the statistic, ... are not defined in this
 * interface to allow parameter sharing between states, and, hence, have to be defined on their own.
 * 
 * @author Jens Keilwagen
 */
public interface TrainableState extends State {

	/**
	 * This method allows to add a certain <code>weight</code> to the sufficient statistic of the parameters that
	 * are used for scoring the specific subsequence(s).
	 * 
	 * @param startPos the start position
	 * @param endPos the end position
	 * @param weight the weight which will be added to the sufficient statistic
	 * @param seq the {@link Sequence}(s)
	 * 
	 * @throws OperationNotSupportedException if the reverse complement of the sequence can not be computed
	 */
	public void addToStatistic( int startPos, int endPos, double weight, Sequence seq ) throws OperationNotSupportedException;
}
