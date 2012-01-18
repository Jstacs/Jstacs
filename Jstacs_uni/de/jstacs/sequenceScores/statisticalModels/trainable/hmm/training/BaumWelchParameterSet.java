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

package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training;

import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;

/**
 * This class implements an {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.HMMTrainingParameterSet} for the Baum-Welch training of an {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM}.
 * 
 * @author Jens Keilwagen
 */
public class BaumWelchParameterSet extends MultiThreadedTrainingParameterSet {

	/**
	 * This is the empty constructor that can be used to fill the parameters after creation.
	 */
	public BaumWelchParameterSet() {}

	/**
	 * This constructor can be used to create an instance with specified parameters.
	 * 
	 * @param starts the number of different starts
	 * @param tc the termination condition for stopping the algorithm
	 * @param threads the number of threads that should be used during optimization
	 * 
	 * @throws Exception if the {@link BaumWelchParameterSet} could not be created
	 */
	public BaumWelchParameterSet( int starts, AbstractTerminationCondition tc, int threads ) throws Exception {
		super( starts, tc, threads );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link BaumWelchParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link BaumWelchParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public BaumWelchParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	public boolean hasDefaultOrIsSet() {
		try {
			return super.hasDefaultOrIsSet() && getTerminationCondition().isSimple();
		} catch (NotInstantiableException e) {
			throw new RuntimeException( e );
		}
	}
}
