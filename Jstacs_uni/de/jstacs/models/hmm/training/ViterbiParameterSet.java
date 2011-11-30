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

package de.jstacs.models.hmm.training;

import de.jstacs.NonParsableException;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;

/**
 * This class implements an {@link de.jstacs.parameters.ParameterSet} for the viterbi training of an {@link de.jstacs.models.hmm.AbstractHMM}.
 * 
 * @author Jens Keilwagen
 */
public class ViterbiParameterSet extends MultiThreadedTrainingParameterSet {

	/**
	 * This is the empty constructor that can be used to fill the parameters after creation.
	 */
	public ViterbiParameterSet() {}

	/**
	 * This constructor can be used to create an instance with specified parameters.
	 * 
	 * @param starts the number of different starts
	 * @param tc the termination condition for stopping the algorithm
	 * @param threads the number of threads that should be used during optimization
	 * 
	 * @throws Exception if the {@link ViterbiParameterSet} could not be created
	 */
	public ViterbiParameterSet( int starts, AbstractTerminationCondition tc, int threads ) throws Exception {
		super( starts, tc, threads );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ViterbiParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ViterbiParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public ViterbiParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	public boolean hasDefaultOrIsSet() {
		try {
			return super.hasDefaultOrIsSet() && getTerminantionCondition().isSimple();
		} catch (NotInstantiableException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return false;
		}
	}
}
