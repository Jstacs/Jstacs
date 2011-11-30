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

package de.jstacs.models.hmm.training;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * This class is the super class for any {@link MaxHMMTrainingParameterSet} that
 * is used for a multi-threaded maximizing training algorithm of a hidden Markov model.
 *  
 * @author Jens Keilwagen
 */
public abstract class MultiThreadedTrainingParameterSet extends MaxHMMTrainingParameterSet {

	/**
	 * This is the empty constructor that can be used to fill the parameters after creation.
	 */
	protected MultiThreadedTrainingParameterSet() {}

	/**
	 * This constructor can be used to create an instance with specified parameters.
	 * 
	 * @param starts the number of different starts
	 * @param tc the termination condition for stopping the algorithm
	 * @param threads the number of threads that should be used during optimization
	 * 
	 * @throws Exception if this {@link MultiThreadedTrainingParameterSet} could not be created
	 */
	protected MultiThreadedTrainingParameterSet(int starts, AbstractTerminationCondition tc, int threads) throws Exception {
		super(starts, tc);
		parameters.get(2).setDefault( threads );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link MultiThreadedTrainingParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MultiThreadedTrainingParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	protected MultiThreadedTrainingParameterSet(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.models.hmm.training.MaxHMMTrainingParameterSet#loadParameters()
	 */
	@Override
	protected void loadParameters() throws Exception {
		super.loadParameters();
		parameters.add( new SimpleParameter( DataType.INT, "Threads", "the number of threads that should be used during optimization", true, new NumberValidator<Integer>(1, Integer.MAX_VALUE), 1 ) );
	}
	
	/**
	 * This method returns the number of threads that should be used during optimization.
	 * 
	 * @return the number of threads that should be used during optimization
	 */
	public int getNumberOfThreads() {
		return (Integer)getParameterAt( 2 ).getValue();
	}
}
