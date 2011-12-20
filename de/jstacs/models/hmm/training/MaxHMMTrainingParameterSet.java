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

import de.jstacs.NonParsableException;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition.AbstractTerminationConditionParameterSet;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;
import de.jstacs.models.hmm.HMMTrainingParameterSet;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.utils.SubclassFinder;

/**
 * This class is the super class for any {@link HMMTrainingParameterSet} that
 * is used for a maximizing training algorithm of a hidden Markov model.
 *  
 * @author Jens Keilwagen
 */
public abstract class MaxHMMTrainingParameterSet extends HMMTrainingParameterSet {

	/**
	 * This is the empty constructor that can be used to fill the parameters after creation.
	 */
	protected MaxHMMTrainingParameterSet() {
		addParameters();
	}

	/**
	 * This constructor can be used to create an instance with specified parameters.
	 * 
	 * @param starts the number of different starts
	 * @param tc the termination condition for stopping the algorithm
	 * 
	 * @throws Exception if this {@link MaxHMMTrainingParameterSet} could not be created
	 */
	protected MaxHMMTrainingParameterSet( int starts, AbstractTerminationCondition tc ) throws Exception {
		super( starts );
		addParameters();
		parameters.get(1).setValue(tc.getCurrentParameterSet());
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link MaxHMMTrainingParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MaxHMMTrainingParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	protected MaxHMMTrainingParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}


	private void addParameters() {
		try{
			parameters.add(
					SubclassFinder.getSelectionParameter(
							AbstractTerminationConditionParameterSet.class,
							AbstractTerminationCondition.class.getPackage().getName(),
							"termination condition",
							"the terminantion condition for stopping the training algorithm",
							true
					)
			);
		}catch(Exception hopefullydoesnothappen){
			RuntimeException ex = new RuntimeException( hopefullydoesnothappen );
			ex.setStackTrace( hopefullydoesnothappen.getStackTrace() );
			throw ex;
		}
	}

	/**
	 * This method returns the {@link AbstractTerminationCondition} for stopping the training, e.g., if the
	 * difference of the scores between two iterations is smaller than a given
	 * threshold the training is stopped.
	 * 
	 * @return the {@link AbstractTerminationCondition} for stopping the training
	 * @throws NotInstantiableException if the {@link AbstractTerminationCondition} could not be created from its {@link de.jstacs.parameters.ParameterSet}
	 */
	public AbstractTerminationCondition getTerminantionCondition() throws NotInstantiableException {
		return (AbstractTerminationCondition)(((InstanceParameterSet)getParameterForName( "termination condition" ).getValue()).getInstance());
	}
}
