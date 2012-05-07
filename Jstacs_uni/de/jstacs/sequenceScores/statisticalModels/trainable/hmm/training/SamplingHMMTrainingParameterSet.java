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

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.sampling.AbstractBurnInTest;
import de.jstacs.sampling.AbstractBurnInTestParameterSet;
import de.jstacs.utils.SubclassFinder;

/**
 * This class contains the parameters for training training an {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM} using a sampling strategy.
 * 
 * @author Michael Scharfe, Jens Keilwagen
 */
public class SamplingHMMTrainingParameterSet extends HMMTrainingParameterSet {

	/**
	 * This is the empty constructor that can be used to fill the parameters after creation.
	 */
	public SamplingHMMTrainingParameterSet() {
		super();
		addParameters();
	}	

	/**
	 * This is the main constructor creating an already filled parameter set for training an {@link de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM} using a sampling strategy.
	 * 
	 * @param starts the number of independent starts
	 * @param stepsPerIteration the number of steps to be done for each start before testing for the end of the burn in phase
	 * @param stationarySteps the number of steps in the stationary phase
	 * @param burnInTestParameters the parameter set for the {@link de.jstacs.sampling.BurnInTest}
	 * 
	 * @throws IllegalValueException if some parameter is not permitted
	 */
	public SamplingHMMTrainingParameterSet( int starts, int stepsPerIteration, int stationarySteps, AbstractBurnInTestParameterSet burnInTestParameters ) throws IllegalValueException {
		super( starts );
		addParameters();
		getParameterAt( 1 ).setValue( stepsPerIteration );
		getParameterAt( 2 ).setValue( stationarySteps );
		parameters.get( 3 ).setValue( burnInTestParameters );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link SamplingHMMTrainingParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SamplingHMMTrainingParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public SamplingHMMTrainingParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	private void addParameters() {
		try{
			parameters.add( new SimpleParameter( DataType.INT,
					"stepsPerIteration",
					"the number of steps during one iteration",
					true,
					new NumberValidator<Integer>( 1, Integer.MAX_VALUE ) ) );

			parameters.add( new SimpleParameter( DataType.INT,
					"steps",
					"the total number of steps after burn-in",
					true,
					new NumberValidator<Integer>( 1, Integer.MAX_VALUE ) ) );	

			parameters.add(
					SubclassFinder.getSelectionParameter(
							AbstractBurnInTestParameterSet.class,
							AbstractBurnInTest.class.getPackage().getName(),
							"Burn in test parameters",
							"the parameters used to create a burn in test",
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
	 * This method returns the number of steps to be done in each start before testing for the end of the burn in phase (again).
	 * 
	 * @return the number of steps to be done in each start before testing for the end of the burn in phase (again)
	 */
	public int getNumberOfStepsPerIteration() {
		return (Integer)getParameterForName( "stepsPerIteration" ).getValue();
	}
	
	/**
	 * The method returns the number of steps to be done in the stationary phase.
	 * 
	 * @return the number of steps to be done in the stationary phase
	 */
	public int getNumberOfStepsInStationaryPhase() {
		return (Integer)getParameterForName( "steps" ).getValue();
	}

	/**
	 * This method return the burn in test to be used during sampling.
	 * 
	 * @return the burn in test to be used
	 * 
	 * @throws NotInstantiableException if the burn in test could not be created 
	 */
	public AbstractBurnInTest getBurnInTest() throws NotInstantiableException {
		return (AbstractBurnInTest)((AbstractBurnInTestParameterSet)getParameterForName( "Burn in test parameters" ).getValue()).getInstance();
	}
}