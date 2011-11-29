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

package de.jstacs.models.hmm;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * This class implements an abstract {@link ParameterSet} that is used for the training of an {@link AbstractHMM}.
 * 
 * @author Jens Keilwagen
 */
public class HMMTrainingParameterSet extends ParameterSet {

	/**
	 * This is the empty constructor that can be used to fill the parameters after creation.
	 */
	protected HMMTrainingParameterSet() {
	}
	
	/**
	 * This constructor can be used to create an instance with a specified number of starts.
	 * 
	 * @param starts the number of different starts
	 * 
	 * @throws IllegalValueException if the parameter can not be set
	 */
	protected HMMTrainingParameterSet( int starts ) throws IllegalValueException {
		getParameterAt( 0 ).setValue( starts );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link HMMTrainingParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link HMMTrainingParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	protected HMMTrainingParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	@Override
	protected void loadParameters() throws Exception {
		parameters = new ParameterList();
		parameters.add( new SimpleParameter( DataType.INT, "starts", "the number of different starts of HMM training", true, new NumberValidator<Integer>(1,Integer.MAX_VALUE)) );
	}

	/**
	 * The method returns the number of different starts.
	 * 
	 * @return the number of starts
	 */
	public int getNumberOfStarts() {
		return (Integer) getParameterAt( 0 ).getValue();
	}
}