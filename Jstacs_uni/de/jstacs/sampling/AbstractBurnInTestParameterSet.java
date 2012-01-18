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

package de.jstacs.sampling;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * Class for the parameters of a {@link AbstractBurnInTest}. This
 * class fulfills the requirements of a {@link InstanceParameterSet} and
 * can be used to create a new {@link AbstractBurnInTest}.
 * 
 * @author Jens Keilwagen
 */
public abstract class AbstractBurnInTestParameterSet extends InstanceParameterSet<AbstractBurnInTest> {

	/**
	 * Creates a new {@link AbstractBurnInTestParameterSet} with empty parameter values.
	 * 
	 * @param instanceClass
	 *            the class to be instantiated
	 *            
	 * @throws IllegalArgumentException if <code>instanceClass</code> is <code>null</code>
	 */
	protected AbstractBurnInTestParameterSet( Class<? extends AbstractBurnInTest> instanceClass ) throws IllegalArgumentException {
		super( instanceClass );
		try {
			parameters.add( new SimpleParameter( DataType.INT, "starts", "the number of Gibbs Sampling starts", true, new NumberValidator<Integer>( 3, Integer.MAX_VALUE ) ) );
		} catch (DatatypeNotValidException doesNotHappen) {
			throw new RuntimeException( doesNotHappen );
		}
	}
	
	/**
	 * Creates a new {@link AbstractBurnInTestParameterSet} with
	 * pre-defined parameter values.
	 * 
	 * @param instanceClass
	 *            the class to be instantiated
	 * @param starts
	 *            the number of runs the Gibbs Sampler will be started
	 *            
	 * @throws IllegalArgumentException if <code>instanceClass</code> is <code>null</code>
	 * @throws IllegalValueException if <code>starts</code> can not be set
	 */
	protected AbstractBurnInTestParameterSet( Class<? extends AbstractBurnInTest> instanceClass, int starts ) throws IllegalArgumentException, IllegalValueException {
		this(instanceClass);
		parameters.get( 0 ).setValue( starts );
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link AbstractBurnInTestParameterSet} out of an XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AbstractBurnInTestParameterSet} could not be
	 *             reconstructed out of the {@link StringBuffer}
	 *             <code>representation</code>
	 */
	protected AbstractBurnInTestParameterSet( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}
	
	/**
	 * Returns the number of starts.
	 * 
	 * @return the number of starts
	 */
	public int getNumberOfStarts() {
		return (Integer) getParameterAt( 0 ).getValue();
	}
	
	/* 
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.ParameterSet#clone()
	 */
	public AbstractBurnInTestParameterSet clone() throws CloneNotSupportedException {
		return (AbstractBurnInTestParameterSet)super.clone();
	}
}
