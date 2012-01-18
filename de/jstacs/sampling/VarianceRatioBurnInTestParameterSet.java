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
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * Class for the parameters of a {@link VarianceRatioBurnInTest}. This
 * class fulfills the requirements of a {@link AbstractBurnInTestParameterSet} and
 * can be used to create a new {@link VarianceRatioBurnInTest}.
 * 
 * @author Jens Keilwagen
 */
public class VarianceRatioBurnInTestParameterSet extends AbstractBurnInTestParameterSet {

	/**
	 * Creates a new {@link VarianceRatioBurnInTestParameterSet} with empty parameter values.
	 * 
	 * @throws IllegalArgumentException if <code>instanceClass</code> is <code>null</code> 
	 */
	public VarianceRatioBurnInTestParameterSet() throws IllegalArgumentException {
		super( VarianceRatioBurnInTest.class );
		try {
			parameters.add( new SimpleParameter( DataType.DOUBLE, "threshold", "the threshold value for testing the end of the burn-in phase" +
					"with the Variance-Ratio burn-in test, the value has to be greater than 1 since the tested potential scale reduction" +
					"factor R converges to 1", true, new NumberValidator<Double>( 1d, Double.MAX_VALUE ), 1.2 ) );
		} catch ( ParameterException doesNotHappen ) {
			throw new RuntimeException( doesNotHappen );
		}	
	}
	
	/**
	 * Creates a new {@link VarianceRatioBurnInTestParameterSet} with
	 * pre-defined parameter values.
	 * 
	 * @param starts
	 *            the number of runs the Gibbs Sampler will be started
	 * @param t
	 *            the threshold for determining the end of the burn-in phase
	 *            
	 * @throws IllegalArgumentException if <code>instanceClass</code> is <code>null</code>
	 * @throws IllegalValueException if <code>t</code> can not be set
	 */
	public VarianceRatioBurnInTestParameterSet( int starts, double t ) throws IllegalArgumentException, IllegalValueException {
		super( VarianceRatioBurnInTest.class, starts );
		parameters.get( 1 ).setValue( t );
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link VarianceRatioBurnInTestParameterSet} out of an XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link VarianceRatioBurnInTestParameterSet} could not be
	 *             reconstructed out of the {@link StringBuffer}
	 *             <code>representation</code>
	 */
	public VarianceRatioBurnInTestParameterSet( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}
	
	/**
	 * Returns the threshold used in the {@link VarianceRatioBurnInTestParameterSet}.
	 * @return the threshold used in the {@link VarianceRatioBurnInTestParameterSet}
	 */
	public double getThreshold() {
		return (Double) getParameterAt( 1 ).getValue();
	}

	public String getInstanceComment() {
		return "The parameter set for a " + VarianceRatioBurnInTest.class.getName() + ".";
	}

	public String getInstanceName() {
		return getClass().getSimpleName();
	}
}
