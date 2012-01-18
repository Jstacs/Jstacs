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

package de.jstacs.algorithms.optimization.termination;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.utils.Time;


/**
 * This class implements a {@link TerminationCondition} that stops the optimization if the elapsed time in seconds is
 * greater than a given value.
 * 
 * @author Jens Keilwagen
 *
 * @see Time
 */
public class TimeCondition extends AbstractTerminationCondition {

	private double seconds;
	
	/**
	 * This constructor creates an instance that does not allow any further iteration after <code>s</code> seconds
	 * 
	 * @param s the threshold of stopping the optimization (in seconds)
	 * 
	 * @throws Exception if the parameter can not be set properly
	 */
	public TimeCondition( double s ) throws Exception {
		this( new TimeConditionParameterSet( s ) );
	}
	
	/**
	 * This is the main constructor creating an instance from a given parameter set. 
	 * 
	 * @param parameter the set of parameters
	 * 
	 * @throws CloneNotSupportedException if <code>parameter</code> can not be cloned properly.
	 */
	public TimeCondition( TimeConditionParameterSet parameter ) throws CloneNotSupportedException {
		super( parameter );
	}
	
	protected void set() {
		this.seconds = (Double) parameter.getParameterAt(0).getValue();
	}
	
	protected String getXmlTag() {
		return "TimeCondition";
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link TimeCondition} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link TimeCondition} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public TimeCondition( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	public boolean doNextIteration( int iteration, double f_last, double f_current, double[] gradient, double[] direction, double alpha,
			Time t ) {
		return t.getElapsedTime() < seconds;
	}

	public boolean isSimple() {
		return true;
	}
	
	/**
	 * This class implements the parameter set for a {@link TimeCondition}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class TimeConditionParameterSet extends AbstractTerminationConditionParameterSet {

		/**
		 * This constructor creates an empty parameter set.
		 * @throws DatatypeNotValidException 
		 */
		public TimeConditionParameterSet() throws DatatypeNotValidException {
			super( TimeCondition.class );
			parameters.add( new SimpleParameter( DataType.DOUBLE,
					"seconds",
					"the number of seconds until stopping the algorithm",
					true,
					new NumberValidator<Double>( new Double( 0 ), new Double( Double.MAX_VALUE ) ) ) );
		}
		
		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Constructs an {@link TimeConditionParameterSet} out of an XML representation.
		 * 
		 * @param xml
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link TimeConditionParameterSet} could not be
		 *             reconstructed out of the {@link StringBuffer}
		 *             <code>representation</code>
		 */
		public TimeConditionParameterSet( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}
		
		/**
		 * This constructor creates a filled instance of a parameters set.
		 * 
		 * @param seconds the number of seconds until stopping the algorithm
		 * 
		 * @throws IllegalArgumentException if parameter can not be set
		 * @throws IllegalValueException if parameter can not be set
		 * @throws DatatypeNotValidException 
		 */
		public TimeConditionParameterSet( double seconds ) throws IllegalArgumentException, IllegalValueException, DatatypeNotValidException {
			this();
			this.getParameterAt( 0 ).setValue( seconds );
		}

		@Override
		public String getInstanceComment() {
			return "a set of parameters for the TimeCondition";
		}

		@Override
		public String getInstanceName() {
			return "TimeConditionParameterSet";
		}

	}
}
