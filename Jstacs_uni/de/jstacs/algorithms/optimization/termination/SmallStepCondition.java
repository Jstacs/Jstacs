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
import de.jstacs.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.utils.Time;

/**
 * This class implements a {@link TerminationCondition} that allows no further iteration in an optimization if the
 * scalar product of the current and the last values of <code>x</code> will be small, i.e.,
 * {@latex.inline $(\\underline{x}_i-\\underline{x}_{i-1})^T (\\underline{x}_i-\\underline{x}_{i-1}) < \\epsilon$}.
 * 
 * @author Jens Keilwagen
 */
public class SmallStepCondition extends AbstractTerminationCondition {

	private double epsilon;
	
	/**
	 * This constructor creates an instance that allows no further iteration in an optimization if the
	 * scalar product of the current and the last values of <code>x</code> is smaller than <code>epsilon</code>.
	 * 
	 * @param epsilon the threshold for stopping, has to be non-negative
	 * 
	 * @throws Exception if the parameter can not be set properly
	 */
	public SmallStepCondition( double epsilon ) throws Exception {
		this( new SmallStepConditionParameterSet( epsilon ) );
	}
	
	/**
	 * This is the main constructor creating an instance from a given parameter set. 
	 * 
	 * @param parameter the set of parameters
	 * 
	 * @throws CloneNotSupportedException if <code>parameter</code> can not be cloned properly.
	 */
	public SmallStepCondition( SmallStepConditionParameterSet parameter ) throws CloneNotSupportedException {
		super( parameter );
	}
	
	protected void set() {
		this.epsilon = (Double) parameter.getParameterAt(0).getValue();
	}
	
	protected String getXmlTag() {
		return "SmallStepCondition";
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link TimeCondition} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SmallStepCondition} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public SmallStepCondition( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	public boolean doNextIteration( int iteration, double f_last, double f_current, double[] gradient, double[] direction, double alpha,
			Time t ) {
		double s = 0;
		for( int counter = 0; counter < direction.length; counter++ ) {
			s += direction[counter] * direction[counter];
		}
		return s * ( alpha * alpha ) >= this.epsilon;
	}
	
	public boolean isSimple() {
		return false;
	}
	
	/**
	 * This class implements the parameter set for a {@link SmallStepCondition}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class SmallStepConditionParameterSet extends AbstractTerminationConditionParameterSet {

		/**
		 * This constructor creates an empty parameter set.
		 */
		public SmallStepConditionParameterSet() {
			super( SmallStepCondition.class );
		}
		
		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Constructs an {@link SmallStepCondition} out of an XML representation.
		 * 
		 * @param xml
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link SmallStepCondition} could not be
		 *             reconstructed out of the {@link StringBuffer}
		 *             <code>representation</code>
		 */
		public SmallStepConditionParameterSet( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}
		
		/**
		 * This constructor creates a filled instance of a parameters set.
		 * 
		 * @param eps the epsilon for the size of the step
		 * 
		 * @throws IllegalArgumentException if parameter can not be set
		 * @throws IllegalValueException if parameter can not be set
		 */
		public SmallStepConditionParameterSet( double eps ) throws IllegalArgumentException, IllegalValueException {
			this();
			this.getParameterAt(0).setValue( eps );
		}

		@Override
		public String getInstanceComment() {
			return "a set of parameters for the SmallStepCondition";
		}

		@Override
		public String getInstanceName() {
			return "SmallStepConditionParameterSet";
		}

		@Override
		protected void loadParameters() throws Exception {
			initParameterList();
			parameters.add( new SimpleParameter( DataType.DOUBLE,
					"epsilon",
					"the epsilon for the size of the step used for deciding whether to stop the algorithm or not",
					true,
					new NumberValidator<Double>( new Double( 0 ), new Double( Double.MAX_VALUE ) ) ) );

		}
	}
}
