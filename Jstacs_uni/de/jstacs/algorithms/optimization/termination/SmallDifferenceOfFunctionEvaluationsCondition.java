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
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.utils.Time;

/**
 * This class implements a {@link TerminationCondition} that stops an optimization
 * if the difference of the current and the last function evaluations will be small, i.e., 
 * {@latex.inline $|f(\\underline{x}_{i-1}) - f(\\underline{x}_i)| < \\epsilon$}.
 * 
 * @author Jens Keilwagen
 */
public class SmallDifferenceOfFunctionEvaluationsCondition extends AbstractTerminationCondition {

	private double eps;
	
	/**
	 * This constructor creates an instance that stops the optimization if the difference of the
	 * current and the last function evaluations is smaller than <code>epsilon</code>.
	 *  
	 * @param epsilon the threshold for stopping, has to be non-negative
	 * 
	 * @throws Exception if the parameter can not be set properly
	 */
	public SmallDifferenceOfFunctionEvaluationsCondition( double epsilon ) throws Exception {
		this( new SmallDifferenceOfFunctionEvaluationsConditionParameterSet( epsilon ) );
	}
	
	/**
	 * This is the main constructor creating an instance from a given parameter set. 
	 * 
	 * @param parameter the set of parameters
	 * 
	 * @throws CloneNotSupportedException if <code>parameter</code> can not be cloned properly.
	 */
	public SmallDifferenceOfFunctionEvaluationsCondition( SmallDifferenceOfFunctionEvaluationsConditionParameterSet parameter ) throws CloneNotSupportedException {
		super( parameter );
	}
	
	protected void set() {
		this.eps = (Double) parameter.getParameterAt(0).getValue();
	}
	
	protected String getXmlTag() {
		return "SmallDifferenceOfFunctionEvaluationsCondition";
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link SmallDifferenceOfFunctionEvaluationsCondition} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SmallDifferenceOfFunctionEvaluationsCondition} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public SmallDifferenceOfFunctionEvaluationsCondition( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
		
	public boolean doNextIteration( int iteration, double f_last, double f_current, double[] gradient, double[] direction, double alpha,
			Time t ) {
		return eps <= Math.abs(f_last-f_current);
	}
	
	public boolean isSimple() {
		return true;
	}
	
	/**
	 * This class implements the parameter set for a {@link SmallDifferenceOfFunctionEvaluationsCondition}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class SmallDifferenceOfFunctionEvaluationsConditionParameterSet extends AbstractTerminationConditionParameterSet {

		/**
		 * This constructor creates an empty parameter set.
		 */
		public SmallDifferenceOfFunctionEvaluationsConditionParameterSet() {
			super( SmallDifferenceOfFunctionEvaluationsCondition.class );
			try {
				parameters.add( new SimpleParameter( DataType.DOUBLE,
						"epsilon",
						"the epsilon for the difference of function evaluations used for deciding whether to stop the algorithm or not",
						true,
						new NumberValidator<Double>( 0d, Double.MAX_VALUE ),
						1E-6 ) );
			} catch( ParameterException doesNotHappen ) {
				throw new RuntimeException( doesNotHappen.getMessage() );
			}
		}
		
		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Constructs an {@link SmallDifferenceOfFunctionEvaluationsCondition} out of an XML representation.
		 * 
		 * @param xml
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link SmallDifferenceOfFunctionEvaluationsCondition} could not be
		 *             reconstructed out of the {@link StringBuffer}
		 *             <code>representation</code>
		 */
		public SmallDifferenceOfFunctionEvaluationsConditionParameterSet( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}
		
		/**
		 * This constructor creates a filled instance of a parameters set.
		 * 
		 * @param eps the epsilon for the difference of function evaluations
		 * 
		 * @throws ParameterException if parameter can not be set
		 */
		public SmallDifferenceOfFunctionEvaluationsConditionParameterSet( double eps ) throws ParameterException {
			this();
			this.getParameterAt( 0 ).setValue( eps );
		}

		@Override
		public String getInstanceComment() {
			return "a set of parameters for the SmallDifferenceOfFunctionEvaluationsCondition";
		}

		@Override
		public String getInstanceName() {
			return "SmallDifferenceOfFunctionEvaluationsConditionParameterSet";
		}
	}
}
