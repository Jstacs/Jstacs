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
 * This class stops an optimization if the value of the function
 * becomes less or equal to an absolute value, i.e.,
 * {@latex.inline $f(\\underline{x}_i) \\le t$}.
 * 
 * @author Jens Keilwagen
 * 
 * @deprecated use of the absolute value condition is not recommended and it may be removed in future releases
 */
public class AbsoluteValueCondition extends AbstractTerminationCondition {

	private double threshold;
	
	/**
	 * This constructor creates an instance that stops an minimization when the value of the function is below the given <code>threshold</code>  
	 * <br>
	 * 
	 * Be careful! If you set the value too low the minimization will not terminate.
	 * 
	 * @param threshold the threshold for stopping the optimization
	 * 
	 * @throws Exception if the parameter can not be set properly
	 */
	public AbsoluteValueCondition( double threshold ) throws Exception {
		this( new AbsoluteValueConditionParameterSet( threshold ) );
	}
	
	/**
	 * This is the main constructor creating an instance from a given parameter set. 
	 * 
	 * @param parameter the set of parameters
	 * 
	 * @throws CloneNotSupportedException if <code>parameter</code> can not be cloned properly.
	 */
	public AbsoluteValueCondition( AbsoluteValueConditionParameterSet parameter ) throws CloneNotSupportedException {
		super( parameter );
	}
		
	protected void set() {
		this.threshold = (Double) parameter.getParameterAt(0).getValue();
	}
	
	protected String getXmlTag() {
		return "AbsoluteValueCondition";
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link AbsoluteValueCondition} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AbsoluteValueCondition} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public AbsoluteValueCondition( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	public boolean doNextIteration( int iteration, double f_last, double f_current, double[] gradient, double[] direction, double alpha,
			Time t ) {
		return threshold < f_current;
	}

	public boolean isSimple() {
		return true;
	}
	
	/**
	 * This class implements the parameter set for a {@link AbsoluteValueCondition}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class AbsoluteValueConditionParameterSet extends AbstractTerminationConditionParameterSet {

		/**
		 * This constructor creates an empty parameter set.
		 */
		public AbsoluteValueConditionParameterSet() {
			super( AbsoluteValueCondition.class );
			try {
				parameters.add( new SimpleParameter( DataType.DOUBLE,
						"absolute value",
						"if the optimized value is smaller than this value the algorithm is stopped",
						true,
						new NumberValidator<Double>( new Double( 0 ), new Double( Double.MAX_VALUE ) ) ) );
			} catch( Exception e ) {
				//does not happen
				RuntimeException re = new RuntimeException(e.getMessage());
				re.setStackTrace(e.getStackTrace());
				throw re;
			}
		}
		
		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Constructs an {@link AbsoluteValueConditionParameterSet} out of an XML representation.
		 * 
		 * @param xml
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link AbsoluteValueConditionParameterSet} could not be
		 *             reconstructed out of the {@link StringBuffer}
		 *             <code>representation</code>
		 */
		public AbsoluteValueConditionParameterSet( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}
		
		/**
		 * This constructor creates a filled instance of a parameters set.
		 * 
		 * @param absValue the number of seconds until stopping the algorithm
		 * 
		 * @throws IllegalArgumentException if parameter can not be set
		 * @throws IllegalValueException if parameter can not be set
		 */
		public AbsoluteValueConditionParameterSet( double absValue ) throws IllegalArgumentException, IllegalValueException {
			this();
			this.getParameterAt( 0 ).setValue( absValue );
		}

		@Override
		public String getInstanceComment() {
			return "a set of parameters for the AbsoluteValueCondition";
		}

		@Override
		public String getInstanceName() {
			return "AbsoluteValueConditionParameterSet";
		}

	}
}
