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
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.utils.Time;

/**
 * This class will stop an optimization if the number of iteration reaches  a given number.
 * 
 * @author Jens Keilwagen
 */
public class IterationCondition extends AbstractTerminationCondition {

	private int maxIter;
	
	/**
	 * This constructor creates an instance that does not allow any further iteration after <code>maxIter</code> iterations.
	 * 
	 * @param maxIter the maximal number of iterations
	 *
	 * @throws Exception if the parameter can not be set properly
	 */
	public IterationCondition( int maxIter ) throws Exception {
		this( new IterationConditionParameterSet( maxIter ) );
	}
	
	/**
	 * This is the main constructor creating an instance from a given parameter set. 
	 * 
	 * @param parameter the set of parameters
	 * 
	 * @throws CloneNotSupportedException if <code>parameter</code> can not be cloned properly.
	 */
	public IterationCondition( IterationConditionParameterSet parameter ) throws CloneNotSupportedException {
		super( parameter );
	}
	
	protected void set() {
		this.maxIter = (Integer) parameter.getParameterAt(0).getValue();
	}
	
	protected String getXmlTag() {
		return "IterationCondition";
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link IterationCondition} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link IterationCondition} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public IterationCondition( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
		
	public boolean doNextIteration( int iteration, double f_last, double f_current, double[] gradient, double[] direction, double alpha,
			Time t ) {
		return iteration < maxIter;
	}
	
	public boolean isSimple() {
		return true;
	}
	
	/**
	 * This class implements the parameter set for a {@link IterationCondition}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class IterationConditionParameterSet extends AbstractTerminationConditionParameterSet {

		/**
		 * This constructor creates an empty parameter set.
		 * @throws DatatypeNotValidException 
		 */
		public IterationConditionParameterSet() throws DatatypeNotValidException {
			super( IterationCondition.class );
			parameters.add( new SimpleParameter( DataType.INT,
					"maximal iteration",
					"the maximal number of iterations for stopping an algorithm",
					true,
					new NumberValidator<Integer>( new Integer( 0 ), new Integer( Integer.MAX_VALUE ) ) ) );
		}
		
		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Constructs an {@link IterationConditionParameterSet} out of an XML representation.
		 * 
		 * @param xml
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link IterationConditionParameterSet} could not be
		 *             reconstructed out of the {@link StringBuffer}
		 *             <code>representation</code>
		 */
		public IterationConditionParameterSet( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}
		
		/**
		 * This constructor creates a filled instance of a parameters set.
		 * 
		 * @param maxIter the maximal number of iterations of the algorithm
		 * 
		 * @throws IllegalArgumentException if parameter can not be set
		 * @throws IllegalValueException if parameter can not be set
		 * @throws DatatypeNotValidException 
		 */
		public IterationConditionParameterSet( int maxIter) throws IllegalArgumentException, IllegalValueException, DatatypeNotValidException {
			this();
			this.getParameterAt( 0 ).setValue( maxIter );
		}

		@Override
		public String getInstanceComment() {
			return "a set of parameters for the IterationCondition";
		}

		@Override
		public String getInstanceName() {
			return "IterationConditionParameterSet";
		}
	}
}
