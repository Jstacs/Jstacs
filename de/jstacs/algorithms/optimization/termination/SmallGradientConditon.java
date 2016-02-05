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
package de.jstacs.algorithms.optimization.termination;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.utils.Time;


/**
 * This class implements a {@link TerminationCondition} that allows no further iteration in an optimization if the
 * the gradient becomes small, i.e.,
 * {@latex.inline $\\sum_i \\left|\\frac{\\partial f(\\underline{x})}{\\partial x_i}\\right| < \\epsilon$}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class SmallGradientConditon extends AbstractTerminationCondition {

	private double eps;
	
	/**
	 * This constructor creates an instance that stops the optimization if the sum of the absolute
	 * values of gradient components is smaller than <code>epsilon</code>.
	 *  
	 * @param epsilon the threshold for stopping, has to be non-negative
	 * 
	 * @throws Exception if the parameter can not be set properly
	 */
	public SmallGradientConditon( double epsilon ) throws Exception {
		this( new SmallGradientConditonParameterSet( epsilon ) );
	}
	
	/**
	 * This is the main constructor creating an instance from a given parameter set. 
	 * 
	 * @param parameter the set of parameters
	 * 
	 * @throws CloneNotSupportedException if <code>parameter</code> can not be cloned properly.
	 */
	public SmallGradientConditon( SmallGradientConditonParameterSet parameter ) throws CloneNotSupportedException {
		super( parameter );
	}
	
	protected void set() {
		this.eps = (Double) parameter.getParameterAt(0).getValue();
	}
	
	protected String getXmlTag() {
		return "SmallGradientConditon";
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link SmallGradientConditon} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SmallGradientConditon} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public SmallGradientConditon( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	public boolean doNextIteration( int iteration, double fLast, double fCurrent, double[] gradient, double[] direction, double alpha,
			Time t ) {
		double sum = 0;
		for(int i=0; i < gradient.length; i++){
			sum += Math.abs( gradient[i] );
		}
		//System.out.println("sum " + sum);
		return sum >= eps;
	}

	public boolean isSimple() {
		return false;
	}
	
	/**
	 * This class implements the parameter set for a {@link SmallStepCondition}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class SmallGradientConditonParameterSet extends AbstractTerminationConditionParameterSet {

		/**
		 * This constructor creates an empty parameter set.
		 */
		public SmallGradientConditonParameterSet() {
			super( SmallGradientConditon.class );
			try{
				parameters.add( new SimpleParameter( DataType.DOUBLE,
						"epsilon",
						"the epsilon for the gradient used for deciding whether to stop the algorithm or not",
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
		public SmallGradientConditonParameterSet( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}
		
		/**
		 * This constructor creates a filled instance of a parameters set.
		 * 
		 * @param eps the epsilon for the gradient
		 * 
		 * @throws IllegalArgumentException if parameter can not be set
		 * @throws IllegalValueException if parameter can not be set
		 */
		public SmallGradientConditonParameterSet( double eps ) throws IllegalArgumentException, IllegalValueException {
			this();
			getParameterAt( 0 ).setValue( eps );
		}

		@Override
		public String getInstanceComment() {
			return "a set of parameters for the SmallGradientConditon";
		}

		@Override
		public String getInstanceName() {
			return "SmallGradientConditonParameterSet";
		}
	}
}
