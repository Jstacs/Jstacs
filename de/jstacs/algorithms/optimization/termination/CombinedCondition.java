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

import java.io.IOException;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.utils.SubclassFinder;
import de.jstacs.utils.Time;

/**
 * This class allows to use many {@link TerminationCondition}s at once.
 * 
 * @author Jens Keilwagen
 */
public class CombinedCondition extends AbstractTerminationCondition {

	private TerminationCondition[] condition;
	private int threshold;
	
	/**
	 * This constructor creates an instance that allows to use many {@link TerminationCondition}s at once. 
	 * 
	 * @param threshold the number of conditions that has to be fulfilled;
	 * 			if all conditions should be fulfilled than <code>threshold=condition.length</code> (equal to AND);
	 * 			if at least one condition should be fulfilled than <code>threshold=1</code> (equal to OR);
	 * @param condition the conditions that are used to create this instance and that are used to determine whether another iteration should be done
	 * 
	 * @throws Exception if the parameter can not be set properly
	 */
	public CombinedCondition( int threshold, AbstractTerminationCondition...  condition ) throws Exception {
		this( new CombinedConditionParameterSet( threshold, condition ) );
	}
	
	/**
	 * This is the main constructor creating an instance from a given parameter set. 
	 * 
	 * @param parameter the set of parameters
	 * 
	 * @throws CloneNotSupportedException if <code>parameter</code> can not be cloned properly.
	 */
	public CombinedCondition( CombinedConditionParameterSet parameter ) throws CloneNotSupportedException {
		super( parameter );
	}
	
	@Override
	protected void set() {
		this.threshold = (Integer) parameter.getParameterAt(0).getValue();
		ExpandableParameterSet eps = (ExpandableParameterSet) parameter.getParameterAt( 1 ).getValue();
		condition = new AbstractTerminationCondition[eps.getNumberOfParameters()];
		for(int i=0;i<condition.length;i++){
			Parameter p = ((ParameterSet) eps.getParameterAt( i ).getValue()).getParameterAt(0);
			try {
				condition[i] = (AbstractTerminationCondition) ((AbstractTerminationConditionParameterSet)p.getValue()).getInstance();
			} catch ( NotInstantiableException e ) {
				RuntimeException ex = new RuntimeException( e.getMessage() );
				ex.setStackTrace( e.getStackTrace());
				throw ex;
			}
		}
	}
	
	protected String getXmlTag() {
		return "CombinedCondition";
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link CombinedCondition} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link CombinedCondition} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public CombinedCondition( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
		
	public boolean doNextIteration( int iteration, double f_last, double f_current, double[] gradient, double[] direction, double alpha,
			Time t ) {
		int positive = 0;
		for( int i = 0; i < condition.length; i++ ) {
			if( condition[i].doNextIteration( iteration, f_last, f_current, gradient, direction, alpha, t ) ){
				positive++;
			}
		}
		return positive >= threshold;		
	}
	
	public boolean isSimple() {
		int i = 0;
		while( i < condition.length && condition[i].isSimple() ) {
			i++;
		}
		return i == condition.length;
	}
	
	/**
	 * This class implements the parameter set for a {@link CombinedCondition}.
	 * 
	 * @author Jan Grau, Jens Keilwagen
	 */
	public static class CombinedConditionParameterSet extends AbstractTerminationConditionParameterSet {

		/**
		 * This constructor creates an empty parameter set.
		 * @throws DatatypeNotValidException 
		 * @throws IOException 
		 * @throws ClassNotFoundException 
		 * @throws IllegalAccessException 
		 * @throws InstantiationException 
		 * @throws CloneNotSupportedException 
		 */
		public CombinedConditionParameterSet() throws DatatypeNotValidException, CloneNotSupportedException, InstantiationException, IllegalAccessException, ClassNotFoundException, IOException {
			super( CombinedCondition.class );
			parameters.add( new SimpleParameter( DataType.INT,
					"threshold",
					"the number of conditions that has to be fulfilled for stopping the algorithm",
					true
				)
			);
			parameters.add( new ParameterSetContainer( "Termination conditions", "The set of termination conditions that shall be combined.", 
					new ExpandableParameterSet( new SimpleParameterSet(
							SubclassFinder.getCollection( AbstractTerminationCondition.class, "de.jstacs", "Termination condition", "Select a termination condition.", true ) )
							, "Termination conditions", "Add termination conditions to the set of conditions." )
			) );
		}
		
		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Constructs an {@link CombinedConditionParameterSet} out of an XML representation.
		 * 
		 * @param xml
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link CombinedConditionParameterSet} could not be
		 *             reconstructed out of the {@link StringBuffer}
		 *             <code>representation</code>
		 */
		public CombinedConditionParameterSet( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}
		
		/**
		 * This constructor creates a filled instance of a parameters set.
		 * 
		 * @param threshold the number of conditions that has to be fulfilled
		 * @param condition the individual conditions for stopping the algorithm
		 * @throws Exception if the {@link CombinedConditionParameterSet} could not be created
		 */
		public CombinedConditionParameterSet( int threshold, AbstractTerminationCondition[] condition ) throws Exception {
			this();
			this.getParameterAt( 0 ).setValue( threshold );
			ExpandableParameterSet eps = (ExpandableParameterSet) this.getParameterAt( 1 ).getValue();
			for(int i=0;i<condition.length;i++){
				if( i != 0 ) {
					eps.addParameterToSet();
				}
				((ParameterSet) eps.getParameterAt( i ).getValue()).getParameterAt(0).setValue( condition[i].parameter );				
			}
		}

		@Override
		public String getInstanceComment() {
			return "a set of parameters for the CombinedCondition";
		}

		@Override
		public String getInstanceName() {
			return "CombinedConditionParameterSet";
		}

	}
}
