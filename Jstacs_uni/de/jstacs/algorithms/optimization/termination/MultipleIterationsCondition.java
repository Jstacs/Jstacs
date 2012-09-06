package de.jstacs.algorithms.optimization.termination;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.utils.SubclassFinder;
import de.jstacs.utils.Time;

/**
 * This {@link TerminationCondition} requires another provided {@link TerminationCondition} to fail a contiguous specified number of times
 * before the optimization is terminated.
 * 
 * @author Jan Grau
 *
 */
public class MultipleIterationsCondition extends AbstractTerminationCondition {

	private int thresh;
	private int number;
	private AbstractTerminationCondition test;
	
	/**
	 * This constructor creates an instance that stops the optimization if provided termination condition
	 * fails a contiguously a specified number of times.
	 *  
	 * @param threshold the number of contiguous iterations the provided condition must fail to stop the optimization
	 * @param condition the condition for stopping the algorithm
	 * 
	 * @throws Exception if the parameter can not be set properly
	 */
	public MultipleIterationsCondition( int threshold, AbstractTerminationCondition condition ) throws Exception {
		this( new MultipleIterationsConditionParameterSet( threshold, condition ) );
	}
	
	/**
	 * This is the main constructor creating an instance from a given parameter set. 
	 * 
	 * @param parameter the set of parameters
	 * 
	 * @throws CloneNotSupportedException if <code>parameter</code> can not be cloned properly.
	 */
	public MultipleIterationsCondition( MultipleIterationsConditionParameterSet parameter ) throws CloneNotSupportedException {
		super( parameter );
	}
	
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link MultipleIterationsCondition} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MultipleIterationsCondition} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public MultipleIterationsCondition( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	
	@Override
	public boolean doNextIteration( int iteration, double f_last, double f_current, double[] gradient, double[] direction, double alpha,
			Time t ) {
		if(test.doNextIteration( iteration, f_last, f_current, gradient, direction, alpha, t )){
			number = 0;
		}else{
			number++;
		}
		if(number>=thresh){
			return false;
		}else{
			return true;
		}
	}

	@Override
	public boolean isSimple() {
		return test.isSimple();
	}

	@Override
	protected String getXmlTag() {
		return "MultipleIterationsCondition";
	}

	@Override
	protected void set() {
		try{
			this.thresh = (Integer) parameter.getParameterAt(0).getValue();
			this.test = ((AbstractTerminationConditionParameterSet)parameter.getParameterAt( 1 ).getValue()).getInstance();
		} catch ( NotInstantiableException e ) {
			RuntimeException ex = new RuntimeException( e.getMessage() );
			ex.setStackTrace( e.getStackTrace());
			throw ex;
		}
	}

	
	/**
	 * This class implements the parameter set for a {@link MultipleIterationsCondition}.
	 * 
	 * @author Jan Grau, Jens Keilwagen
	 */
	public static class MultipleIterationsConditionParameterSet extends AbstractTerminationConditionParameterSet {

		/**
		 * This constructor creates an empty parameter set.
		 * 
		 * @throws Exception if some error occurs
		 */
		public MultipleIterationsConditionParameterSet() throws Exception {
			super( MultipleIterationsCondition.class );
			parameters.add( new SimpleParameter( DataType.INT,
					"threshold",
					"the number of iterations the provided condition must fail to stop the optimization",
					true,
					new NumberValidator<Integer>(1,Integer.MAX_VALUE)
				)
			);
			parameters.add( SubclassFinder.getSelectionParameter( AbstractTerminationConditionParameterSet.class, AbstractTerminationConditionParameterSet.class.getPackage().getName(), "Termination condition", "Select a termination condition.", true ));
		}
		
		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Constructs an {@link MultipleIterationsConditionParameterSet} out of an XML representation.
		 * 
		 * @param xml
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link MultipleIterationsConditionParameterSet} could not be
		 *             reconstructed out of the {@link StringBuffer}
		 *             <code>representation</code>
		 */
		public MultipleIterationsConditionParameterSet( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}
		
		/**
		 * This constructor creates a filled instance of the parameter set.
		 * 
		 * @param threshold the number of contiguous iterations the provided condition must fail to stop the optimization
		 * @param condition the condition for stopping the algorithm
		 * @throws Exception if the {@link MultipleIterationsConditionParameterSet} could not be created
		 */
		public MultipleIterationsConditionParameterSet( int threshold, AbstractTerminationCondition condition ) throws Exception {
			this();
			this.getParameterAt( 0 ).setValue( threshold );
			this.getParameterAt( 1 ).setValue( condition.parameter );
		}

		@Override
		public String getInstanceComment() {
			return "a set of parameters for the MultipleIterationsCondition";
		}

		@Override
		public String getInstanceName() {
			return "MultipleIterationsConditionParameterSet";
		}

	}
	
}
