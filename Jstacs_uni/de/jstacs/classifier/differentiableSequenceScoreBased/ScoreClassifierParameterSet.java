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

package de.jstacs.classifier.differentiableSequenceScoreBased;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition.AbstractTerminationConditionParameterSet;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition.SmallDifferenceOfFunctionEvaluationsConditionParameterSet;
import de.jstacs.classifier.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.AlphabetContainer.AlphabetContainerType;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.SequenceScoringParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.utils.SubclassFinder;

/**
 * A set of {@link de.jstacs.parameters.Parameter}s for any
 * {@link ScoreClassifier}.
 * 
 * @author Jens Keilwagen
 */
public abstract class ScoreClassifierParameterSet extends SequenceScoringParameterSet {

	private static final String[] algorithmStrings = new String[]{	"steepest descent",
																	"conjugate gradients (F., R.)",
																	"conjugate gradients (P., R. positive)",
																	"quasi newton (D., F., P.)",
																	"quasi newton (B., F., G., S.)",
																	"limited memory quasi newton (B., F., G., S.; n=3)",
																	"limited memory quasi newton (B., F., G., S.; n=4)",
																	"limited memory quasi newton (B., F., G., S.; n=5)",
																	"limited memory quasi newton (B., F., G., S.; n=6)",
																	"limited memory quasi newton (B., F., G., S.; n=7)",
																	"limited memory quasi newton (B., F., G., S.; n=8)",
																	"limited memory quasi newton (B., F., G., S.; n=9)",
																	"limited memory quasi newton (B., F., G., S.; n=10)" };

	private static final Byte[] algorithms = new Byte[]{ Optimizer.STEEPEST_DESCENT,
														Optimizer.CONJUGATE_GRADIENTS_FR,
														Optimizer.CONJUGATE_GRADIENTS_PRP,
														Optimizer.QUASI_NEWTON_DFP,
														Optimizer.QUASI_NEWTON_BFGS,
														(byte)3,
														(byte)4,
														(byte)5,
														(byte)6,
														(byte)7,
														(byte)8,
														(byte)9,
														(byte)10 };

	/**
	 * Creates a new {@link ScoreClassifierParameterSet} with empty parameter
	 * values.
	 * 
	 * @param instanceClass
	 *            the class of the instance
	 * @param simple
	 *            indicates whether the {@link AlphabetContainer} shall be
	 *            simple
	 * @param type
	 *            the type of the {@link AlphabetContainer}
	 * @param variableLength
	 *            indicates whether the corresponding classifier can handle
	 *            sequences of arbitrary length
	 * @throws Exception if the parameters could not be loaded
	 * 
	 * @see SequenceScoringParameterSet#SequenceScoringParameterSet(Class,
	 *      de.jstacs.data.AlphabetContainer.AlphabetContainerType, boolean, boolean)
	 */
	public ScoreClassifierParameterSet( Class<? extends ScoreClassifier> instanceClass, boolean simple, AlphabetContainerType type,
										boolean variableLength ) throws Exception {
		super( ScoreClassifier.class, type, simple, variableLength );
		addParameters();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ScoreClassifierParameterSet} out of its XML
	 * representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ScoreClassifierParameterSet} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see SequenceScoringParameterSet#SequenceScoringParameterSet(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public ScoreClassifierParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/**
	 * The constructor for a simple, instantiated parameter set.
	 * 
	 * @param instanceClass
	 *            the class of the instance
	 * @param alphabet
	 *            the alphabet
	 * @param length
	 *            the length of the sequences, 0 for homogeneous
	 * @param algo
	 *            the choice of algorithm
	 * @param eps
	 *            the threshold for stopping the algorithm
	 * @param lineps
	 *            the threshold for stopping the line search in the algorithm
	 * @param startD
	 *            the start distance for the line search in the algorithm
	 * @param free
	 *            indicates whether only the free parameters or all parameters
	 *            should be used
	 * @param kind
	 *            indicates the kind of class parameter initialization
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see SequenceScoringParameterSet#SequenceScoringParameterSet(Class, AlphabetContainer, int, boolean)
	 * @see #ScoreClassifierParameterSet(Class, AlphabetContainer, int, byte, AbstractTerminationCondition, double, double, boolean, de.jstacs.classifier.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter)
	 * @see SmallDifferenceOfFunctionEvaluationsCondition
	 */
	public ScoreClassifierParameterSet( Class<? extends ScoreClassifier> instanceClass, AlphabetContainer alphabet, int length, byte algo,
										double eps, double lineps, double startD, boolean free, KindOfParameter kind ) throws Exception {
		this( instanceClass, alphabet, length, algo, new SmallDifferenceOfFunctionEvaluationsCondition( eps ), lineps, startD, free, kind );
	}
	
	/**
	 * The constructor for a simple, instantiated parameter set.
	 * 
	 * @param instanceClass
	 *            the class of the instance
	 * @param alphabet
	 *            the alphabet
	 * @param length
	 *            the length of the sequences, 0 for homogeneous
	 * @param algo
	 *            the choice of algorithm
	 * @param tc
	 *            the termination condition for stopping the algorithm
	 * @param lineps
	 *            the threshold for stopping the line search in the algorithm
	 * @param startD
	 *            the start distance for the line search in the algorithm
	 * @param free
	 *            indicates whether only the free parameters or all parameters
	 *            should be used
	 * @param kind
	 *            indicates the kind of class parameter initialization
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see SequenceScoringParameterSet#SequenceScoringParameterSet(Class,
	 *      AlphabetContainer, int, boolean)
	 */
	public ScoreClassifierParameterSet( Class<? extends ScoreClassifier> instanceClass, AlphabetContainer alphabet, int length, byte algo,
										AbstractTerminationCondition tc, double lineps, double startD, boolean free, KindOfParameter kind ) throws Exception {
		super( instanceClass, alphabet, length, length == 0 );
		addParameters();
		parameters.get( 0 ).setValue( algorithmStrings[getIndex( algorithmStrings, algorithms, algo, false )] );
		parameters.get( 1 ).setValue( tc.getCurrentParameterSet() );
		parameters.get( 2 ).setValue( lineps );
		parameters.get( 3 ).setValue( startD );
		parameters.get( 4 ).setValue( free );
		parameters.get( 5 ).setValue( kind );
	}

	private void addParameters() throws Exception {
		parameters.add( new SelectionParameter( DataType.BYTE,
				algorithmStrings,
				algorithms,
				"algorithm",
				"the algorithm that should be used for numerical optimization",
				true ) );
		parameters.get(0).setDefault( algorithmStrings[4] );
		parameters.add(
				SubclassFinder.getSelectionParameter(
						AbstractTerminationConditionParameterSet.class,
						AbstractTerminationCondition.class.getPackage().getName(),
						"termination condition",
						"the terminantion condition for stopping the training algorithm",
						true
				)
		);
		parameters.get(1).setDefault( SmallDifferenceOfFunctionEvaluationsConditionParameterSet.class );
		parameters.add( new SimpleParameter( DataType.DOUBLE,
				"line epsilon",
				"the threshold for stopping the line search in the numerical training",
				true,
				new NumberValidator<Double>( 0d, Double.MAX_VALUE ),
				1E-9 ) );
		parameters.add( new SimpleParameter( DataType.DOUBLE,
				"start distance",
				"the start distance for the line search in the numerical training",
				true,
				new NumberValidator<Double>( 0d, Double.MAX_VALUE ),
				1d ) );
		parameters.add( new SimpleParameter( DataType.BOOLEAN,
				"free parameters",
				"Indicates whether only the free parameters or all parameters should be used.",
				true,
				new Boolean( false ) ) );
		parameters.add( new EnumParameter( KindOfParameter.class,
				"Indicates whether special plugIn parameters or the zero vector should be used as start parameters. For non-concave problems it is highly recommended to use plugIn parameters.",
				true, KindOfParameter.PLUGIN.name() ) );
	}

	/**
	 * This method indicates if only the free parameters shall be used.
	 * 
	 * @return <code>true</code> if only the free parameters shall be used,
	 *         <code>false</code> otherwise
	 */
	public boolean useOnlyFreeParameter() {
		return (Boolean)getParameterForName( "free parameters" ).getValue();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return getClass().getSimpleName();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
	 */
	@Override
	public String getInstanceComment() {
		return "holds the parameters for a score classifier";
	}
	
	/**
	 * This method returns the {@link AbstractTerminationCondition} for stopping the training, e.g., if the
	 * difference of the scores between two iterations is smaller than a given
	 * threshold the training is stopped.
	 * 
	 * @return the {@link AbstractTerminationCondition} for stopping the training
	 * @throws NotInstantiableException if the {@link AbstractTerminationCondition} could not be created from its {@link de.jstacs.parameters.ParameterSet}
	 */
	public AbstractTerminationCondition getTerminantionCondition() throws NotInstantiableException {
		return (AbstractTerminationCondition)(((InstanceParameterSet)getParameterForName("termination condition").getValue()).getInstance());
	}
}
