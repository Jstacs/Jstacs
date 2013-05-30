package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters;

import de.jstacs.DataType;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.ConstraintManager.Decomposition;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.MEMTools;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.MEManager;

/**
 * The ParameterSet for any MEManager.
 * 
 * @author Jens Keilwagen
 */
public abstract class MEManagerParameterSet extends IDGTrainSMParameterSet
{
	private static final String[] algorithmStrings = new String[]{ "BGIS (p-space)", "SGIS (p-space)", "GIS", "BGIS",
			"SGIS", "steepest descent", "conjugate gradients (F., R.)", "conjugate gradients (P., R. positive)",
			"quasi newton (D., F., P.)", "quasi newton (B., F., G., S.)",
			"limited memory quasi newton (B., F., G., S.; n=3)", "limited memory quasi newton (B., F., G., S.; n=4)",
			"limited memory quasi newton (B., F., G., S.; n=5)", "limited memory quasi newton (B., F., G., S.; n=6)",
			"limited memory quasi newton (B., F., G., S.; n=7)", "limited memory quasi newton (B., F., G., S.; n=8)",
			"limited memory quasi newton (B., F., G., S.; n=9)", "limited memory quasi newton (B., F., G., S.; n=10)" };

	private static final Byte[] algorithms = new Byte[]{ new Byte( MEMTools.BGIS_P ), new Byte( MEMTools.SGIS_P ),
			new Byte( MEMTools.GIS ), new Byte( MEMTools.BGIS ), new Byte( MEMTools.SGIS ), new Byte( Optimizer.STEEPEST_DESCENT ),
			new Byte( Optimizer.CONJUGATE_GRADIENTS_FR ), new Byte( Optimizer.CONJUGATE_GRADIENTS_PRP ),
			new Byte( Optimizer.QUASI_NEWTON_DFP ), new Byte( Optimizer.QUASI_NEWTON_BFGS ), new Byte( (byte) 3 ),
			new Byte( (byte) 4 ), new Byte( (byte) 5 ), new Byte( (byte) 6 ), new Byte( (byte) 7 ),
			new Byte( (byte) 8 ), new Byte( (byte) 9 ), new Byte( (byte) 10 ) };

	/**
	 * The constructor for the {@link de.jstacs.Storable} interface.
	 * 
	 * @param s
	 *            the StringBuffer
	 * 
	 * @throws NonParsableException
	 *             if the StringBuffer is not parsable
	 */
	public MEManagerParameterSet( StringBuffer s ) throws NonParsableException
	{
		super( s );
	}

	/**
	 * The constructor an empty constructor of extended class.
	 * 
	 * @param instanceClass
	 *            the (sub-)class
	 *            
	 * @throws Exception if some parameters template could not be created properly
	 */
	public MEManagerParameterSet( Class<? extends MEManager> instanceClass ) throws Exception
	{
		super( instanceClass );
		addParameters();
	}

	/**
	 * The fast constructor.
	 * 
	 * @param instanceClass
	 *            the (sub-)class
	 * @param alphabet
	 *            the alphabet
	 * @param length
	 *            the length of the modeled sequences
	 * @param ess
	 *            the ess
	 * @param description
	 *            the description
	 * @param decomposition
	 *            the kind of decomposition
	 * @param reduce
	 *            whether the constraints should be reduced
	 * @param algorithm
	 *            the choice of algorithm
	 * @param epsilon
	 *            the threshold for stopping the numerical algorithms
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see de.jstacs.models.discrete.ConstraintManager.Decomposition#DECOMPOSE_NOTHING
	 * @see de.jstacs.models.discrete.ConstraintManager.Decomposition#DECOMPOSE_UNCONNECTED
	 * @see de.jstacs.models.discrete.ConstraintManager.Decomposition#DECOMPOSE_LESS_CONNECTED
	 */
	public MEManagerParameterSet( Class<? extends MEManager> instanceClass, AlphabetContainer alphabet, int length,
			double ess, String description, Decomposition decomposition, boolean reduce, byte algorithm, double epsilon )
			throws Exception
	{
		super( instanceClass, alphabet, length, ess, description );
		addParameters();
		parameters.get( 2 ).setValue( decomposition );
		parameters.get( 3 ).setValue( new Boolean( reduce ) );
		parameters.get( 4 ).setValue(
				algorithmStrings[getIndex( algorithmStrings, algorithms, new Byte( algorithm ), false )] );
		parameters.get( 5 ).setValue( new Double( epsilon ) );
	}

	/**
	 * Adds the parameter template in the constructor.
	 * 
	 * @throws Exception if some parameters could not be created properly
	 */
	protected void addParameters() throws Exception
	{
		parameters.add( new EnumParameter( Decomposition.class, "the kind the model should be decomposed", true ) );
		parameters.get( 2 ).setDefault( Decomposition.DECOMPOSE_LESS_CONNECTED );
		parameters.add( new SimpleParameter( DataType.BOOLEAN, "reduce",
				"whether the constraints should be reduced or not", true, new Boolean( true ) ) );
		parameters.add( new SelectionParameter( DataType.BYTE, algorithmStrings, algorithms, "algorithm",
				"the algorithm that should be used for numerical optimization", true ) );
		parameters.get( 4 ).setDefault( "BGIS (p-space)" );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "epsilon",
				"the bound for stopping the numercal optimization algorithm", true, new NumberValidator<Double>(
						new Double( 0 ), new Double( Double.MAX_VALUE ) ), new Double( 1E-6 ) ) );
	}
}
