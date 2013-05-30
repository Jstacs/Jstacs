package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous;

import de.jstacs.data.DataSet;
import de.jstacs.io.NonParsableException;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.ConstraintManager.Decomposition;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.DGTrainSMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.ConstraintParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.FSMEMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.IDGTrainSMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.MEManagerParameterSet;

/**
 * This class can be used for any discrete fixed structure maximum entropy model (FSMEM).
 * 
 * @author Jens Keilwagen
 */
public class FSMEManager extends MEManager
{
	/**
	 * Creates a new {@link MEManager} from a given
	 * {@link MEManagerParameterSet}.
	 * 
	 * @param params
	 *            the given parameter set
	 * 
	 * @throws CloneNotSupportedException
	 *             if the parameter set could not be cloned
	 * @throws IllegalArgumentException
	 *             if the parameter set is not instantiated
	 * @throws NonParsableException
	 *             if the parameter set is not parsable
	 * 
	 * @see MEManager#MEManager(MEManagerParameterSet)
	 */
	public FSMEManager( FSMEMParameterSet params ) throws CloneNotSupportedException, IllegalArgumentException,
			NonParsableException
	{
		super( params );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link FSMEManager} out of its XML representation.
	 * 
	 * @param stringBuff
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link FSMEManager} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see MEManager#MEManager(StringBuffer)
	 */
	public FSMEManager( StringBuffer stringBuff ) throws NonParsableException
	{
		super( stringBuff );
	}

	public String getInstanceName()
	{
		return "FSMEM(" + getConstraintString() + ")";
	}

	public void train( DataSet data, double[] weights ) throws Exception
	{
		//System.out.println( getConstraintString() );
		trainFactors( data, weights );
		trained = true;
	}

	private String getConstraintString()
	{
		Object constrParam = params.getParameterAt( 6 ).getValue();
		if( constrParam instanceof String )
		{
			return (String) constrParam;
		}
		else
		{
			return (String) ((ConstraintParameterSet) constrParam).getParameterAt( 0 ).getValue();
		}		
	}

	private static final String XML_TAG = "FSMEManager";

	@Override
	protected String getXMLTag()
	{
		return XML_TAG;
	}

	protected void set( DGTrainSMParameterSet params, boolean trained ) throws CloneNotSupportedException,
			NonParsableException
	{
		super.set( params, trained );
		if( !trained )
		{
			factors = getFactors( getConstraintString(), (Boolean) params.getParameterAt( 3 ).getValue(), (Decomposition) params
					.getParameterAt( 2 ).getValue() );
		}
	}
}
