package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous;

import java.text.NumberFormat;
import java.util.ArrayList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.NotTrainedException;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.ConstraintManager;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.ConstraintManager.Decomposition;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.MEManagerParameterSet;
import de.jstacs.utils.DiscreteInhomogenousDataSetEmitter;

/**
 * This class is the super class for all maximum entropy models
 * 
 * @author Jens Keilwagen
 */
public abstract class MEManager extends InhomogeneousDGTrainSM
{
	/**
	 * The independent maximum entropy models.
	 */
	protected MEM[] factors;

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
	 * @see InhomogeneousDGTrainSM#InhomogeneousDGTrainSM(de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.IDGTrainSMParameterSet)
	 */
	public MEManager( MEManagerParameterSet params ) throws CloneNotSupportedException, IllegalArgumentException,
			NonParsableException
	{
		super( params );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link MEManager} out of its XML representation.
	 * 
	 * @param stringBuff
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MEManager} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see InhomogeneousDGTrainSM#InhomogeneousDGTrainSM(StringBuffer)
	 */
	public MEManager( StringBuffer stringBuff ) throws NonParsableException
	{
		super( stringBuff );
	}

	public MEManager clone() throws CloneNotSupportedException
	{
		MEManager clone = (MEManager) super.clone();
		if( factors != null )
		{
			clone.factors = new MEM[factors.length];
			for( int i = 0; i < factors.length; i++ )
			{
				clone.factors[i] = factors[i].clone();
			}
		}
		return clone;
	}
	
	public DataSet emitDataSet( int n, int... lengths ) throws NotTrainedException, Exception
	{
		if( !(lengths == null || lengths.length == 0) )
		{
			throw new Exception( "This is an inhomogeneous model. Please check parameter lengths." );
		}
		return DiscreteInhomogenousDataSetEmitter.emitDataSet( this, n );
	}

	public double getLogPriorTerm() throws Exception
	{
		if( trained ){
			double ess = getESS();
			if( ess == 0 ){
				return 0;
			} else {
				if( (Decomposition) params.getParameterAt( 2 ).getValue() != Decomposition.DECOMPOSE_LESS_CONNECTED ){
					double res = 0;
					for( int i = 0; i < factors.length; i++ ){
						res += factors[i].getLogPriorPart( ess );
					}
					return res;
				} else {
					// TODO theory?
					throw new OperationNotSupportedException( "The prior can not be computed if the MEM is decomposed uding DECOMPOSE_LESS_CONNECTED." );
				}			
			}
		}
		else
		{
			throw new NotTrainedException();
		}
	}
	
	public double getLogProbFor( Sequence sequence, int startpos, int endpos ) throws NotTrainedException, Exception
	{
		check( sequence, startpos, endpos );
		double res = factors[0].getLogScoreFor( sequence, startpos );
		for( int i = 1; i < factors.length; i++ )
		{
			res += factors[i].getLogScoreFor( sequence, startpos );
		}
		return res;
	}

	public NumericalResultSet getNumericalCharacteristics()
	{
		/*
		if( factors != null )
		{
			int[] sum = { 0, 0 };
			int counter1 = 0, counter2;
			for( ; counter1 < factors.length; counter1++ )
			{
				sum[0] += factors[counter1].constraints.length;
				for( counter2 = 0; counter2 < factors[counter1].constraints.length; counter2++ )
				{
					sum[1] += factors[counter1].constraints[counter2].getNumberOfSpecificConstraints();
				}
			}
			return new NumericalResultSet( new NumericalResult[]{
					new NumericalResult( "number of constraints", "the number of used cliques", sum[0] ),
					new NumericalResult( "number of specific constraints", "the number of used parameters", sum[1] ) } );
		}
		else
		{
			throw new NotTrainedException();
		}
		*/
		return null;
	}

	public String getStructure() throws NotTrainedException
	{
		if( factors == null )
		{
			throw new NotTrainedException();
		}
		int counter1 = 0, counter2;
		String erg = "";
		while( counter1 < factors.length )
		{
			erg += "factor " + counter1 + ": " + factors[counter1].toString() + "\n";
			for( counter2 = 0; counter2 < factors[counter1].constraints.length; counter2++ )
			{
				erg += factors[counter1].constraints[counter2] + "\n";
			}
			counter1++;
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.DiscreteGraphicalTrainSM#toString(java.text.NumberFormat)
	 */
	@Override
	public String toString( NumberFormat nf )
	{
		String erg = super.toString(nf);
		if( factors != null && factors.length > 1 )
		{
			erg += "\n\n";
			for( int i = 0; i < factors.length; i++ )
			{
				erg += "MEM " + i + ": " + factors[i].toString(nf) + "\n";
			}
		}
		return erg;
	}

	/**
	 * This method returns an array of independent maximum entropy models parsed from the given constraints.
	 * 
	 * @param constraints the constraints to build the maximum entropy model
	 * @param reduce a switch whether redundant constraint should be removed
	 * @param decomposition a switch how to decompose the complete model if possible 
	 * 
	 * @return  an array of independent maximum entropy models
	 */
	protected MEM[] getFactors( String constraints, boolean reduce, Decomposition decomposition )
	{
		return getFactors( ConstraintManager.extract( length, constraints ), reduce, decomposition );
	}
	
	/**
	 * This method returns an array of independent maximum entropy models parsed from the given constraints.
	 * 
	 * @param list a list of positions arrays that build the constraints
	 * @param reduce a switch whether redundant constraint should be removed
	 * @param decomposition a switch how to decompose the complete model if possible 
	 * 
	 * @return  an array of independent maximum entropy models
	 */
	protected MEM[] getFactors( ArrayList<int[]> list, boolean reduce, Decomposition decomposition )
	{
		if( reduce )
		{
			ConstraintManager.reduce( list );
		}
		return ConstraintManager.disconnect( list, alphabetLength, decomposition );
	}

	@Override
	protected StringBuffer getFurtherModelInfos()
	{
		StringBuffer b = new StringBuffer( 10000 );
		XMLParser.appendObjectWithTags( b, sostream.doesNothing(), "outStream_does_nothing" );
		if( trained )
		{
			XMLParser.appendObjectWithTags( b, factors, "mem" );
		}
		return b;
	}

	/**
	 * This method trains the internal {@link MEM} array,
	 * i.e., it optimizes the parameters of the underlying {@link MEMConstraint}s.
	 *  
	 * @param data the data
	 * @param weights the weights for the data, can be <code>null</code>
	 * @throws Exception if some error occurs in the training process 
	 */
	protected void trainFactors( DataSet data, double[] weights ) throws Exception
	{
		int i = 0, offset = 0, l;
		// generate MEMConstraints-array
		while( i < factors.length )
		{
			offset += factors[i++].constraints.length;
		}
		MEMConstraint[] c = new MEMConstraint[offset];
		offset = 0;
		for( i = 0; i < factors.length; i++ )
		{
			l = factors[i].constraints.length;
			System.arraycopy( factors[i].constraints, 0, c, offset, l );
			offset += l;
		}
		// count
		ConstraintManager.countInhomogeneous( alphabets, length, data, weights, true, c );
		// estimate freq
		ConstraintManager.computeFreqs( getESS(), c );
		
		// optimize
		SequenceIterator s = new SequenceIterator( length );
		byte algo = (Byte) params.getParameterAt( 4 ).getValue();
		SmallDifferenceOfFunctionEvaluationsCondition eps = new SmallDifferenceOfFunctionEvaluationsCondition( (Double) params.getParameterAt( 5 ).getValue() );
		for( i = 0; i < factors.length; i++ )
		{
			factors[i].train( s, algo, eps, sostream );
		}
		System.gc();
	}

	@Override
	protected void setFurtherModelInfos( StringBuffer xml ) throws NonParsableException
	{
		if( XMLParser.extractObjectForTags( xml, "outStream_does_nothing", boolean.class ) )
		{
			setOutputStream( null );
		}
		else
		{
			setOutputStream( DEFAULT_STREAM );
		}
		if( trained )
		{
			this.factors = XMLParser.extractObjectForTags( xml, "mem", MEM[].class );
		}
	}
}
