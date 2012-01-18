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
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.trainable;

import java.io.IOException;
import java.util.Random;

import de.jstacs.NotTrainedException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.ArbitrarySequence;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResultSet;

/**
 * This class represents a uniform model. Sometimes it's also called uninformed model. It can be used if nothing is
 * known about a statistical process.
 * 
 * @author Jens Keilwagen
 */
public class UniformTrainSM extends AbstractTrainSM
{
	private static final long serialVersionUID = 1L;

	private static final String XML_TAG = "DiscreteUniformModel";

	private double p;
	
	/**
	 * Creates a new {@link UniformTrainSM} using a given {@link AlphabetContainer}.
	 * 
	 * @param alphabet the alphabets used in the model
	 */
	public UniformTrainSM( AlphabetContainer alphabet )
	{
		super( alphabet, alphabet.getPossibleLength() );
		init();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link UniformTrainSM} out of a {@link StringBuffer}.
	 * 
	 * @param stringBuff  the {@link StringBuffer} to be parsed
	 * 
	 * @throws NonParsableException  if the {@link StringBuffer} is not parsable
	 */
	public UniformTrainSM( StringBuffer stringBuff ) throws NonParsableException
	{
		super( stringBuff );
		init();
	}
	
	private void init() {
		if( length > 0 )
		{
			p = 0;
			for( int i = 0; i < length; i++ )
			{
				p -= Math.log( alphabets.getAlphabetLengthAt( i ) );
			}
		}
		else
		{
			p = -Math.log( alphabets.getAlphabetLengthAt( 0 ) );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#clone()
	 */
	@Override
	public UniformTrainSM clone() throws CloneNotSupportedException
	{
		return (UniformTrainSM) super.clone();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getProbFor(de.jstacs.data.Sequence, int, int)
	 */
	public double getLogProbFor( Sequence sequence, int startpos, int endpos ) throws IllegalArgumentException,
			WrongAlphabetException
	{
		if( !alphabets.checkConsistency( sequence.getAlphabetContainer() ) )
		{
			throw new WrongAlphabetException( "This sequence and model doesnot match." );
		}
		else if( 0 > startpos || endpos > sequence.getLength() )
		{
			throw new IllegalArgumentException( "Attention: 0 <= startpos, endpos < sequence.length()" );
		}
		else if( endpos - startpos + 1 < 0 )
		{
			throw new IllegalArgumentException(
					"The sequence does not have a correct length. The length has to be non-negative." );
		}
		else
		{
			if( length > 0 )
			{
				if( endpos - startpos + 1 == length )
				{
					return p;
				}
				else
				{
					throw new IllegalArgumentException( "The sequence does not have a correct length (" + length + ")." );
				}
			}
			else
			{
				return p * ( endpos - startpos + 1 );
			}
		}
	}
	
	/**
	 * Returns <code>true</code> if the model is trained, <code>false</code> otherwise.
	 * 
	 * @return <code>true</code> if the model is trained, <code>false</code> otherwise
	 */
	public boolean isInitialized()
	{
		return true;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#fromXML(java.lang.StringBuffer)
	 */
	@Override
	public void fromXML( StringBuffer representation ) throws NonParsableException
	{
		alphabets = new AlphabetContainer( XMLParser.extractForTag( representation, XML_TAG ) );
		length = alphabets.getPossibleLength();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML()
	{
		StringBuffer xml = alphabets.toXML();
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}

	/**
	 * Returns the String &quot;&quot;.
	 */
	@Override
	public String toString()
	{
		return "";
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#train(de.jstacs.data.DataSet, double[])
	 */
	@Deprecated
	public void train( DataSet data, double[] weights ) throws IOException
	{
		// nothing
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#emitDataSet(int, int[])
	 */
	@Override
	public DataSet emitDataSet( int n, int... lengths ) throws Exception
	{
		Sequence[] seq;
		if( length == 0 )
		{
			if( lengths.length == 1 )
			{
				seq = getSequences( n, lengths[0] );
			}
			else
			{
				seq = new Sequence[n];
				for( int i = 0; i < n; i++ )
				{
					seq[i] = getSequences( 1, lengths[i] )[0];
				}
			}
		}
		else
		{
			if( !(lengths == null || lengths.length == 0) )
			{
				throw new Exception( "This is an inhomogeneous model. Please check parameter lengths." );
			}
			double[] content = new double[length];
			seq = new Sequence[n];
			for( int j, i = 0; i < n; i++ )
			{
				for( j = 0; j < length; j++ )
				{
					content[j] = alphabets.getMin( j ) + r.nextDouble() * alphabets.getAlphabetLengthAt( j );
					if( alphabets.isDiscreteAt( j ) )
					{
						content[j] = (int) content[j];
					}
				}
				seq[i] = new ArbitrarySequence( alphabets, content );
			}
		}
		return new DataSet( "sampled from " + getInstanceName(), seq );
	}
	
	private static final Random r = new Random();

	private Sequence[] getSequences( int n, int length ) throws WrongAlphabetException,
			WrongSequenceTypeException
	{
		Sequence[] seqs = new Sequence[n];
		if( alphabets.isDiscrete() )
		{
			int[] seq = new int[length];
			int i, j = 0, l = (int) alphabets.getAlphabetLengthAt( 0 );
			while( j < n )
			{
				for( i = 0; i < length; i++ )
				{
					seq[i] = r.nextInt( l );
				}
				seqs[j++] = new IntSequence( alphabets, seq );
			}
		}
		else
		{
			double[] seq = new double[length];
			double m = alphabets.getMin( 0 ), l = alphabets.getAlphabetLengthAt( 0 );
			int i, j = 0;
			while( j < n )
			{
				for( i = 0; i < length; i++ )
				{
					seq[i] = m + r.nextDouble() * l;
				}
				seqs[j++] = new ArbitrarySequence( alphabets, seq );
			}
		}
		return seqs;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorTerm() throws Exception
	{
		return 0;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#getMaximalMarkovOrder()
	 */
	@Override
	public byte getMaximalMarkovOrder() throws UnsupportedOperationException
	{
		return (byte) 0;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getNumericalCharacteristics()
	 */
	public NumericalResultSet getNumericalCharacteristics() throws Exception
	{
		return null;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getInstanceName()
	 */
	public String getInstanceName()
	{
		return "uniform";
	}
}
