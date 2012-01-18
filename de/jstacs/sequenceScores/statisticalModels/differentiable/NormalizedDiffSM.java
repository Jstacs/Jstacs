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

package de.jstacs.sequenceScores.statisticalModels.differentiable;

import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.Mutable;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.AbstractMixtureDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.StrandDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class makes an unnormalized {@link DifferentiableStatisticalModel} to a normalized {@link DifferentiableStatisticalModel}.
 * However, the class allows to use only {@link DifferentiableStatisticalModel} that do not implement {@link VariableLengthDiffSM}.
 * This class should be used only in cases when it is not possible to avoid its usage.
 * 
 * @author Jens Keilwagen
 */
public final class NormalizedDiffSM extends AbstractDifferentiableStatisticalModel implements Mutable
{
	/**
	 * This method returns a normalized version of a DifferentiableStatisticalModel. Either it returns a clone of
	 * <code>nsf</code> or an instance of NormalizedDiffSM using <code>nsf</code> and <code>starts</code>.
	 * 
	 * @param nsf
	 *            the DifferentiableStatisticalModel to be normalized
	 * @param starts
	 *            the number of recommended starts for a NormalizedDiffSM
	 * 
	 * @return a normalized scoring function
	 * @throws Exception if <code>nsf</code> could not be cloned

	 */
	public static final DifferentiableStatisticalModel getNormalizedVersion( DifferentiableStatisticalModel nsf, int starts ) throws Exception
	{
		if( nsf.isNormalized() )
		{
			return (DifferentiableStatisticalModel) nsf.clone();
		}
		else
		{
			return new NormalizedDiffSM( nsf, starts );
		}
	}

	private DifferentiableStatisticalModel nsf;

	private int starts;

	private double logNorm;

	private double[] proportion;

	/**
	 * Creates a new instance using a given DifferentiableStatisticalModel.
	 * 
	 * @param nsf
	 *            the function to be used internal
	 * @param starts
	 *            the number of recommended starts ({@link de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore#getNumberOfRecommendedStarts()})
	 * @throws Exception is <code>nsf</code> could not be cloned or some error occurred during computation of some values
	 * 
	 */
	public NormalizedDiffSM( DifferentiableStatisticalModel nsf, int starts ) throws Exception
	{
		super( nsf.getAlphabetContainer(), nsf.getLength() );
		if( nsf instanceof VariableLengthDiffSM ) {
			throw new IllegalArgumentException();
		}
		if( starts <= 0 )
		{
			throw new IllegalArgumentException( "The number of starts has to be positive." );
		}
		this.starts = Math.max( starts, nsf.getNumberOfRecommendedStarts() );
		this.nsf = (DifferentiableStatisticalModel) nsf.clone();
		precomputeIfPossible();
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}.
	 * 
	 * @param xml the xml representation
	 * 
	 * @throws NonParsableException if the representation could not be parsed.
	 */
	public NormalizedDiffSM( StringBuffer xml ) throws NonParsableException
	{
		super( xml );
		try {
			precomputeIfPossible();
		} catch ( Exception e ) {
			NonParsableException n = new NonParsableException( e.getMessage() );
			n.setStackTrace( e.getStackTrace() );
			throw n;
		}
	}

	public NormalizedDiffSM clone() throws CloneNotSupportedException
	{
		NormalizedDiffSM clone = (NormalizedDiffSM) super.clone();
		clone.nsf = (DifferentiableStatisticalModel) nsf.clone();
		return clone;
	}

	public int getSizeOfEventSpaceForRandomVariablesOfParameter( int index )
	{
		return nsf.getSizeOfEventSpaceForRandomVariablesOfParameter( index );
	}

	public double getLogNormalizationConstant()
	{
		return 0;
	}

	public double getLogPartialNormalizationConstant( int parameterIndex ) throws Exception
	{
		return Double.NEGATIVE_INFINITY;
	}

	public double getESS()
	{
		return nsf.getESS();
	}

	public double getLogPriorTerm()
	{
		return nsf.getLogPriorTerm() - nsf.getESS() * logNorm;
	}

	public void addGradientOfLogPriorTerm( double[] grad, int start ) throws Exception
	{
		nsf.addGradientOfLogPriorTerm( grad, start );
		double e = nsf.getESS();
		for( int i = 0; i < proportion.length; i++ )
		{
			grad[start + i] -= e * proportion[i];
		}
	}

	public void initializeFunction( int index, boolean freeParams, DataSet[] data, double[][] weights ) throws Exception
	{
		nsf.initializeFunction( index, freeParams, data, weights );
		precompute();
	}

	public void initializeFunctionRandomly( boolean freeParams ) throws Exception
	{
		nsf.initializeFunctionRandomly( freeParams );
		precompute();
	}

	protected void fromXML( StringBuffer xml ) throws NonParsableException
	{
		StringBuffer b = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		nsf = XMLParser.extractObjectForTags( b, "function", DifferentiableStatisticalModel.class );
		alphabets = nsf.getAlphabetContainer();
		length = nsf.getLength();
		starts = XMLParser.extractObjectForTags( b, "starts", int.class );
	}

	public String getInstanceName()
	{
		return "normalized " + nsf.getInstanceName();
	}

	public double getLogScoreFor( Sequence seq, int start )
	{
		return nsf.getLogScoreFor( seq, start ) - logNorm;
	}

	public double getLogScoreAndPartialDerivation( Sequence seq, int start, IntList indices, DoubleList partialDer )
	{
		double score = nsf.getLogScoreAndPartialDerivation( seq, start, indices, partialDer ) - logNorm;
		for( int i = 0; i < proportion.length; i++ )
		{
			indices.add( i );
			partialDer.add( -proportion[i] );
		}
		return score;
	}

	public int getNumberOfParameters()
	{
		return nsf.getNumberOfParameters();
	}

	public double[] getCurrentParameterValues() throws Exception
	{
		return nsf.getCurrentParameterValues();
	}

	public void setParameters( double[] params, int start )
	{
		nsf.setParameters( params, start );
		try
		{
			precompute();
		}
		catch( Exception e )
		{
			RuntimeException r = new RuntimeException( e.getMessage() );
			r.setStackTrace( e.getStackTrace() );
			throw r;
		}
	}

	/**
	 * Trys to precompute the proportion of each parameter and the log of the normalization of nsf.
	 * 
	 * @throws Exception if the partial normalization constants can not be computed {@link DifferentiableStatisticalModel#getLogPartialNormalizationConstant(int)}
	 */
	private void precomputeIfPossible() throws Exception
	{
		if( nsf.isInitialized() )
		{
			precompute();
		}
		else
		{
			logNorm = Double.NEGATIVE_INFINITY;
			proportion = null;
		}
	}

	/**
	 * Precomputes the proportion of each parameter and the log of the normalization of nsf.
	 * 
	 * @throws Exception
	 */
	private void precompute() throws Exception
	{
		if( proportion == null )
		{
			proportion = new double[nsf.getNumberOfParameters()];
		}
		logNorm = nsf.getLogNormalizationConstant();
		for( int i = 0; i < proportion.length; i++ )
		{
			proportion[i] = Math.exp( nsf.getLogPartialNormalizationConstant( i ) - logNorm );
		}
	}

	public boolean isInitialized()
	{
		return nsf.isInitialized();
	}

	public StringBuffer toXML()
	{
		StringBuffer xml = new StringBuffer( 100000 );
		XMLParser.appendObjectWithTags( xml, nsf, "function" );
		XMLParser.appendObjectWithTags( xml, starts, "starts" );
		XMLParser.addTags( xml, getClass().getSimpleName() );
		return xml;
	}

	public int getNumberOfRecommendedStarts()
	{
		return starts;
	}

	public boolean isNormalized()
	{
		return true;
	}

	public String toString()
	{
		return "normalized variante of\n" + nsf.toString();
	}

	/**
	 * This method returns the internal function.
	 * 
	 * @return the internal function
	 * 
	 * @throws CloneNotSupportedException
	 *             if the internal function could not be cloned
	 */
	public DifferentiableStatisticalModel getFunction() throws CloneNotSupportedException
	{
		return (DifferentiableStatisticalModel) nsf.clone();
	}
	
	public boolean modify( int offsetLeft, int offsetRight ) {
		if( nsf instanceof Mutable ) {
			boolean res = ((Mutable) nsf).modify( offsetLeft, offsetRight );
			if( res ) {
				proportion = null;
				length = nsf.getLength();
				try
				{
					precompute();
				}
				catch( Exception e )
				{
					RuntimeException r = new RuntimeException( e.getMessage() );
					r.setStackTrace( e.getStackTrace() );
					throw r;
				}
			}
			return res;
		} else {
			return false;
		}
	}
	
	/**
	 * This method returns <code>true</code> if the internal {@link DifferentiableStatisticalModel} is a {@link StrandDiffSM} otherwise <code>false</code>.
	 * 
	 * @return <code>true</code> if the internal {@link DifferentiableStatisticalModel} is a {@link StrandDiffSM} otherwise <code>false</code>
	 */
	public boolean isStrandModel()
	{
		if( nsf instanceof NormalizedDiffSM )
		{
			return ((NormalizedDiffSM)nsf).isStrandModel();
		}
		else
		{
			return nsf instanceof StrandDiffSM;
		}
	}
	
	/**
	 * This method return the preferred {@link Strand} for a {@link Sequence} beginning at <code>startPos</code>.
	 * 
	 * @param seq the sequence
	 * @param startPos the start position
	 *  
	 * @return the preferred {@link Strand}
	 */
	public Strand getStrand( Sequence seq, int startPos )
	{
		if( nsf instanceof NormalizedDiffSM )
		{
			return ((NormalizedDiffSM)nsf).getStrand( seq, startPos );
		}
		else
		{
			if( nsf instanceof StrandDiffSM )
			{
				return ((StrandDiffSM)nsf).getStrand( seq, startPos );
			}
			else
			{
				return Strand.FORWARD;
			}
		}
	}
	
	/**
	 * This method initializes the hidden parameters of the internal {@link DifferentiableStatisticalModel} uniformly if it is a {@link AbstractMixtureDiffSM}.
	 */
	public void initializeHiddenUniformly()
	{
		if( nsf instanceof NormalizedDiffSM )
		{
			((NormalizedDiffSM)nsf).initializeHiddenUniformly();
		}
		else if( nsf instanceof AbstractMixtureDiffSM )
		{
			((AbstractMixtureDiffSM)nsf).initializeHiddenUniformly();
		}
	}
}
