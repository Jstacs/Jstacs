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

package de.jstacs.scoringFunctions;

import de.jstacs.NonParsableException;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.Mutable;
import de.jstacs.scoringFunctions.mix.AbstractMixtureScoringFunction;
import de.jstacs.scoringFunctions.mix.StrandScoringFunction;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class makes an unnormalized {@link NormalizableScoringFunction} to a normalized {@link NormalizableScoringFunction}.
 * However, the class allows to use only {@link NormalizableScoringFunction} that do not implement {@link VariableLengthScoringFunction}.
 * This class should be used only in cases when it is not possible to avoid its usage.
 * 
 * @author Jens Keilwagen
 */
public final class NormalizedScoringFunction extends AbstractNormalizableScoringFunction implements Mutable
{
	/**
	 * This method returns a normalized version of a NormalizableScoringFunction. Either it returns a clone of
	 * <code>nsf</code> or an instance of NormalizedScoringFunction using <code>nsf</code> and <code>starts</code>.
	 * 
	 * @param nsf
	 *            the NormalizableScoringFunction to be normalized
	 * @param starts
	 *            the number of recommended starts for a NormalizedScoringFunction
	 * 
	 * @return a normalized scoring function
	 * @throws Exception if <code>nsf</code> could not be cloned

	 */
	public static final NormalizableScoringFunction getNormalizedVersion( NormalizableScoringFunction nsf, int starts ) throws Exception
	{
		if( nsf.isNormalized() )
		{
			return (NormalizableScoringFunction) nsf.clone();
		}
		else
		{
			return new NormalizedScoringFunction( nsf, starts );
		}
	}

	private NormalizableScoringFunction nsf;

	private int starts;

	private double logNorm;

	private double[] proportion;

	/**
	 * Creates a new instance using a given NormalizableScoringFunction.
	 * 
	 * @param nsf
	 *            the function to be used internal
	 * @param starts
	 *            the number of recommended starts ({@link ScoringFunction#getNumberOfRecommendedStarts()})
	 * @throws Exception is <code>nsf</code> could not be cloned or some error occurred during computation of some values
	 * 
	 */
	public NormalizedScoringFunction( NormalizableScoringFunction nsf, int starts ) throws Exception
	{
		super( nsf.getAlphabetContainer(), nsf.getLength() );
		if( nsf instanceof VariableLengthScoringFunction ) {
			throw new IllegalArgumentException();
		}
		if( starts <= 0 )
		{
			throw new IllegalArgumentException( "The number of starts has to be positive." );
		}
		this.starts = Math.max( starts, nsf.getNumberOfRecommendedStarts() );
		this.nsf = (NormalizableScoringFunction) nsf.clone();
		precomputeIfPossible();
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}.
	 * 
	 * @param xml the xml representation
	 * 
	 * @throws NonParsableException if the representation could not be parsed.
	 */
	public NormalizedScoringFunction( StringBuffer xml ) throws NonParsableException
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

	public NormalizedScoringFunction clone() throws CloneNotSupportedException
	{
		NormalizedScoringFunction clone = (NormalizedScoringFunction) super.clone();
		clone.nsf = (NormalizableScoringFunction) nsf.clone();
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

	public void initializeFunction( int index, boolean freeParams, Sample[] data, double[][] weights ) throws Exception
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
		nsf = XMLParser.extractObjectForTags( b, "function", NormalizableScoringFunction.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
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
	 * @throws Exception if the partial normalization constants can not be computed {@link NormalizableScoringFunction#getLogPartialNormalizationConstant(int)}
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
	public NormalizableScoringFunction getFunction() throws CloneNotSupportedException
	{
		return (NormalizableScoringFunction) nsf.clone();
	}

	/*
	public int[] determineNotSignificantPositions( double samples, double[] weightsLeft, double[] weightsRight, double[][][][] contrastLeft, double[][][][] contrastRight, double sign )
	{
		if( nsf instanceof Mutable )
		{
			return ((Mutable)nsf).determineNotSignificantPositions( samples, weightsLeft, weightsRight, contrastLeft, contrastRight, sign );
		}
		else
		{
			return new int[]{0,0};
		}
	}

	public boolean modify( double[] weightsLeft, double[] weightsRight, double[][][][] replacementLeft, double[][][][] replacementRight, int offsetLeft, int offsetRight )
	{
		if( nsf instanceof Mutable )
		{
			boolean modified = ((Mutable)nsf).modify( weightsLeft, weightsRight, replacementLeft, replacementRight, offsetLeft, offsetRight );
			if( modified )
			{
				if( proportion.length != nsf.getNumberOfParameters() )
				{
					proportion = new double[nsf.getNumberOfParameters()];
				}
				try
				{
					precompute();
				}
				catch( Exception e )
				{
					RuntimeException r = new RuntimeException( e.getMessage() );
					r.setStackTrace(e.getStackTrace());
					throw r;
				}
				length = nsf.getLength();
			}
			return modified;
		}
		else
		{
			return false;
		}
	}
	*/
	
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
	 * This method returns <code>true</code> if the internal {@link NormalizableScoringFunction} is a {@link StrandScoringFunction} otherwise <code>false</code>.
	 * 
	 * @return <code>true</code> if the internal {@link NormalizableScoringFunction} is a {@link StrandScoringFunction} otherwise <code>false</code>
	 */
	public boolean isStrandScoringFunction()
	{
		if( nsf instanceof NormalizedScoringFunction )
		{
			return ((NormalizedScoringFunction)nsf).isStrandScoringFunction();
		}
		else
		{
			return nsf instanceof StrandScoringFunction;
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
		if( nsf instanceof NormalizedScoringFunction )
		{
			return ((NormalizedScoringFunction)nsf).getStrand( seq, startPos );
		}
		else
		{
			if( nsf instanceof StrandScoringFunction )
			{
				return ((StrandScoringFunction)nsf).getStrand( seq, startPos );
			}
			else
			{
				return Strand.FORWARD;
			}
		}
	}
	
	/**
	 * This method initializes the hidden parameters of the internal {@link NormalizableScoringFunction} uniformly if it is a {@link AbstractMixtureScoringFunction}.
	 */
	public void initializeHiddenUniformly()
	{
		if( nsf instanceof NormalizedScoringFunction )
		{
			((NormalizedScoringFunction)nsf).initializeHiddenUniformly();
		}
		else if( nsf instanceof AbstractMixtureScoringFunction )
		{
			((AbstractMixtureScoringFunction)nsf).initializeHiddenUniformly();
		}
	}
}
