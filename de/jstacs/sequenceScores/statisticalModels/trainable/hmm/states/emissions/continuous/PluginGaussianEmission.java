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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.continuous;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * Basic Gaussian emission distribution without random initialization of parameters.
 * 
 * @author Michael Seifert
 *
 */
public class PluginGaussianEmission extends GaussianEmission
{
	/**
	 * Initial mean value.
	 */
	protected double meanValue;
	
	/**
	 * Initial standard deviation.
	 */
	protected double standardDeviation;
	
	private static final String TAG = "GaussianEmission_NoRandomInit";
	
	/**
	 * Creates a Gaussian emission density with mean <code>mean</code> and standard deviation <code>sd</code>.
	 * 
	 * @param mean mean value
	 * 
	 * @param sd standard deviation (initially handled as precision)
	 * 
	 * @param ess scale factor of a priori mean (epsilon)
	 * 
	 * @param priorMu a priori mean (nhi)
	 * 
	 * @param priorAlpha shape parameter (r)
	 * 
	 * @param priorBeta scale parameter (alpha)
	 */
	public PluginGaussianEmission( double mean, double sd, double ess, double priorMu, double priorAlpha, double priorBeta )
	{
		super( new AlphabetContainer( new ContinuousAlphabet() ), ess, priorMu, priorAlpha, priorBeta, false );
		
		this.meanValue = mean;
		this.standardDeviation = sd;
				
		//Set mean and log-precision of the Gaussian emission density
		this.setParameter( new double[]{ this.meanValue, Math.log( 1 / Math.pow( this.standardDeviation, 2 ) ) }, 0 );
				
		//Pre-computations
		this.precompute();
	}
	
	/**
	 * Creates a {@link PluginGaussianEmission} from its XML representation.
	 * @param xml the XML representation.
	 * @throws NonParsableException if the XML representation could not be parsed
	 */
	public PluginGaussianEmission( StringBuffer xml ) throws NonParsableException
	{
		super( xml );
	}
	
	/*
	 * Non-random initialization with pre-defined mean and pre-defined precision.
	 * 
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.GaussianEmission#initializeFunctionRandomly()
	 */
	public void initializeFunctionRandomly()
	{
		this.setParameter( new double[]{ this.meanValue, Math.log( 1 / Math.pow( this.standardDeviation, 2 ) ) }, 0 );
		precompute();
	}
	
	public String toString()
	{
		double[] parameters = new double[ 2 ];
		this.fillCurrentParameter( parameters ); //get mean and log-precision
		return( "   - Mean = " + parameters[ 0 ] + "\t Sd = " + Math.pow( 1 / Math.exp( parameters[ 1 ] ), 0.5 ) );
	}
	
	
	public StringBuffer toXML()
	{
		StringBuffer buf = super.toXML(); 
		
		XMLParser.appendObjectWithTags( buf, this.meanValue, "MeanValue" );
		XMLParser.appendObjectWithTags( buf, this.standardDeviation, "StandardDeviation" );
		
		//XMLParser.addTags( buf, getClass().getSimpleName() );
		XMLParser.addTags( buf, TAG );		
		
		return buf;
	}
	
	protected void fromXML( StringBuffer xml ) throws NonParsableException
	{		
		//xml = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		xml = XMLParser.extractForTag( xml, TAG );
		
		this.meanValue         = XMLParser.extractObjectForTags( xml, "MeanValue", double.class );
		this.standardDeviation = XMLParser.extractObjectForTags( xml, "StandardDeviation", double.class );
		
		super.fromXML( xml );		
	}
}
