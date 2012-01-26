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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.elements;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;


/**
 * Basic transition element without random initialization of parameters.
 * 
 * @author Michael Seifert
 */
public class BasicPluginTransitionElement extends BasicTransitionElement
{
	/**
	 * Probabilities assigned to the initial transition element.
	 */
	private double[] probs; 
	
	/**
	 * Constructor creating a new instance with given context, descendant states, and hyper parameters.
	 *  
	 * @param context previously visited state indices; last entry corresponds to the last state visited
	 * 
	 * @param states all possible states that can be reached from the context; if <code>null</code> then no transition allowed
	 * 
	 * @param probs initial probabilities 
	 * 
	 * @param hyperParameters hyper-parameters for transitions; if <code>null</code> then no prior is used
	 *
	 * @see #BasicPluginTransitionElement(int[], int[], double[], double[], double[])
	 */
	public BasicPluginTransitionElement( int[] context, int[] states, double[] probs, double[] hyperParameters )
	{
		this( context, states, probs, hyperParameters, null );
	}
	
	/**
	 * Constructor creating a new instance with given context, descendant states, and hyper parameters.
	 *  
	 * @param context previously visited state indices; last entry corresponds to the last state visited
	 * 
	 * @param states all possible states that can be reached from the context; if <code>null</code> then no transition allowed
	 * 
	 * @param probs initial probabilities 
	 * 
	 * @param hyperParameters hyper-parameters for transitions; if <code>null</code> then no prior is used
	 *
	 * @param weight the weight for plotting the edge in Graphviz, enables to modify the edge length, larger weights imply shorter edges (default: 1)
	 */
	public BasicPluginTransitionElement( int[] context, int[] states, double[] probs, double[] hyperParameters, double[] weight )
	{
		super( context, states, hyperParameters, weight );
		this.probs = probs.clone();
		this.initializeRandomly();		
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link BasicPluginTransitionElement} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link BasicPluginTransitionElement} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 *             
	 */
	public BasicPluginTransitionElement( StringBuffer xml ) throws NonParsableException
	{
		super( xml );
	}
	
	/*
	 * Non-random initialization by setting the probabilities to <code>probs</code>.
	 * 
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.hmm.transitions.BasicHigherOrderTransition.BasicTransitionElement#initializeRandomly()
	 */
	public void initializeRandomly()
	{		
		for( int i = 0; i < probs.length; i++ )
		{
			this.parameters[ i ] = Math.log( this.probs[ i ] );
		}
		this.logNorm = 0;
	}
	
	protected void appendFurtherInformation( StringBuffer xml )
	{
		XMLParser.appendObjectWithTags( xml, this.probs, "InitialProbabilities" );
	}

	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException
	{
		this.probs = (double[]) XMLParser.extractObjectForTags( xml, "InitialProbabilities" );
	}
	
	protected String getXMLTag()
	{
		return XML_TAG;
	}
	
	private static final String XML_TAG = "BasicPluginTransitionElement";

}
