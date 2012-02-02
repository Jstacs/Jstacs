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
package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.io.NonParsableException;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.FSDAGModelForGibbsSampling;

/**
 * The class for the parameters of a {@link FSDAGModelForGibbsSampling}.
 * 
 * @author Jens Keilwagen
 * 
 * @see FSDAGModelForGibbsSampling
 */
public class FSDAGModelForGibbsSamplingParameterSet extends FSDAGTrainSMParameterSet {

	/**
	 * The constructor for the {@link de.jstacs.Storable} interface. Creates a
	 * new {@link FSDAGModelForGibbsSamplingParameterSet} out of its XML
	 * representation.
	 * 
	 * @param s
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} is not parsable
	 */
	public FSDAGModelForGibbsSamplingParameterSet( StringBuffer s ) throws NonParsableException {
		super( s );
	}

	/**
	 * This is the constructor that creates an empty parameter set for a
	 * {@link FSDAGModelForGibbsSampling}.
	 */
	public FSDAGModelForGibbsSamplingParameterSet() {
		super( FSDAGModelForGibbsSampling.class );
	}

	/**
	 * This is the constructor that creates a filled parameter set.
	 * 
	 * @param alphabet
	 *            the alphabet container that is used in the model
	 * @param length
	 *            the length of the model (has to be positive)
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize, has to
	 *            be positive)
	 * @param description
	 *            a short description of the model (used for a better handling
	 *            of the object by the user)
	 * @param graph
	 *            the graph description string, encodes in XML-like manner the
	 *            parents of each node &quot;&lt;parents
	 *            node=i&gt;j,k,l&lt;/parents&gt;&quot;
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see FSDAGTrainSMParameterSet#encode(int[][])
	 */
	public FSDAGModelForGibbsSamplingParameterSet( AlphabetContainer alphabet, int length, double ess, String description, String graph )
																																			throws Exception {
		super( FSDAGModelForGibbsSampling.class, alphabet, length, ess, description, graph );
	}

}
