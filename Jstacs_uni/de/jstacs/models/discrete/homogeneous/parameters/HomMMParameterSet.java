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

package de.jstacs.models.discrete.homogeneous.parameters;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.models.discrete.homogeneous.HomogeneousMM;

/**
 * This class implements a container for all parameters of a homogeneous Markov
 * model.
 * 
 * @author Jens Keilwagen
 */
public class HomMMParameterSet extends HomogeneousModelParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link HomMMParameterSet} out of its XML representation.
	 * 
	 * @param s
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link HomMMParameterSet} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see HomogeneousModelParameterSet#HomogeneousModelParameterSet(StringBuffer)
	 */
	public HomMMParameterSet( StringBuffer s ) throws NonParsableException {
		super( s );
	}

	/**
	 * An empty constructor. Creates a new {@link HomMMParameterSet}.
	 */
	public HomMMParameterSet() {
		super( HomogeneousMM.class );
	}

	/**
	 * Creates a new {@link HomMMParameterSet} with {@link AlphabetContainer},
	 * ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize), description and order
	 * of the homogeneous Markov model.
	 * 
	 * @param alphabet
	 *            the {@link AlphabetContainer}
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * @param description
	 *            the description
	 * @param order
	 *            the order of the Markov model
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see HomogeneousModelParameterSet#HomogeneousModelParameterSet(Class,
	 *      AlphabetContainer, double, String, byte)
	 */
	public HomMMParameterSet( AlphabetContainer alphabet, double ess, String description, byte order ) throws Exception {
		super( HomogeneousMM.class, alphabet, ess, description, order );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
	 */
	@Override
	public String getInstanceComment() {
		return "a homogeneous Markov model with user-specified order";
	}
}
