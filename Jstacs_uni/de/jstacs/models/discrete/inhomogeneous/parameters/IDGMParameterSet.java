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

package de.jstacs.models.discrete.inhomogeneous.parameters;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.models.discrete.DGMParameterSet;
import de.jstacs.models.discrete.inhomogeneous.InhomogeneousDGM;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.ModelType;

/**
 * This is the abstract container of parameters that is a root container for all
 * inhomogeneous discrete graphical model parameter containers.
 * 
 * @author Jens Keilwagen
 */
public abstract class IDGMParameterSet extends DGMParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link IDGMParameterSet} out of its XML representation.
	 * 
	 * @param s
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link IDGMParameterSet} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see DGMParameterSet#DGMParameterSet(StringBuffer)
	 */
	protected IDGMParameterSet( StringBuffer s ) throws NonParsableException {
		super( s );
	}

	/**
	 * This constructor creates an empty {@link IDGMParameterSet} instance from
	 * the class that can be instantiated using this {@link IDGMParameterSet}.
	 * 
	 * @param instanceClass
	 *            the instance class
	 * 
	 * @see DGMParameterSet#DGMParameterSet(Class, boolean, boolean)
	 */
	protected IDGMParameterSet( Class<? extends InhomogeneousDGM> instanceClass ) {
		super( instanceClass, false, false );
	}

	/**
	 * This constructor creates an {@link IDGMParameterSet} instance for the
	 * specified class. It sets the {@link AlphabetContainer}, the length, the
	 * ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize) and the model
	 * description.
	 * 
	 * @param instanceClass
	 *            the instance class
	 * @param alphabet
	 *            the {@link AlphabetContainer} for the model
	 * @param length
	 *            the length of the model
	 * @param ess
	 *            the ess of the model
	 * @param description
	 *            the model description
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected IDGMParameterSet( Class<? extends InhomogeneousDGM> instanceClass, AlphabetContainer alphabet, int length, double ess,
								String description ) throws Exception {
		super( instanceClass, alphabet, length, ess, description );
	}

	/**
	 * This method returns a short textual representation of the model instance.
	 * 
	 * @param model
	 *            the type of the model
	 * @param order
	 *            the order of the model
	 * @param method
	 *            the learning method
	 * @param ess
	 *            the used ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * 
	 * @return a short textual representation of the model instance
	 * 
	 * @see ModelType
	 * @see LearningType
	 */
	public static String getModelInstanceName( ModelType model, byte order, LearningType method, double ess ) {
		return model.name() + "(" + order + ") " + ( method == LearningType.ML_OR_MAP ? ( ess == 0 ? "ML" : "MAP" ) : method.name() );
	}
}
