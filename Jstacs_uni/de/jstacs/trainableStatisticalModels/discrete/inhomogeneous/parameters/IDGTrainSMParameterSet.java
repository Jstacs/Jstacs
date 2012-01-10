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

package de.jstacs.trainableStatisticalModels.discrete.inhomogeneous.parameters;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.trainableStatisticalModels.discrete.DGTrainSMParameterSet;
import de.jstacs.trainableStatisticalModels.discrete.inhomogeneous.InhomogeneousDGTrainSM;
import de.jstacs.trainableStatisticalModels.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.trainableStatisticalModels.discrete.inhomogeneous.StructureLearner.ModelType;

/**
 * This is the abstract container of parameters that is a root container for all
 * inhomogeneous discrete graphical model parameter containers.
 * 
 * @author Jens Keilwagen
 */
public abstract class IDGTrainSMParameterSet extends DGTrainSMParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link IDGTrainSMParameterSet} out of its XML representation.
	 * 
	 * @param s
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link IDGTrainSMParameterSet} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see DGTrainSMParameterSet#DGMParameterSet(StringBuffer)
	 */
	protected IDGTrainSMParameterSet( StringBuffer s ) throws NonParsableException {
		super( s );
	}

	/**
	 * This constructor creates an empty {@link IDGTrainSMParameterSet} instance from
	 * the class that can be instantiated using this {@link IDGTrainSMParameterSet}.
	 * 
	 * @param instanceClass
	 *            the instance class
	 * 
	 * @see DGTrainSMParameterSet#DGMParameterSet(Class, boolean, boolean)
	 */
	protected IDGTrainSMParameterSet( Class<? extends InhomogeneousDGTrainSM> instanceClass ) {
		super( instanceClass, false, false );
	}

	/**
	 * This constructor creates an {@link IDGTrainSMParameterSet} instance for the
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
	protected IDGTrainSMParameterSet( Class<? extends InhomogeneousDGTrainSM> instanceClass, AlphabetContainer alphabet, int length, double ess,
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
