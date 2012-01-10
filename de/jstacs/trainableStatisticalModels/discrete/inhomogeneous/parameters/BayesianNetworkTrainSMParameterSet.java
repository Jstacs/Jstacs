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

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.trainableStatisticalModels.discrete.inhomogeneous.BayesianNetworkTrainSM;
import de.jstacs.trainableStatisticalModels.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.trainableStatisticalModels.discrete.inhomogeneous.StructureLearner.ModelType;

/**
 * The {@link de.jstacs.parameters.ParameterSet} for the class
 * {@link BayesianNetworkTrainSM}.
 * 
 * @author Jens Keilwagen
 */
public class BayesianNetworkTrainSMParameterSet extends IDGTrainSMParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link BayesianNetworkTrainSMParameterSet} out of its XML
	 * representation.
	 * 
	 * @param s
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link BayesianNetworkTrainSMParameterSet} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see IDGTrainSMParameterSet#IDGTrainSMParameterSet(StringBuffer)
	 */
	public BayesianNetworkTrainSMParameterSet( StringBuffer s ) throws NonParsableException {
		super( s );
	}

	/**
	 * The simple constructor for an empty
	 * {@link BayesianNetworkTrainSMParameterSet} for a
	 * {@link BayesianNetworkTrainSM}.
	 * 
	 * @see IDGTrainSMParameterSet#IDGTrainSMParameterSet(Class)
	 */
	public BayesianNetworkTrainSMParameterSet() {
		super( BayesianNetworkTrainSM.class );
		addParameters();
	}

	/**
	 * This is the constructor of a filled
	 * {@link BayesianNetworkTrainSMParameterSet} for a
	 * {@link BayesianNetworkTrainSM}.
	 * 
	 * @param alphabet
	 *            the {@link AlphabetContainer} that is used in the model
	 * @param length
	 *            the length of the model (has to be positive)
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize) of the
	 *            model (has to be positive)
	 * @param description
	 *            a short description of the model (for a better handling of the
	 *            object by the user)
	 * @param model
	 *            the type of model: IMM, PMM or BN
	 * @param order
	 *            the order of the model
	 * @param method
	 *            the method how to learn the structure (only relevant for PMM,
	 *            BN): ML_OR_MAP or BMA
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see ModelType
	 * @see ModelType#IMM
	 * @see ModelType#PMM
	 * @see ModelType#BN
	 * @see LearningType
	 * @see LearningType#ML_OR_MAP
	 * @see LearningType#BMA
	 * @see IDGTrainSMParameterSet#IDGTrainSMParameterSet(Class, AlphabetContainer, int,
	 *      double, String)
	 */
	public BayesianNetworkTrainSMParameterSet( AlphabetContainer alphabet, int length, double ess, String description, ModelType model,
												byte order, LearningType method ) throws Exception {
		super( BayesianNetworkTrainSM.class, alphabet, length, ess, description );
		addParameters();
		parameters.get( 2 ).setValue( model );
		parameters.get( 3 ).setValue( order );
		parameters.get( 4 ).setValue( method );
	}

	private void addParameters() {
		try {
			parameters.add( new EnumParameter( ModelType.class, "the standard model that should be learned", true ) );

			parameters.add( new SimpleParameter( DataType.BYTE,
					"Markov order",
					"the used markov order is the number of used dependencies for (each) random variable",
					true,
					new NumberValidator<Byte>( (byte)0, Byte.MAX_VALUE ) ) );
			parameters.add( new EnumParameter( LearningType.class, "the learning method for the parameters of the model", true ) );
		} catch ( ParameterException doesnothappen ) { }
	}

	/**
	 * This method allows a simple change of the model type.
	 * 
	 * @param modelType
	 *            the type of the model, one of &quot;iMM&quot;, &quot;pMM&quot;
	 *            or &quot;BN&quot;
	 * 
	 * @throws IllegalValueException
	 *             if the <code>modelType</code> is illegal
	 */
	public void setModelType( String modelType ) throws IllegalValueException {
		parameters.get( 2 ).setValue( modelType );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
	 */
	@Override
	public String getInstanceComment() {
		return "a Bayesian Network model (i.e. inhomogeneous Markov model (iMM), permuted Markov model (pMM) or Bayesian network (BN)) with user-specified order";
	}

	/**
	 * This method returns a short description of the model.
	 * 
	 * @return a short description of the model
	 */
	public String getModelInstanceName() {
		return getModelInstanceName( (ModelType)parameters.get( 2 ).getValue(),
				(Byte)parameters.get( 3 ).getValue(),
				(LearningType)parameters.get( 4 ).getValue(),
				(Double)parameters.get( 0 ).getValue() );
	}
}
