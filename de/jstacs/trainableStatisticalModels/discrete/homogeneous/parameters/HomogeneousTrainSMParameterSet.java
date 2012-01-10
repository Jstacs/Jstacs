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

package de.jstacs.trainableStatisticalModels.discrete.homogeneous.parameters;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.AlphabetContainerParameterSet;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.trainableStatisticalModels.discrete.DGTrainSMParameterSet;
import de.jstacs.trainableStatisticalModels.discrete.homogeneous.HomogeneousTrainSM;

/**
 * This class implements a container for all parameters of any homogeneous
 * model.
 * 
 * @author Jens Keilwagen
 */
public abstract class HomogeneousTrainSMParameterSet extends DGTrainSMParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link HomogeneousTrainSMParameterSet} out of its XML
	 * representation.
	 * 
	 * @param s
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link HomogeneousTrainSMParameterSet} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see DGTrainSMParameterSet#DGMParameterSet(StringBuffer)
	 */
	protected HomogeneousTrainSMParameterSet( StringBuffer s ) throws NonParsableException {
		super( s );
	}

	/**
	 * This is the constructor that creates an empty
	 * {@link HomogeneousTrainSMParameterSet} from the class that can be
	 * instantiated using this {@link HomogeneousTrainSMParameterSet}.
	 * 
	 * @param instanceClass
	 *            the (sub-)class
	 * @see de.jstacs.trainableStatisticalModels.discrete.homogeneous.HomogeneousMM
	 * @see HomogeneousTrainSM
	 * @see DGTrainSMParameterSet#DGMParameterSet(Class, boolean, boolean)
	 */
	protected HomogeneousTrainSMParameterSet( Class<? extends HomogeneousTrainSM> instanceClass ) {
		super( instanceClass, true, true );
		try {
			parameters.add( new SimpleParameter( DataType.BYTE,
					"order",
					"the order of the model specifies the number of used ancestors of a random variable that are used to determine its propability",
					true,
					new NumberValidator<Byte>( (byte)0, Byte.MAX_VALUE ) ) );
		} catch ( DatatypeNotValidException doesnothappen ) { }
	}

	/**
	 * Creates a new {@link HomogeneousTrainSMParameterSet} with
	 * {@link AlphabetContainer}, ess (<b>e</b>quivalent <b>s</b>ample
	 * <b>s</b>ize), description and order of the homogeneous Markov model. This
	 * constructor is for models that can handle variable lengths.
	 * 
	 * @param instanceClass
	 *            the (sub-)class
	 * @param alphabet
	 *            the {@link AlphabetContainer}
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * @param description
	 *            the description
	 * @param order
	 *            the order of the model
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see DGTrainSMParameterSet#DGMParameterSet(Class, AlphabetContainer, double,
	 *      String)
	 */
	protected HomogeneousTrainSMParameterSet( Class<? extends HomogeneousTrainSM> instanceClass, AlphabetContainer alphabet, double ess,
											String description, byte order ) throws Exception {
		super( instanceClass, alphabet, ess, description );
		parameters.add( new SimpleParameter( DataType.BYTE,
				"order",
				"the order of the model specifies the number of used ancestors of a random variable that are used to determine its propability",
				true,
				new NumberValidator<Byte>( (byte)0, Byte.MAX_VALUE ) ) );
		parameters.get( 2 ).setValue( new Byte( order ) );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.discrete.DGTrainSMParameterSet#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		if( super.hasDefaultOrIsSet() ) {
			//return getAlphabet().isSimple();
			return ( (AlphabetContainerParameterSet)alphabet.getValue() ).isSimple();
		} else {
			return false;
		}
	}
}
