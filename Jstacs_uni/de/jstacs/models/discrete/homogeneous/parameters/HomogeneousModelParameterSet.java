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

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.AlphabetContainerParameterSet;
import de.jstacs.models.discrete.DGMParameterSet;
import de.jstacs.models.discrete.homogeneous.HomogeneousModel;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * This class implements a container for all parameters of any homogeneous
 * model.
 * 
 * @author Jens Keilwagen
 */
public abstract class HomogeneousModelParameterSet extends DGMParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link HomogeneousModelParameterSet} out of its XML
	 * representation.
	 * 
	 * @param s
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link HomogeneousModelParameterSet} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see DGMParameterSet#DGMParameterSet(StringBuffer)
	 */
	protected HomogeneousModelParameterSet( StringBuffer s ) throws NonParsableException {
		super( s );
	}

	/**
	 * This is the constructor that creates an empty
	 * {@link HomogeneousModelParameterSet} from the class that can be
	 * instantiated using this {@link HomogeneousModelParameterSet}.
	 * 
	 * @param instanceClass
	 *            the (sub-)class
	 * @see de.jstacs.models.discrete.homogeneous.HomogeneousMM
	 * @see HomogeneousModel
	 * @see DGMParameterSet#DGMParameterSet(Class, boolean, boolean)
	 */
	protected HomogeneousModelParameterSet( Class<? extends HomogeneousModel> instanceClass ) {
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
	 * Creates a new {@link HomogeneousModelParameterSet} with
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
	 * @see DGMParameterSet#DGMParameterSet(Class, AlphabetContainer, double,
	 *      String)
	 */
	protected HomogeneousModelParameterSet( Class<? extends HomogeneousModel> instanceClass, AlphabetContainer alphabet, double ess,
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
	 * @see de.jstacs.models.discrete.DGMParameterSet#hasDefaultOrIsSet()
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
