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

package de.jstacs.models.discrete;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.AlphabetContainerParameterSet;
import de.jstacs.data.AlphabetContainer.AlphabetContainerType;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SequenceScoringParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * The super {@link de.jstacs.parameters.ParameterSet} for any parameter set of
 * a {@link DiscreteGraphicalModel}.
 * 
 * @author Jens Keilwagen
 */
public abstract class DGMParameterSet extends SequenceScoringParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link DiscreteGraphicalModel} out of its XML
	 * representation.
	 * 
	 * @param s
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DiscreteGraphicalModel} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see SequenceScoringParameterSet#SequenceScoringParameterSet(StringBuffer)
	 */
	protected DGMParameterSet( StringBuffer s ) throws NonParsableException {
		super( s );
	}

	/**
	 * An empty constructor. Creates a new {@link DGMParameterSet} from the
	 * class that can be instantiated using this {@link DGMParameterSet}.
	 * 
	 * @param instanceClass
	 *            the (sub-)class
	 * @param simple
	 *            indicates whether the alphabet should be simple or not
	 * @param variableLength
	 *            indicates whether the model can handle sequences of variable
	 *            length
	 * 
	 * @see SequenceScoringParameterSet#SequenceScoringParameterSet(Class,
	 *      AlphabetContainer.AlphabetContainerType, boolean, boolean)
	 */
	protected DGMParameterSet( Class<? extends DiscreteGraphicalModel> instanceClass, boolean simple, boolean variableLength ) {
		super( instanceClass, AlphabetContainerType.DISCRETE, simple, variableLength );
		addParameters();
	}

	/**
	 * The constructor for models that can handle variable lengths. Creates a
	 * new {@link DGMParameterSet} with a <code>description</code> from the
	 * class that can be instantiated using this {@link DGMParameterSet}, a
	 * given {@link AlphabetContainer} and a given <b>e</b>quivalent
	 * <b>s</b>ample <b>s</b>ize (ess).
	 * 
	 * @param instanceClass
	 *            the (sub-)class
	 * @param alphabet
	 *            the alphabet
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * @param description
	 *            the description
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected DGMParameterSet( Class<? extends DiscreteGraphicalModel> instanceClass, AlphabetContainer alphabet, double ess,
								String description ) throws Exception {
		this( instanceClass, alphabet, 0, true, ess, description );
	}

	/**
	 * The constructor for models that can handle only sequences of fixed length
	 * given by <code>length</code>. Creates a new {@link DGMParameterSet} with
	 * a <code>description</code> from the class that can be instantiated using
	 * this {@link DGMParameterSet}, a given {@link AlphabetContainer} and a
	 * given <b>e</b>quivalent <b>s</b>ample <b>s</b>ize (ess).
	 * 
	 * @param instanceClass
	 *            the (sub-)class
	 * @param alphabet
	 *            the alphabet
	 * @param length
	 *            the length of the modeled sequences
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * @param description
	 *            the description
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected DGMParameterSet( Class<? extends DiscreteGraphicalModel> instanceClass, AlphabetContainer alphabet, int length, double ess,
								String description ) throws Exception {
		this( instanceClass, alphabet, length, false, ess, description );
	}

	private DGMParameterSet( Class<? extends DiscreteGraphicalModel> instanceClass, AlphabetContainer alphabet, int length,
								boolean variableLength, double ess, String description ) throws Exception {
		super( instanceClass, alphabet, length, variableLength );
		addParameters();
		setEss( ess );
		if( description != null ) {
			parameters.get( 1 ).setValue( description );
		}
	}

	private void addParameters(){
		try{
			parameters.add( new SimpleParameter( DataType.DOUBLE,
					"ESS",
					"the equivalent sample size",
					true,
					new NumberValidator<Double>( new Double( 0 ), new Double( Double.MAX_VALUE ) ) ) );
			parameters.add( new SimpleParameter( DataType.STRING,
					"description",
					"a textual description or comment for the model",
					false,
					"none" ) );
			}catch(ParameterException doesnothappen){ }
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.SequenceScoringParameterSet#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		if( super.hasDefaultOrIsSet() ) {
			// return getAlphabet().isDiscrete();
			return ( (AlphabetContainerParameterSet)alphabet.getValue() ).isDiscrete();
		} else {
			return false;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return getInstanceClass().getSimpleName();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.SequenceScoringParameterSet#clone()
	 */
	@Override
	public DGMParameterSet clone() throws CloneNotSupportedException {
		return (DGMParameterSet)super.clone();
	}

	/**
	 * This method can be used to set the ess (<b>e</b>quivalent <b>s</b>ample
	 * <b>s</b>ize) of this parameter set.
	 * 
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * 
	 * @throws IllegalValueException
	 *             if the ess is negative
	 */
	public void setEss( double ess ) throws IllegalValueException {
		getParameterForName( "ESS" ).setValue( ess );
	}
}
