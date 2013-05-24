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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete;

import java.text.NumberFormat;

import de.jstacs.InstantiableFromParameterSet;
import de.jstacs.NotTrainedException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.AbstractTrainableStatisticalModel;

/**
 * This is the main class for all <b>d</b>iscrete <b>g</b>raphical <b>m</b>odels
 * (DGM).
 * 
 * @author Jens Keilwagen
 * 
 * @see DGTrainSMParameterSet
 */
public abstract class DiscreteGraphicalTrainSM extends AbstractTrainableStatisticalModel implements InstantiableFromParameterSet {

	private static final long serialVersionUID = 1L;

	/**
	 * The current parameter set of the model.
	 */
	protected DGTrainSMParameterSet params;

	/**
	 * Indicates whether the model is trained or not.
	 */
	protected boolean trained;

	/**
	 * The default constructor. Creates a new {@link DiscreteGraphicalTrainSM}
	 * from a given {@link DGTrainSMParameterSet}.
	 * 
	 * @param params
	 *            the given parameter set
	 * 
	 * @throws CloneNotSupportedException
	 *             if the parameter set could not be cloned
	 * @throws IllegalArgumentException
	 *             if the parameter set is not instantiated
	 * @throws NonParsableException
	 *             if the parameter set is not parsable
	 * 
	 * @see AbstractTrainableStatisticalModel#AbstractTrainableStatisticalModel(de.jstacs.data.AlphabetContainer, int)
	 */
	public DiscreteGraphicalTrainSM( DGTrainSMParameterSet params ) throws CloneNotSupportedException, IllegalArgumentException,
															NonParsableException {
		super( params.getAlphabetContainer(), params.getLength() );
		fromParameterSet( params );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link DiscreteGraphicalTrainSM} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link DiscreteGraphicalTrainSM} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see AbstractTrainableStatisticalModel#AbstractTrainableStatisticalModel(StringBuffer)
	 */
	public DiscreteGraphicalTrainSM( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#clone()
	 */
	@Override
	public DiscreteGraphicalTrainSM clone() throws CloneNotSupportedException {
		DiscreteGraphicalTrainSM clone = (DiscreteGraphicalTrainSM)super.clone();
		clone.params = params.clone();
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.InstantiableFromParameterSet#getCurrentParameterSet()
	 */
	public final DGTrainSMParameterSet getCurrentParameterSet() throws Exception {
		return params.clone();
	}

	/**
	 * Returns a short description of the model that was given by the user in
	 * the parameter set.
	 * 
	 * @return a short description of the model
	 */
	public final String getDescription() {
		return (String)params.getParameterAt( 1 ).getValue();
	}

	/**
	 * This method returns the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * that is used in this model.
	 * 
	 * @return the ess
	 */
	public final double getESS() {
		return (Double)params.getParameterAt( 0 ).getValue();
	}

	private void fromParameterSet( ParameterSet parameters ) throws CloneNotSupportedException,
			IllegalArgumentException,
			NonParsableException {
		DGTrainSMParameterSet p = (DGTrainSMParameterSet)parameters;
		if( !p.hasDefaultOrIsSet() ) {
			throw new IllegalArgumentException( "Parameters were not correct instantiated" );
		}
		set( p, false );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected final void fromXML( StringBuffer representation ) throws NonParsableException {
		StringBuffer erg = XMLParser.extractForTag( representation, getXMLTag() );
		try {
			set( XMLParser.extractObjectForTags( erg, "ParameterSet", DGTrainSMParameterSet.class ), XMLParser.extractObjectForTags( erg, "trained", boolean.class ) );
		} catch ( CloneNotSupportedException e ) {
			NonParsableException n = new NonParsableException( e.getMessage() );
			n.setStackTrace( e.getStackTrace() );
			throw n;
		}
		setFurtherModelInfos( erg );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#isTrained()
	 */
	public final boolean isInitialized() {
		return trained;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.SequenceScore#toString(java.text.NumberFormat)
	 */
	@Override
	public String toString( NumberFormat nf ) {
		return getDescription();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public final StringBuffer toXML() {
		StringBuffer erg = new StringBuffer();
		XMLParser.appendObjectWithTags( erg, params, "ParameterSet" );
		XMLParser.appendObjectWithTags( erg, trained, "trained" );
		StringBuffer help = getFurtherModelInfos();
		if( help != null ) {
			erg.append( help );
		}
		XMLParser.addTags( erg, getXMLTag() );
		return erg;
	}

	/**
	 * Checks some conditions on a {@link Sequence}. These are in general
	 * conditions on the {@link de.jstacs.data.AlphabetContainer} of a (sub)
	 * {@link Sequence} between <code>startpos</code> und <code>endpos</code>.
	 * 
	 * @param sequence
	 *            the {@link Sequence}
	 * @param startpos
	 *            the startposition
	 * @param endpos
	 *            the endposition
	 * 
	 * @throws NotTrainedException
	 *             if the model is not trained
	 * @throws IllegalArgumentException
	 *             if some constraints are not fulfilled
	 */
	protected void check( Sequence sequence, int startpos, int endpos ) throws NotTrainedException, IllegalArgumentException {
		if( !trained ) {
			throw new NotTrainedException();
		} else if( !alphabets.checkConsistency( sequence.getAlphabetContainer().getSubContainer( startpos, endpos - startpos + 1 ) ) ) {
			throw new IllegalArgumentException( "This sequence is not possible with the given alphabet." );
		} else if( startpos < 0 ) {
			throw new IllegalArgumentException( "This startposition is impossible. Try: 0 <= startposition" );
		}
	}

	/**
	 * Returns further model information as a {@link StringBuffer}.
	 * 
	 * @return further model information like parameters of the distribution
	 *         etc. in XML format
	 * 
	 * @see DiscreteGraphicalTrainSM#toXML()
	 */
	protected abstract StringBuffer getFurtherModelInfos();

	/**
	 * Returns the XML tag that is used for this model in
	 * {@link #fromXML(StringBuffer)} and {@link #toXML()}.
	 * 
	 * @return the XML tag that is used in {@link #fromXML(StringBuffer)} and
	 *         {@link #toXML()}
	 * 
	 * @see DiscreteGraphicalTrainSM#fromXML(StringBuffer)
	 * @see DiscreteGraphicalTrainSM#toXML()
	 */
	protected abstract String getXMLTag();

	/**
	 * This method replaces the internal model information with those from a
	 * {@link StringBuffer}.
	 * 
	 * @param xml
	 *            contains the model information like parameters of the
	 *            distribution etc. in XML format
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} could not be parsed
	 * 
	 * @see DiscreteGraphicalTrainSM#fromXML(StringBuffer)
	 */
	protected abstract void setFurtherModelInfos( StringBuffer xml ) throws NonParsableException;

	/**
	 * Sets the parameters as internal parameters and does some essential
	 * computations. Used in <code>fromParameterSet</code>-methods.
	 * 
	 * @param params
	 *            the new {@link ParameterSet}
	 * @param trained
	 *            indicates if the model is trained or not
	 * 
	 * @throws CloneNotSupportedException
	 *             if the parameter set could not be cloned
	 * @throws NonParsableException
	 *             if the parameters of the model could not be parsed
	 */
	protected void set( DGTrainSMParameterSet params, boolean trained ) throws CloneNotSupportedException, NonParsableException {
		if( !params.getInstanceClass().equals( this.getClass() ) ) {
			throw new RuntimeException( "The parameters set can not be used for this class." );
		}
		this.params = params.clone();
		alphabets = params.getAlphabetContainer();
		length = params.getLength();
		this.trained = trained;
	}
}
