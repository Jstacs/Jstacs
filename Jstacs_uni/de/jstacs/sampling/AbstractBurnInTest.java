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

package de.jstacs.sampling;

import de.jstacs.InstantiableFromParameterSet;
import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;

/**
 * This abstract class implements some of the methods of {@link BurnInTest} to
 * alleviate the implementation of efficient and new burn-in tests.
 * 
 * @author Jens Keilwagen
 */
public abstract class AbstractBurnInTest implements BurnInTest, InstantiableFromParameterSet {

	/**
	 * This array contains all values that will be set via
	 * {@link BurnInTest#setValue(double)}. The index of the array element
	 * denotes the sampling index.
	 */
	protected DoubleList[] values;

	/**
	 * This value indicates the current sampling index.
	 */
	private int currentSamplingIndex;

	/**
	 * This switch indicates whether the burn-in length has been computed.
	 */
	private boolean computed;

	/**
	 * The length of the burn-in phase.
	 */
	private int burnInLength;
	
	/**
	 * The set of parameters.
	 */
	private AbstractBurnInTestParameterSet parameters;

	/**
	 * This is the main constructor that creates a burn-in test given a specified
	 * set of parameters
	 * 
	 * @param parameters
	 *            set of parameters
	 *            
	 * @throws CloneNotSupportedException if the parameters can not be cloned            
	 */
	protected AbstractBurnInTest( AbstractBurnInTestParameterSet parameters ) throws CloneNotSupportedException {
		this.parameters = parameters.clone();
		int starts = this.parameters.getNumberOfStarts();
		values = new DoubleList[starts];
		for( int i = 0; i < starts; i++ ) {
			values[i] = new DoubleList();
		}
		resetAllValues();
	}
	
	/**
	 * This is the constructor for the {@link de.jstacs.Storable} interface.
	 * 
	 * @param rep
	 *            the XML representation
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} could not be parsed
	 */
	protected AbstractBurnInTest( StringBuffer rep ) throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag( rep, getXMLTag() );
		parameters = XMLParser.extractObjectForTags( xml, "parameters", AbstractBurnInTestParameterSet.class );
		computed = XMLParser.extractObjectForTags( xml, "computed", boolean.class );		burnInLength = XMLParser.extractObjectForTags( xml, "burnInLength", int.class );
		currentSamplingIndex = XMLParser.extractObjectForTags( xml, "currentSamplingIndex", int.class );
		values = XMLParser.extractObjectForTags( xml, "values", DoubleList[].class );
		setFurtherInformation( xml );
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	public AbstractBurnInTest clone() throws CloneNotSupportedException {
		return (AbstractBurnInTest)super.clone();
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.InstantiableFromParameterSet#getCurrentParameterSet()
	 */
	public AbstractBurnInTestParameterSet getCurrentParameterSet() throws CloneNotSupportedException {
		return parameters.clone();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.gibbssampling.BurnInTest#resetAllValues()
	 */
	public final void resetAllValues() {
		currentSamplingIndex = -1;
		for( int i = 0, end = values.length; i < end; i++ ) {
			values[i].clear();
		}
		computed = false;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.gibbssampling.BurnInTest#setCurrentSamplingIndex(int)
	 */
	public final void setCurrentSamplingIndex( int index ) {
		currentSamplingIndex = index;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.gibbssampling.BurnInTest#setValue(double)
	 */
	public final void setValue( double val ) {
		values[currentSamplingIndex].add( val );
		computed = false;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public final StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( values.length * 10000 );
		XMLParser.appendObjectWithTags( xml, parameters, "parameters" );
		XMLParser.appendObjectWithTags( xml, computed, "computed" );
		XMLParser.appendObjectWithTags( xml, burnInLength, "burnInLength" );
		XMLParser.appendObjectWithTags( xml, currentSamplingIndex, "currentSamplingIndex" );
		XMLParser.appendObjectWithTags( xml, values, "values" );
		xml.append( getFurtherInformation() );
		XMLParser.addTags( xml, getXMLTag() );
		return xml;
	}

	/**
	 * This method returns the XML tag that is used in
	 * {@link AbstractBurnInTest#toXML()} and
	 * {@link AbstractBurnInTest#AbstractBurnInTest(StringBuffer)}.
	 * 
	 * @return the XML tag that is used in {@link AbstractBurnInTest#toXML()}
	 *         and {@link AbstractBurnInTest#AbstractBurnInTest(StringBuffer)}
	 */
	protected abstract String getXMLTag();

	/**
	 * This method returns further information for the
	 * {@link AbstractBurnInTest}. It enables to store test specific data via
	 * the method {@link AbstractBurnInTest#toXML()}. <br>
	 * <br>
	 * 
	 * This method should only be used in {@link AbstractBurnInTest#toXML()}.
	 * 
	 * @return further information in XML format
	 */
	protected abstract StringBuffer getFurtherInformation();

	/**
	 * This method sets further information for the {@link AbstractBurnInTest}.
	 * It enables to load test specific data in the constructor
	 * {@link AbstractBurnInTest#AbstractBurnInTest(StringBuffer)}. <br>
	 * <br>
	 * 
	 * This method should only be used in
	 * {@link AbstractBurnInTest#AbstractBurnInTest(StringBuffer)}.
	 * 
	 * @param xml
	 *            contains further information in XML format
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} could not be parsed
	 */
	protected abstract void setFurtherInformation( StringBuffer xml ) throws NonParsableException;

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.gibbssampling.BurnInTest#getLengthOfBurnIn()
	 */
	public final int getLengthOfBurnIn() {
		if( !computed ) {
			burnInLength = computeLengthOfBurnIn();
			computed = true;
		}
		return burnInLength;
	}

	/**
	 * Computes and returns the length of the burn-in phase using the values
	 * from {@link BurnInTest#setValue(double)}. This method is used by
	 * {@link AbstractBurnInTest#getLengthOfBurnIn()}. The result is stored in
	 * an internal variable to avoid multiple meaningless assessments.
	 * 
	 * @return the length of the burn-in phase
	 */
	protected abstract int computeLengthOfBurnIn();
}
