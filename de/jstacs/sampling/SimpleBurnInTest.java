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

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * This is a very simple test for the length of the burn-in phase. It returns
 * for all samplings the length given in the constructor.
 * 
 * @author Jens Keilwagen
 * 
 * @deprecated
 */
public class SimpleBurnInTest implements BurnInTest {

	private int burnInLength;

	/**
	 * This is the main constructor that creates an instance of
	 * {@link SimpleBurnInTest} with fixed burn-in length.
	 * 
	 * @param burnInLength
	 *            the fixed length of the burn-in
	 */
	public SimpleBurnInTest( int burnInLength ) {
		if( burnInLength >= 0 ) {
			this.burnInLength = burnInLength;
		} else {
			throw new IllegalArgumentException( "The length of the burn-in phase has to be non-negative" );
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SimpleBurnInTest} out of its XML representation.
	 * 
	 * @param xml
	 *            the {@link StringBuffer} containing the model as XML
	 *            representation
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} can not be parsed
	 */
	public SimpleBurnInTest( StringBuffer xml ) throws NonParsableException {
		burnInLength = XMLParser.extractObjectForTags( XMLParser.extractForTag( xml, getClass().getSimpleName() ), "burnInLength", int.class );
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	public SimpleBurnInTest clone() throws CloneNotSupportedException {
		return (SimpleBurnInTest)super.clone();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.gibbssampling.BurnInTest#setCurrentSamplingIndex(int)
	 */
	public void setCurrentSamplingIndex( int index ) {}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.gibbssampling.BurnInTest#setValue(double)
	 */
	public void setValue( double val ) {}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.gibbssampling.BurnInTest#resetAllValues()
	 */
	public void resetAllValues() {}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.gibbssampling.BurnInTest#getLengthOfBurnIn()
	 */
	public int getLengthOfBurnIn() {
		return burnInLength;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 2000 );
		XMLParser.appendObjectWithTags( xml, burnInLength, "burnInLength" );
		XMLParser.addTags( xml, getClass().getSimpleName() );
		return xml;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.gibbssampling.BurnInTest#getInstanceName()
	 */
	public String getInstanceName() {
		return "simple burn in test (length=" + burnInLength + ")";
	}
}
