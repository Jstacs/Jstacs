/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */
package de.jstacs.sequenceScores.differentiable.logistic;

import de.jstacs.NonParsableException;
import de.jstacs.data.Sequence;
import de.jstacs.io.XMLParser;

/**
 * 
 * @author Jens Keilwagen
 */
public class ProductConstraint implements LogisticConstraint {

	private int[] pos;
	
	public ProductConstraint( int... pos ) {
		//TODO
		this.pos = pos.clone();
	}

	public ProductConstraint( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, XML_TAG );
		pos = XMLParser.extractObjectForTags( xml, "pos", int[].class );
	}
	
	private static final String XML_TAG = ProductConstraint.class.getSimpleName();

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, pos, "pos" );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}
	
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	public ProductConstraint clone() throws CloneNotSupportedException {
		ProductConstraint clone = (ProductConstraint) super.clone();
		clone.pos = pos.clone();
		return clone;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.differentiable.logistic.LogisticConstraint#getValue(de.jstacs.data.Sequence, int)
	 */
	@Override
	public double getValue( Sequence seq, int start ) {
		double res = 1;
		for( int i = 0; i < pos.length; i++ ) {
			res *= seq.continuousVal( start + pos[i] );
		}
		return res;
	}
}
