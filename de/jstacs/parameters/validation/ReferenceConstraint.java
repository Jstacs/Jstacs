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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.parameters.validation;

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.Parameter;

/**
 * Class for a {@link Constraint} that defines a condition on one
 * {@link Parameter} (the one containing the {@link ConstraintValidator}) with
 * respect to another {@link Parameter}.
 * 
 * @author Jan Grau
 */
public abstract class ReferenceConstraint implements Constraint {

	/**
	 * The reference to the {@link Parameter} that is part of the condition.
	 */
	protected Parameter constraintParameter;
	/**
	 * The message of the last error or <code>null</code>.
	 */
	protected String errorMessage;

	/**
	 * Creates a new {@link ReferenceConstraint} with respect to the
	 * {@link Parameter} <code>constraintParameter</code>.
	 * 
	 * @param constraintParameter
	 *            the parameter
	 */
	public ReferenceConstraint(Parameter constraintParameter) {
		this.constraintParameter = constraintParameter;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ReferenceConstraint} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation
	 * 
	 * @throws NonParsableException
	 *             if <code>representation</code> could not be parsed
	 */
	public ReferenceConstraint(StringBuffer representation)
			throws NonParsableException {
		fromXML(representation);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public abstract ReferenceConstraint clone()
			throws CloneNotSupportedException;

	/**
	 * Sets the standard fields of <code>clone</code> to the values (or clones
	 * of them) of this instance.
	 * 
	 * @param clone
	 *            the clone
	 * 
	 * @throws CloneNotSupportedException
	 *             if the contained {@link Parameter} could not be cloned
	 */
	protected void fillWithStandardFieldForClone(ReferenceConstraint clone)
			throws CloneNotSupportedException {
		clone.constraintParameter = constraintParameter.clone();
		clone.errorMessage = errorMessage;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.validation.Constraint#check(java.lang.Object)
	 */
	public abstract boolean check(Object value);

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.validation.Constraint#getErrorMessage()
	 */
	public String getErrorMessage() {
		return errorMessage;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags(buf, constraintParameter, "parameter");
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		XMLParser.addTags(buf, "referenceConstraint");
		return buf;
	}

	/**
	 * Parses a {@link ReferenceConstraint} from the XML representation as
	 * returned by {@link ReferenceConstraint#toXML()}.
	 * 
	 * @param representation
	 *            the XML representation
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag(representation,
				"referenceConstraint");
		constraintParameter = XMLParser.extractObjectForTags( representation, "parameter", Parameter.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		errorMessage = XMLParser.extractObjectForTags(representation, "errorMessage", String.class );
	}

}
