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

import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * Class for a {@link ParameterValidator} that is based on {@link Constraint}s.
 * Each instance of a {@link ConstraintValidator} may contain a set of
 * constraints that are all fulfilled, if {@link #checkValue(Object)} returns
 * <code>true</code>.
 * 
 * @author Jan Grau
 */
public class ConstraintValidator implements ParameterValidator {

	private LinkedList<Constraint> constraints;
	private String errorMessage;

	/**
	 * Creates a new {@link ConstraintValidator} having an empty list of
	 * {@link Constraint}s, i.e. the constraints of this
	 * {@link ConstraintValidator} are always fulfilled before additional
	 * {@link Constraint}s are added using {@link #addConstraint(Constraint)}.
	 */
	public ConstraintValidator() {
		constraints = new LinkedList<Constraint>();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ConstraintValidator} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if <code>representation</code> could not be parsed
	 */
	public ConstraintValidator(StringBuffer representation)
			throws NonParsableException {
		fromXML(representation);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public ConstraintValidator clone() throws CloneNotSupportedException {
		ConstraintValidator clone = new ConstraintValidator();
		Iterator<Constraint> constIt = constraints.iterator();
		while (constIt.hasNext()) {
			clone.constraints.add(constIt.next().clone());
		}
		clone.errorMessage = errorMessage;
		return clone;
	}

	/**
	 * Adds an additional {@link Constraint} to the list of {@link Constraint}s.
	 * 
	 * @param c
	 *            the {@link Constraint} to be added
	 */
	public void addConstraint(Constraint c) {
		constraints.add(c);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.parameters.validation.ParameterValidator#checkValue(java.lang
	 * .Object)
	 */
	public boolean checkValue(Object value) {
		Iterator<Constraint> it = constraints.iterator();
		boolean checkPassed = true;
		errorMessage = null;
		while (it.hasNext()) {
			if (!it.next().check(value)) {
				if (errorMessage == null) {
					errorMessage = it.next().getErrorMessage();
				} else {
					errorMessage += "; " + it.next().getErrorMessage();
				}
				checkPassed = false;
			}
		}
		return checkPassed;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.validation.ParameterValidator#getErrorMessage()
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
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		XMLParser.appendObjectWithTags(buf, constraints.size(), "size");
		Iterator<Constraint> it = constraints.iterator();
		StringBuffer buf2 = new StringBuffer();
		while (it.hasNext()) {
			XMLParser.appendObjectWithTags(buf2, it.next(), "constraint");
		}
		XMLParser.addTags(buf2, "constraints");
		buf.append(buf2);
		XMLParser.addTags(buf, "referenceValidator");
		return buf;

	}

	/**
	 * Parses a {@link ConstraintValidator} from the XML representation as
	 * returned by {@link ConstraintValidator#toXML()}.
	 * 
	 * @param representation
	 *            the XML representation
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag( representation, "referenceValidator" );
		errorMessage = XMLParser.extractObjectForTags( representation, "errorMessage", String.class );
		int size = XMLParser.extractObjectForTags( representation, "size", int.class );
		constraints = new LinkedList<Constraint>();
		representation = XMLParser.extractForTag( representation, "constraints" );
		for (int i = 0; i < size; i++) {
			constraints.add(XMLParser.extractObjectForTags( representation, "constraint", Constraint.class ));
		}

	}

}
