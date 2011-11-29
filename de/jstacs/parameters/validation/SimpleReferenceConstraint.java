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
import de.jstacs.parameters.SimpleParameter;

/**
 * Class for a {@link ReferenceConstraint} that checks for &quot;simple&quot;
 * conditions as defined in the interface {@link Constraint}.
 * 
 * @author Jan Grau
 */
public class SimpleReferenceConstraint extends ReferenceConstraint {

	private int comparisonOperator;

	/**
	 * Creates a new {@link SimpleReferenceConstraint} from a reference
	 * {@link SimpleParameter} and a comparison operator, which is one of the
	 * values defined in the {@link Constraint} interface.
	 * 
	 * @param constraintParameter
	 *            the reference {@link SimpleParameter}
	 * @param comparisonOperator
	 *            the comparison operator
	 */
	public SimpleReferenceConstraint(SimpleParameter constraintParameter,
			int comparisonOperator) {
		super(constraintParameter);
		this.comparisonOperator = comparisonOperator;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SimpleReferenceConstraint} from its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if <code>representation</code> could not be parsed
	 */
	public SimpleReferenceConstraint(StringBuffer representation)
			throws NonParsableException {
		super(representation);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.validation.ReferenceConstraint#clone()
	 */
	@Override
	public SimpleReferenceConstraint clone() throws CloneNotSupportedException {
		SimpleReferenceConstraint clone = new SimpleReferenceConstraint(null,
				comparisonOperator);
		this.fillWithStandardFieldForClone(clone);
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.parameters.validation.ReferenceConstraint#check(java.lang.Object
	 * )
	 */
	@SuppressWarnings("unchecked")
	@Override
	public boolean check(Object value) {
		if (value == null) {
			errorMessage = "Value was null";
			return false;
		}
		if (!constraintParameter.hasDefaultOrIsSet()) {
			errorMessage = "Value of parameter "
					+ constraintParameter.getName() + " must be set first";
			return false;
		}

		Object constValue = constraintParameter.getValue();
		if (comparisonOperator == EQUALS) {
			if (!constraintParameter.getValue().equals(value)) {
				errorMessage = "Value must be equal to value of "
						+ constraintParameter.getName();
				return false;
			}
		} else if (comparisonOperator == LT || comparisonOperator == GT
				|| comparisonOperator == LEQ || comparisonOperator == GEQ) {
			if (value instanceof Comparable && constValue instanceof Comparable) {
				if (comparisonOperator == LT
						&& !(((Comparable) value).compareTo(constValue) < 0)) {
					errorMessage = "Value must be less than value of "
							+ constraintParameter.getName();
					return false;
				} else if (comparisonOperator == GT
						&& !(((Comparable) value).compareTo(constValue) > 0)) {
					errorMessage = "Value must be greater than value of "
							+ constraintParameter.getName();
					return false;
				} else if (comparisonOperator == LEQ
						&& ((Comparable) value).compareTo(constValue) > 0) {
					errorMessage = "Value must be less than or equal to value of "
							+ constraintParameter.getName();
					return false;
				} else if (comparisonOperator == GEQ
						&& ((Comparable) value).compareTo(constValue) < 0) {
					errorMessage = "Value must be greater than or equal to value of "
							+ constraintParameter.getName();
					return false;
				}
			} else {
				errorMessage = "Both values must be comparable";
				return false;
			}
		} else {
			errorMessage = "Wrong comparison operator";
			return false;
		}
		errorMessage = null;
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.validation.ReferenceConstraint#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags(buf, "superConstraint");
		XMLParser.appendObjectWithTags(buf, comparisonOperator,
				"comparisonOperator");
		XMLParser.addTags(buf, "simpleReferenceConstraint");
		return buf;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.parameters.validation.ReferenceConstraint#fromXML(java.lang
	 * .StringBuffer)
	 */
	@Override
	public void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag(representation,
				"simpleReferenceConstraint");
		super.fromXML(XMLParser
				.extractForTag(representation, "superConstraint"));
		comparisonOperator = XMLParser.extractObjectForTags(representation, "comparisonOperator", int.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
	}

}
