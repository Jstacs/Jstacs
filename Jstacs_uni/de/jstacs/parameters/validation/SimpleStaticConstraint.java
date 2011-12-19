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

/**
 * Class for a {@link Constraint} that checks values against static values using
 * the comparison operators defined in the interface {@link Constraint}. As long
 * as just an upper and lower bound shall be checked, {@link NumberValidator}
 * has the same functionality and is easier to use. Thus it is recommended to
 * use this class for more complex constraints, e.g. multiple allowed ranges or
 * combinations of static and reference-based constraints.
 * 
 * @author Jan Grau
 * 
 * @see de.jstacs.parameters.validation.NumberValidator
 */
public class SimpleStaticConstraint implements Constraint {

	private int comparisonOperator;
	private Object reference;
	private String errorMessage;

	/**
	 * Creates a new {@link SimpleStaticConstraint} from a {@link Number}
	 * -reference and a comparison operator as defined in {@link Constraint}.
	 * 
	 * @param reference
	 *            the {@link Number}-reference
	 * @param comparisonOperator
	 *            the comparison operator
	 */
	public SimpleStaticConstraint(Number reference, int comparisonOperator) {
		this.comparisonOperator = comparisonOperator;
		this.reference = reference;
	}

	/**
	 * Creates a new {@link SimpleStaticConstraint} from a {@link String}
	 * -reference and a comparison operator as defined in {@link Constraint}.
	 * 
	 * @param reference
	 *            the {@link String}-reference
	 * @param comparisonOperator
	 *            the comparison operator
	 */
	public SimpleStaticConstraint(String reference, int comparisonOperator) {
		this.comparisonOperator = comparisonOperator;
		this.reference = reference;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SimpleStaticConstraint} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if <code>representation</code> could not be parsed
	 */
	public SimpleStaticConstraint(StringBuffer representation)
			throws NonParsableException {
		fromXML(representation);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public SimpleStaticConstraint clone() throws CloneNotSupportedException {
		if (reference instanceof String) {
			SimpleStaticConstraint clone = new SimpleStaticConstraint(
					(String) reference, comparisonOperator);
			clone.errorMessage = errorMessage;
			return clone;
		} else if (reference instanceof Number) {
			SimpleStaticConstraint clone = new SimpleStaticConstraint(
					(Number) reference, comparisonOperator);
			clone.errorMessage = errorMessage;
			return clone;
		} else {
			throw new CloneNotSupportedException(
					"Reference was of non-expected type!");
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.validation.Constraint#check(java.lang.Object)
	 */
	@SuppressWarnings("unchecked")
	public boolean check(Object value) {
		if (value == null) {
			errorMessage = "Value was null";
			return false;
		}

		if (comparisonOperator == EQUALS) {
			if (!reference.equals(value)) {
				errorMessage = "Value must be equal to " + reference.toString();
				return false;
			}
		} else if (comparisonOperator == LT || comparisonOperator == GT
				|| comparisonOperator == LEQ || comparisonOperator == GEQ) {
			if (value instanceof Comparable) {
				if (comparisonOperator == LT
						&& !(((Comparable) value).compareTo(reference) < 0)) {
					errorMessage = "Value must be less than "
							+ reference.toString();
					return false;
				} else if (comparisonOperator == GT
						&& !(((Comparable) value).compareTo(reference) > 0)) {
					errorMessage = "Value must be greater than "
							+ reference.toString();
					return false;
				} else if (comparisonOperator == LEQ
						&& ((Comparable) value).compareTo(reference) > 0) {
					errorMessage = "Value must be less than or equal to "
							+ reference.toString();
					return false;
				} else if (comparisonOperator == GEQ
						&& ((Comparable) value).compareTo(reference) < 0) {
					errorMessage = "Value must be greater than or equal to "
							+ reference.toString();
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
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags( buf, reference.getClass().getName(),	"className" );
		XMLParser.appendObjectWithTags( buf, reference.toString(), "value");
		XMLParser.addTags( buf, "reference" );
		XMLParser.appendObjectWithTags( buf, comparisonOperator, "comparisonOperator" );
		XMLParser.appendObjectWithTags( buf, errorMessage, "errorMessage" );
		XMLParser.addTags( buf, "simpleStaticConstraint" );
		return buf;
	}

	/**
	 * Parses a {@link SimpleStaticConstraint} from the XML representation as
	 * returned by {@link SimpleStaticConstraint#toXML()}.
	 * 
	 * @param representation
	 *            the XML representation
	 * 
	 * @throws NonParsableException
	 *             the XML code could not be parsed
	 */
	public void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag(representation,
				"simpleStaticConstraint");
		StringBuffer ref = XMLParser.extractForTag(representation, "reference");
		try {
			Class refClass = Class.forName(XMLParser.extractObjectForTags(ref, "className", String.class));
			reference = refClass.getConstructor(new Class[] { String.class })
					.newInstance(
							new Object[] { XMLParser.extractObjectForTags(ref, "value", String.class) });
		} catch (Exception e) {
			throw new NonParsableException(e.getMessage());
		}
		comparisonOperator = XMLParser.extractObjectForTags(representation, "comparisonOperator", int.class );
		errorMessage = XMLParser.extractObjectForTags(representation, "errorMessage", String.class );
	}

	public String getErrorMessage() {
		return errorMessage;
	}

}
